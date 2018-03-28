########################################################################
### Jae (Muschen Lab) Collaboration Cell Surface Proteomics on Jurkats/Jekos
########################################################################

############################################################
# Step 1 - Set working directory
############################################################
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Jae_surfaceomics/analysis/")



############################################################
# Step 2 - Read in protein groups CSV data file
############################################################
jeko.raw = read.delim("proteinGroups-jeko.txt", header = TRUE, stringsAsFactors = FALSE)
jurk.raw = read.delim("proteinGroups-jurkat.txt", header = TRUE, stringsAsFactors = FALSE)



############################################################
# Step 3 - Data clean-up
############################################################
cleanData = function(df) {
  # df = data frame containing raw data proteomics data
  
  ## Remove reverse proteins, contaminations, and only-identified-by-site
  df = df[df$Reverse != "+", ]
  df = df[df$Potential.contaminant != "+", ]
  df = df[df$Only.identified.by.site != "+", ]
  
  ## Log2-transform LFQ data
  LFQ.names = grep("^LFQ", names(df), value = T)
  log2.names = sub("^LFQ.intensity", "LOG2", LFQ.names)
  df[log2.names] = log2(df[LFQ.names])
  df[log2.names] = lapply(df[log2.names], function(x) {  # convert missing values to NA
    x[!is.finite(x)] = NA
    return(x)
  })

  ## Extract UniProtID from "Majority.protein.IDs"
  df$UniProtID = sapply(df$Majority.protein.IDs,   # Iterate through protein groups
                        function(IDs) {
                          IDs = strsplit(IDs, ";")[[1]]
                          regex = regexpr("(?<=\\|).*(?=\\|)", IDs, perl = TRUE)
                          out = rep(NA, length(IDs))
                          out[regex != -1] = regmatches(IDs, regex)
                          paste(out, collapse = ";")
                        })
  
  ## Extract gene symbols (NOT HGNC) from "Majority.protein.IDs"
  df$Protein = sapply(df$Majority.protein.IDs,   # Iterate through protein groups
                        function(IDs) {
                          IDs = strsplit(IDs, ";")[[1]]
                          regex = regexpr("(?<=\\|.{6}\\|).*(?=_)", IDs, perl = TRUE)
                          out = rep(NA, length(IDs))
                          out[regex != -1] = regmatches(IDs, regex)
                          paste(out, collapse = ";")
                        })
  
  ## Create Entry column from "UniProtID"
  df$Entry = sapply(df$UniProtID, function(IDs) strsplit(IDs, ";")[[1]][1])
  
  ## Extract gene names from first element of "Fasta.headers"
  df$Name = sapply(df$Fasta.headers,   # Iterate through protein groups
                   function(IDs) {
                     ID = strsplit(IDs, ";")[[1]][1]
                     regex = regexpr("(?<=HUMAN).*(?=OS)", ID, perl = TRUE)
                     out = NA
                     out[regex != -1] = regmatches(IDs, regex)
                     out
                   })
  
  ## Organize columns
  razorUniquePep.names = grep("^Razor...unique.peptides", names(df), value = TRUE)
  keep = c("Entry", "UniProtID", "Protein", "Name", "Number.of.proteins", log2.names, razorUniquePep.names)
  
  return(df[, keep])
}
jeko = cleanData(jeko.raw)
jurk = cleanData(jurk.raw)



############################################################
# Step 4 - Data Filtering
############################################################
filter_valids = function(df, conditions, min_count, at_least_one = FALSE,
                         advanced = FALSE, advanced_cutoff = c(1,1)) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  # advanced = logical indicating whether to perform an additional filter to remove cases 
  #     where the mean of one group with "few" observations > mean of another group with
  #     "min_count" observations). Designed for conditions = 2
  # advanced_cutoff = a numeric vector of length 2 dictating the threshold for "few" observation 
  #     in the two defined conditions
  require(dplyr)
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 columns
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
    
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  if (advanced) {
    DISCARD = logical(nrow(df))
    cond1 = as.matrix(df[cond.names[[1]]])   # Extract group 1
    cond2 = as.matrix(df[cond.names[[2]]])   # Extract group 2
    
    cond1[!is.finite(cond1)] = NA   # Convert invalid values to NA
    cond2[!is.finite(cond2)] = NA   # Convert invalid values to NA
    
    valid1 = rowSums(is.finite(cond1))   # Count number of valids in group 1
    valid2 = rowSums(is.finite(cond2))   # Count number of valids in group 2
    
    for (i in 1:nrow(df)) {   # Iterate through each row
      if ((valid1[i] > 0) && (valid1[i] <= advanced_cutoff[1]) && df$KEEP[i]) {  # Mean group 1 > 2?
        if ((mean(cond1[i, ], na.rm = TRUE) > mean(cond2[i, ], na.rm = TRUE))) DISCARD[i] = TRUE
      } else if ((valid2[i] > 0) && (valid2[i] <= advanced_cutoff[2]) && df$KEEP[i]) {  # Mean group 2 > 1?
        if ((mean(cond2[i, ], na.rm = TRUE) > mean(cond1[i, ], na.rm = TRUE))) DISCARD[i] = TRUE
      }
    }
    df$DISCARD = DISCARD
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}
jekoF = filter_valids(jeko, c("JEKO.KO", "JEKO.WT"), c(2, 2),
                      at_least_one = TRUE, advanced = TRUE, advanced_cutoff = c(1, 1))
jurkF = filter_valids(jurk, c("JURKAT.KO", "JURKAT.WT"), c(2, 2),
                      at_least_one = TRUE, advanced = TRUE, advanced_cutoff = c(1, 1))



############################################################
# Step 5 - Data Normalization
############################################################
hybrid_impute = function(df, conditions, use_keep = TRUE, use_discard = TRUE) {
  # df = data frame containing filtered 
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to imputation (WILL AFFECT IMPUTATION OUTCOME)
  # use_discard = filter out rows using DISCARD column prior to imputation (WILL AFFECT IMPUTATION OUTCOME)
  require(imputeLCMD)
  require(dplyr)
  
  # Apply KEEP filter
  if (use_keep) df = df[df$KEEP, ]
  if (use_discard) df = df[!df$DISCARD, ]
  
  # Group column names by condition
  log2DF = select(df, starts_with("LOG2"))   # Extract LOG2 dataframe
  DF = select(df, -starts_with("LOG2"))
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, names(log2DF), value = TRUE, perl = TRUE))
  
  # Create new columns indicating whether the values are imputed
  impute.names = sub("^LOG2", "IMP", unlist(cond.names))
  DF[impute.names] = lapply(unlist(cond.names), function(x) !is.finite(log2DF[[x]]))
  
  # Imputation
  set.seed(1)
  imputeDF = lapply(cond.names,   # Impute each group separately
                    function(x) {
                      tempDF = as.matrix(log2DF[x])
                      modelS = model.Selector(tempDF)
                      impute.MAR.MNAR(tempDF, modelS, method.MAR = "MLE", method.MNAR = "MinProb")
                    })
  imputeDF = do.call(cbind, imputeDF)   # Combine a list of data frames into a single data frame
  
  return(cbind(DF, imputeDF))
}
jekoFI = hybrid_impute(jekoF, c("JEKO.KO", "JEKO.WT"))
jurkFI = hybrid_impute(jurkF, c("JURKAT.KO", "JURKAT.WT"))



############################################################
# Step 6 - Data Analysis
############################################################
## Annotate data with surface protein information
membraneDF = read.csv("membrane-proteins-annotation.csv", stringsAsFactors = FALSE)
jekoFIA = left_join(jekoFI, membraneDF, by = "Entry")
jurkFIA = left_join(jurkFI, membraneDF, by = "Entry")


## Filter for membrane proteins
jekoFIA = filter(jekoFIA, membrane == "+")
jurkFIA = filter(jurkFIA, membrane == "+")


## Perform Welch's t-test (copied from MaxQuant_report_Christine.Rmd)
welch_test = function(df, diff = c("JEKO.KO", "JEKO.WT"), bioRep = c("1$", "2$", "3$")) {
  # Function only works for pairwise comparisons
  # df = data frame containing log2 intensity data (AFTER EXCLUDING OR IMPUTING MISSING VALUES)
  # diff = find the Ratio by computing element1 - element2
  # bioRep = biological replicate identifiers
  log2.names = grep("^LOG2", names(df), value = TRUE)   #Extract all LOG2. names
  group1 = grep(diff[[1]], log2.names, value = TRUE)
  group2 = grep(diff[[2]], log2.names, value = TRUE)
  group1 = sapply(bioRep, function(x) grep(x, group1, value = TRUE))
  group2 = sapply(bioRep, function(x) grep(x, group2, value = TRUE))
  
  # Perform Welch's two-sided paired two-sample t-test
  result = data.frame(Pvalue = numeric(0), MEAN.RATIO = numeric(0))
  for (i in 1:nrow(df)) {
    stats = t.test(unlist(df[i, group1]), unlist(df[i, group2]), alternative = "two.sided",
                   paired = TRUE, var.equal = FALSE)
    result[i, ] = c(stats$p.value, stats$estimate)
  }
  result$LOG.Pvalue = -log10(result$Pvalue)
  
  return(cbind(df, result))
}
jekoFIAT = welch_test(jekoFIA, diff = c("JEKO.KO", "JEKO.WT"))
jurkFIAT = welch_test(jurkFIA, diff = c("JURKAT.KO", "JURKAT.WT"))


## FDR correction
pValue_correction = function(df, adj.method = "BH") {
  # df = data frame containing PValue column (from welch_test function)
  # adj.method = character of length 1 indicating the method (see p.adjust.methods)
  df$adj.Pvalue = p.adjust(df$Pvalue, method = adj.method)
  df$LOG.adj.Pvalue = -log10(df$adj.Pvalue)
  return(df)
}
jekoFIAT = pValue_correction(jekoFIAT)
jurkFIAT = pValue_correction(jurkFIAT)



############################################################
# Step 7 - Data Visualization
############################################################
## Pairwise scatter plot
plot_pairs = function(df, pattern, use_keep = FALSE) {
  # df = data frame carrying data of interest
  # pattern = regex pattern to select column of interest
  # use_keep = TRUE means to construct plot on filtered values; FALSE uses all available values
  require(gpairs)
  require(scales)
  if (use_keep) {
    plot.df = df[df$KEEP, grep(pattern, names(df), value = TRUE, perl = TRUE)]
  } else {
    plot.df = df[grep(pattern, names(df), value = TRUE, perl = TRUE)]
  }
  # AUTOMATICALLY EXCLUDE -Inf VALUES
  plot.df = as.matrix(plot.df)
  plot.df[!is.finite(plot.df)] = NA
  
  gpairs(plot.df,
         upper.pars = list(scatter = "lm"),
         scatter.pars = list(pch = 20,
                             col = alpha("black", 0.3)),
         lower.pars = list(scatter = "stats"),
         stat.pars = list(verbose = FALSE, fontsize = 15))
}
#plot_pairs(jekoFIAT, "^LOG2")
#plot_pairs(jurkFIAT, "^LOG2")


## Histogram on imputed data
hist_impute = function(df, use_keep = TRUE) {
  # df = data frame containing imputed data
  # use_keep = logical indicating whether to use KEEP column to filter rows
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(stringr)
  
  if (use_keep) df = filter(df, KEEP)
  
  LOG2.df = dplyr::select(df, starts_with("LOG2"))
  impute.df = dplyr::select(df, starts_with("IMP"))
  
  # Reshape data into key-value pairs
  LOG2.df = gather(LOG2.df, "sample", "intensity")
  impute.df = gather(impute.df, "sample", "impute")
  
  # Combine data
  combine.df = bind_cols(LOG2.df, impute.df["impute"])
  
  # Create labels
  combine.df = mutate(combine.df, sample = sub("^LOG2\\.", "", sample)) %>%
    mutate(replicate = str_extract(sample, ".$")) %>%
    mutate(sample = sub(".$", "", sample))
  
  ggplot(combine.df, aes(x = intensity, fill = impute)) +
    geom_histogram(alpha = 0.3, binwidth = 0.4, position = "identity") +
    labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
    facet_grid(replicate ~ sample) +
    scale_fill_discrete(name = "Imputed",
                        breaks = c("FALSE", "TRUE"),
                        labels = c("-", "+"))
}
#hist_impute(jurkFIA, use_keep = TRUE)
#hist_impute(jekoFIA, use_keep = TRUE)


## Volcano plot
plot_volcano = function(df, ratio_cutoff = 1, P_cutoff = 0.05, use_adjusted = FALSE, labeling = TRUE) {
  # df = dataframe containing plot data
  # ratio_cutoff = threshold of biological significance
  # P_cutoff = threshold of statistical significance
  # use_adjusted = boolean indicating whether to use Pvalue or FDR-adjusted Pvalue (from pValue_correction function)
  # labeling = boolean indicating whether to include point labels
  if (use_adjusted) {
    df$Pvalue = df$adj.Pvalue
    df$LOG.Pvalue = df$LOG.adj.Pvalue
  }
  
  df$color = "black"
  df$color[df$MEAN.RATIO > ratio_cutoff & df$Pvalue < P_cutoff] = "darkgreen"
  df$color[df$MEAN.RATIO < -ratio_cutoff & df$Pvalue < P_cutoff] = "red"
  df$trans = 0.1   # point transparency
  df$trans[df$color != "black"] = 0.7
  df$color = factor(df$color)
  
  require(ggplot2)
  require(ggrepel)
  fig = ggplot(df) +
          geom_point(aes(MEAN.RATIO, LOG.Pvalue, colour = color), alpha = df$trans) +
          scale_colour_manual(values = levels(df$color)) + 
          geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "gray50") +
          geom_vline(xintercept = c(ratio_cutoff, -ratio_cutoff), linetype = "dashed", color = "gray50") +
          scale_x_continuous(
            name = expression("log"[2]*"( fold change over control )"),
            breaks = seq(-10, 10, 1)) +
          theme(legend.position="none")
          #annotate("text", x = 1.5, y = -log10(P_cutoff)-0.1, label = paste0("P-value = ", P_cutoff), 
          #         size = 3, colour = "gray50")
  
  if (labeling) {
    fig = fig + geom_text_repel(
      data = subset(df, df$color != "black"),
      aes(MEAN.RATIO, LOG.Pvalue, label = Protein, colour = color),
      size = 3)
  }
  
  if (use_adjusted) {
    print(fig + ylab(expression("- log"[10]*"( adjusted P-value )")))
  } else {
    print(fig + ylab(expression("- log"[10]*"( P-value )")))
  }
}
#plot_volcano(jurkFIAT, use_adjusted = FALSE)
#plot_volcano(jekoFIAT, use_adjusted = FALSE)



############################################################
# Step 8 - Data Table Organization
############################################################
organize_data = function(df) {
  # df = data frame containing final data for reporting (tailored for Jae_surfaceome_analysis.R)
  select(df, c(UniProtID, Protein, Name, MEAN.RATIO, Pvalue, adj.Pvalue, starts_with("LOG2"), 
               starts_with("IMP"), Number.of.proteins, starts_with("Razor"), Passes, singlepass, 
               multipass, GPI.anchored))
}
jekoFIATO = organize_data(jekoFIAT)
jurkFIATO = organize_data(jurkFIAT)


organize_data_table = function(df) {
  # df = data frame containing data processed by "organize_data" function
  select(df, c(UniProtID, Protein, Name, Razor...unique.peptides, MEAN.RATIO, Pvalue, adj.Pvalue)) %>%
    rename(Peptides = Razor...unique.peptides, log2FC = MEAN.RATIO) %>%
    mutate(log2FC = round(log2FC, 3), Pvalue = round(Pvalue, 3), adj.Pvalue = round(adj.Pvalue, 3))
}
jekoFIATO2 = organize_data_table(jekoFIATO)
jurkFIATO2 = organize_data_table(jurkFIATO)

# Save workspace
#save.image("Jae_workspace.RData")

