### ShinyApp accompanying Jae_surfaceome_analysis.R

#############################################
# Step 1 - Open workspace
#############################################
#setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Jae_surfaceomics/analysis/IFITM3")
#load("./Data/Jae_workspace.RData")   # Data analysis performed in Jae_surfaceome_analysis.R (issues when deploying app)



############################################################
# Step 2 - Read in protein groups CSV data file
############################################################
jeko.raw = read.delim("./Data/proteinGroups-jeko.txt", header = TRUE, stringsAsFactors = FALSE)
jurk.raw = read.delim("./Data/proteinGroups-jurkat.txt", header = TRUE, stringsAsFactors = FALSE)



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
membraneDF = read.csv("./Data/membrane-proteins-annotation.csv", stringsAsFactors = FALSE)
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


## Produce summary table
summary_stat = function(filteredDF, imputedDF, annotateDF) {
  # filteredDF = data frame after performing missing value filter
  # imputedDF = data frame after performing missing imputation
  # annotateDF = data frame after annotating for membrane proteins
  
  require(dplyr)
  df = select(filteredDF, starts_with("LOG2"))
  df2 = select(imputedDF, starts_with("IMP"))
  df3 = select(annotateDF, starts_with("IMP"))
  
  data.frame(
    Sample = sub("^LOG2.", "", names(df)),
    Raw = sapply(df, function(x) {
      paste0(sum(is.finite(x)), "/", length(x), " (", round(sum(is.finite(x))/length(x)*100, 1), "%)")
    }),
    Filtered = sapply(df2, function(x) {
      paste0(sum(!x), "/", length(x), " (", round(sum(!x)/length(x)*100, 1), "%)")
    }),
    Membrane = rep(
      paste0(nrow(df3), "/", nrow(df2), " (", round(nrow(df3)/nrow(df2)*100, 1), "%)"), ncol(df)
    ),
    Imputed = sapply(df3, function(x) {
      paste0(sum(x), "/", length(x), " (", round(sum(x)/length(x)*100, 1), "%)")
    })
  )
}


## Produce Popup Barplot
barplot_popup = function(df, index) {
  # df = data frame containing LOG2 and IMP columns
  # index = numeric vector of length 1 indicating the row to plot
  require(tidyr)
  require(plyr)   # Adjusting factor levels
  
  LOG2 = select(df, starts_with("LOG2")) %>%   # Extract LOG2 columns
    slice(index) %>%   # Extract specific row corresponding to selection
    gather("sample", "intensity") %>%   # Organize into columns
    mutate(sample = sub("LOG2.", "", sample))
  
  IMP = select(df, starts_with("IMP")) %>%   # Extract LOG2 columns
    slice(index) %>%   # Extract specific row corresponding to selection
    gather("sample", "Imputed") %>%   # Organize into columns
    mutate(sample = sub("IMP.", "", sample))
  
  DF = full_join(LOG2, IMP, by = "sample") %>%
    mutate(Replicates = str_extract(sample, ".$")) %>%
    mutate(sample = sub(".$", "", sample)) %>%
    mutate(Imputed = factor(Imputed, levels = c(F, T))) %>%   # Reordering factor levels
    mutate(Imputed = revalue(Imputed, c("TRUE" = "+", "FALSE" = "-")))   # Switching level labels
  
  ggplot(DF, aes(x = Replicates, y = intensity, fill = Imputed)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    scale_y_continuous(name = expression("log"[2]*"( intensity )"), 
                       limits = c(15, NA),
                       oob = rescale_none) +
    facet_grid(. ~ sample) +
    ggtitle(paste0("LOG2 Fold Change = ", round(df$MEAN.RATIO[index], 2)), 
            subtitle = paste0("P-value: ", round(df$Pvalue[index], 3), "    ", 
                              "Adj. P-value:", round(df$adj.Pvalue[index], 3)))
}



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
    dplyr::rename(Peptides = Razor...unique.peptides, log2FC = MEAN.RATIO) %>%
    mutate(log2FC = round(log2FC, 3), Pvalue = round(Pvalue, 3), adj.Pvalue = round(adj.Pvalue, 3))
}
jekoFIATO2 = organize_data_table(jekoFIATO)
jurkFIATO2 = organize_data_table(jurkFIATO)



#############################################
# Step 9 - ShinyApp code
#############################################
library(shinythemes)
library(shiny)
library(DT)  # For interactive data tables

# USER INTERFACEs ------------------------------------------
ui = navbarPage(
  "Cell Surface Proteomics",
  #theme = shinytheme("sandstone"),
  theme = shinytheme("simplex"),
  
  # ABOUT PAGE ----------------------------------------------
  tabPanel("About",
           h2("Surface Proteome Analysis on IFITM3-knockouts in JeKo-1 and Jurkat Cells"),
           hr(),
           fluidRow(
             column(8,
               h3("Motivation"),
               p(HTML("The aim of this experiment is to identify membrane proteins whose abundances are 
                      affected as a result of IFITM3 (Interferon-induced transmembrane protein 3) 
                      deletion. IFITM3 is a single-pass transmembrane protein localized in the plasma 
                      membrane and endosomes. It has been shown to play a role in maintaining 
                      cholesterol distribution in the cell membrane for the formation of lipid rafts, 
                      which in turn serve as anchoring platforms for many membrane proteins. Using 
                      the biotin-hydrazide-based cell surfaceome method, we hope to characterize surface
                      receptors that are lost when IFITM3 is knocked out.")),
               br(),
               h3("Experimental Procedure"),
               p(HTML("Two cell lines were used in this experiment: <b>JeKo-1</b> (B cell lymphoma) and
                      <b>Jurkat</b> (T cell lymphoma). IFITM3-knockout and control lines were 
                      generated from each cell type, resulting a total of four lines: JeKo-1 Cas9
                      control <b>(JEKO.WT)</b>, JeKo-1 Cas9 IFITM3 -/- <b>(JEKO.KO)</b>, Jurkat Cas9
                      control <b>(JURKAT.WT)</b>, and Jurkat Cas9 IFITM3 -/- <b>(JURKAT.KO)</b>.")),
               p(HTML("<ol><li><b>Cell Culture:</b>
                      Cells were received on ice and allowed to recover in culture for a week. After 
                      recovery, cells were harvested three to four days apart for biological
                      replication. Three total replicates were obtained for each line.</li>
                      <li><b>Surface Proteomics:</b>                   
                      During harvest, cells were washed in cold PBS to remove residual media prior to
                      surface proteins labeling with biocytin hydrazide. Samples were snap frozen
                      and stored in -80 degree C until all replicates were collected. Next, the cells
                      were lysed and the biotinylated proteins captured on NeutrAvidin beads. The beads 
                      were washed extensively to remove non-specific binders. To facilitate 
                      protein digestion, reduction and alkylation steps were performed before the 
                      overnight incubation with trypsin-LysC. After 16 hours, the peptides were separated
                      from the beads and desalted by reversed-phase chromatography on a C18 column. 
                      Finally, peptides were dried down by SpeedVac and stored in -80 degree C.</li>
                      <li><b>Liquid Chromatography-Mass Spectrometry (LC-MS) Analysis:</b>
                      Dried peptides were reconstituted in 2% acetonitrile/0.1% formic acid. The 
                      concentration was adjusted to 0.2 ug/uL for a 5 uL injection (1 ug) onto the 
                      LC-MS. Each sample was analyzed on a 4-hour graident.</li></ol>")),
               br(),
               h3("Data Analysis"),
               p(HTML("The acquired raw data were analyzed by the proteomics software MaxQuant to
                      match the identified spectra to peptides and to derive quantitative measures
                      of protein abundance. A search was performed against the UniProt/SwissProt database 
                      containing ~20,000 protein sequence. Match-between-run and LFQ settings were enabled. 
                      For more information on MaxQuant, please see the "),
                 a("online documentation", href = "https://www.nature.com/articles/nprot.2016.136", target="_blank"),
                 HTML(".")),
               p(HTML("The MaxQuant output was processed in a series of steps in preparation for 
                      visualization."),
                 HTML("<ol>
                        <li><b>Cleaning</b>: Contaminants and reverse proteins were excluded. 
                            Verbose columns were removed.</li>
                        <li><b>Transformation</b>: log<sub>2</sub>-transformation was applied to 
                            intensity data to normalize the distribution.</li>
                        <li><b>Filtering</b>: Proteins with less than 2 quantified observations in both WT and
                            KO groups were omitted from the analysis. Additionally, proteins carrying a single 
                            observation in one group that is greater than the mean of two or more 
                            observations in the other group were removed.</li>
                        <li><b>Imputation</b>: A hybrid imputation method was employed to assign substitute
                            values to missing data.</li>
                        <li><b>Annotation</b>: External annotation was used to classify detected proteins as
                            membrane-bound.</li>
                        <li><b>Statistical testing</b>: Welch's two-sample paired t-test was performed 
                            to evaluate the statistical significance of change in surface marker abundance
                            between WT and KO.</li>
                        <li><b>Visualization</b>: Dynamic visualization was implemented using the ShinyApp
                            package in R.</li>
                      </ol>")),
                 hr(),
                 helpText("Please direct any technical issues to Tony Lin (lin.yu.hsiu@ucsf.edu).")
               )
             )
           ),
  
  
  
  
  # JEKO PAGE ----------------------------------------------
  tabPanel("JEKO",
           fluidRow(
             column(6,
                    wellPanel(
                      h4("Summary"),
                      p(HTML("Summary statistics are displayed below for each sample. <b>Raw</b> shows the 
                        ratio of the number of quantified proteins to those identified. <b>Filtered</b> 
                        shows the ratios after applying a missing value filter. <b>Membrane</b> represents
                        the fraction of surface proteins in the filtered list. <b>Imputed</b> indicates 
                        the percentage of membrane proteins with imputed values.")),
                      tableOutput("summary_jeko")),
                      
                      wellPanel(
                        h4("Volcano Plot"),
                        p(HTML("Each point is a membrane protein plotted along the x- and y-axes, which 
                               represent the scales of biological and statistical significance, respectively.
                               The x-axis shows the log<sub>2</sub> fold change of KO over WT. Red indicates 
                               a statistically significant decrease in proteins on the membrane of the
                               IFITM3-depleted line, while green suggests an increase. Option to switch 
                               to the P-value adjusted according to the Benjamini-Hochberg procedure
                               to account for multiple testing.")),
                        fluidRow(
                          column(9, 
                                 sliderInput("volcano_jeko_ratioCutoff", label = "Fold Change Cutoff",
                                             min = 0, max = 3, step = 0.05, value = 1),
                                 sliderInput("volcano_jeko_pCutoff", label = "P-value Cutoff",
                                             min = 0.001, max = 1, step = 0.001, value = 0.05)
                          ),
                          column(3,
                                 radioButtons("volcano_jeko_pvalue", label = "Display Type",
                                              choices = list("P-value" = FALSE, "Adjusted P-value" = TRUE), 
                                              selected = FALSE),
                                 radioButtons("volcano_jeko_labels", label = "Display Labels",
                                              choices = list("On" = TRUE, "Off" = FALSE), 
                                              selected = TRUE)
                          )),
                        plotOutput(("volcano_jeko"))
               )
             ),
             
             column(6,
                    wellPanel(
                      h4("Pairwise Scatter Plot"),
                      p(HTML("The histograms display the distribution of
                             log<sub>2</sub>-transformed membrane protein intensities along the diagonal.
                             In addition, scatter plots between any two 
                             samples are shown at the upper right corner (lines of best fit 
                             in red). The corresponding Pearson's correlations, r, are found at the lower left 
                             corner. Proteins quantified in one sample but not the other are counted as \"missing\".")),
                      fluidRow(
                        column(6,
                          selectInput('pairwise_jeko_type', 'Samples:',
                                      c("WT", "KO", "WT + KO"), selectize = FALSE)),
                        column(6, 
                          selectInput('pairwise_jeko_display', 'Display Type:',
                                      c("Before imputation", "After imputation"), selectize = FALSE))
                      ),
                      plotOutput(("pairwise_jeko"))
                    ),
                    
                    wellPanel(
                      h4("Histogram of Imputed Values"),
                      p(HTML("Distribution of log<sub>2</sub>-transformed intensity values for membrane
                             proteins is depicted below for each sample. Replicates are in
                             rows and the groups (WT and KO) in columns. The histogram of imputed values
                            are and quantified measurements are superimposed.")),
                      plotOutput(("hist_jeko"))
                      )
             )

           ),
           wellPanel(
             h4("Data Table"),
             p(HTML("Measures of fold change, P-value, FDR-adjusted P-value, and the number of peptides
                    used for quantification are found below. The data table has been processed and 
                    contains only membrane proteins. Rows can be selected to reveal the measured
                    or imputed values used for statistical testing. Complete data tables are
                    available for download.")),
             downloadButton("DL_jeko_processed", "Download Processed CSV"),
             downloadButton("DL_jeko_raw", "Download MaxQuant Output"),
             hr(),
             DT::dataTableOutput('dataTable_jeko'))
           ),
  

  # JURK PAGE ----------------------------------------------
  tabPanel("JURKAT",
           fluidRow(
             column(6,
                    wellPanel(
                      h4("Summary"),
                      p(HTML("Summary statistics are displayed below for each sample. <b>Raw</b> shows the 
                             ratio of the number of quantified proteins to those identified. <b>Filtered</b> 
                             shows the ratios after applying a missing value filter. <b>Membrane</b> represents
                             the fraction of surface proteins in the filtered list. <b>Imputed</b> indicates 
                             the percentage of membrane proteins with imputed values.")),
                      tableOutput("summary_jurk")
                      ),
                    
                    wellPanel(
                     h4("Volcano Plot"),
                     p(HTML("Each point is a membrane protein plotted along the x- and y-axes, which 
                            represent the scales of biological and statistical significance, respectively.
                            The x-axis shows the log<sub>2</sub> fold change of KO over WT. Red indicates 
                            a statistically significant decrease in proteins on the membrane of the
                            IFITM3-depleted line, while green suggests an increase. Option to switch 
                            to the P-value adjusted according to the Benjamini-Hochberg procedure
                            to account for multiple testing.")),
                     fluidRow(
                       column(9, 
                              sliderInput("volcano_jurk_ratioCutoff", label = "Fold Change Cutoff",
                                          min = 0, max = 3, step = 0.05, value = 1),
                              sliderInput("volcano_jurk_pCutoff", label = "P-value Cutoff",
                                          min = 0.001, max = 1, step = 0.001, value = 0.05)
                       ),
                       column(3,
                              radioButtons("volcano_jurk_pvalue", label = "Display Type",
                                           choices = list("P-value" = FALSE, "Adjusted P-value" = TRUE), 
                                           selected = FALSE),
                              radioButtons("volcano_jurk_labels", label = "Display Labels",
                                           choices = list("On" = TRUE, "Off" = FALSE), 
                                           selected = TRUE)
                        )),
                     plotOutput(("volcano_jurk"))
                     )
              ),
             
             column(6,
                    wellPanel(
                      h4("Pairwise Scatter Plot"),
                      p(HTML("The histograms display the distribution of
                             log<sub>2</sub>-transformed membrane protein intensities along the diagonal.
                             In addition, scatter plots between any two 
                             samples are shown at the upper right corner (lines of best fit 
                             in red). The corresponding Pearson's correlations, r, are found at the lower left 
                             corner. Proteins quantified in one sample but not the other are counted as \"missing\".")),
                      fluidRow(
                        column(6,
                               selectInput('pairwise_jurk_type', 'Samples:',
                                           c("WT", "KO", "WT + KO"), selectize = FALSE)),
                        column(6, 
                               selectInput('pairwise_jurk_display', 'Display Type:',
                                           c("Before imputation", "After imputation"), selectize = FALSE))
                      ),
                      plotOutput(("pairwise_jurk"))
                    ),
                    
                    wellPanel(
                     h4("Histogram of Imputed Values"),
                     p(HTML("Distribution of log<sub>2</sub>-transformed intensity values for membrane
                            proteins is depicted below for each sample. Replicates are in
                            rows and the groups (WT and KO) in columns. The histogram of imputed values
                            are and quantified measurements are superimposed.")),
                     plotOutput(("hist_jurk"))
                     )
                             )
             ),
           wellPanel(
             h4("Data Table"),
             p(HTML("Measures of fold change, P-value, FDR-adjusted P-value, and the number of peptides
                    used for quantification are found below. The data table has been processed and 
                    contains only membrane proteins. Rows can be selected to reveal the measured
                    or imputed values used for statistical testing. Complete data tables are
                    available for download.")),
             downloadButton("DL_jurk_processed", "Download Processed CSV"),
             downloadButton("DL_jurk_raw", "Download MaxQuant Output"),
             hr(),
             DT::dataTableOutput('dataTable_jurk'))
             )
)



# SERVER ------------------------------------------------------------------
server = function(input, output, session) {
  
  # JEKO PAGE --------------------------------------------------
  # JEKO: Output summary table statistics
  output$summary_jeko <- renderTable({ summary_stat(jekoF, jekoFI, jekoFIA) })
  
  # JEKO: Output pairwise plot
  output$pairwise_jeko <- renderPlot({ 
    if (input$pairwise_jeko_type == "WT + KO") {
      pattern = "^LOG2"
    } else if (input$pairwise_jeko_type == "WT") {
      pattern = "^LOG2.JEKO.WT"
    } else (
      pattern = "^LOG2.JEKO.KO"
    )
    if (input$pairwise_jeko_display == "Before imputation") {   # Restore NA to imputed data set
      df = select(jekoFIA, starts_with("LOG2"))
      for (i in 1:ncol(df)) df[jekoFIA[[sub("LOG2", "IMP", names(df)[i])]], i] = NA
      plot_pairs(df, pattern)
    } else plot_pairs(jekoFIA, pattern)
  })
  
  # JEKO: Histogram of Imputed Value
  output$hist_jeko <- renderPlot({ hist_impute(jekoFIA, use_keep = TRUE) })
  
  # JEKO: Volcano plot
  output$volcano_jeko <- renderPlot({ plot_volcano(jekoFIAT, 
                                                   ratio_cutoff = input$volcano_jeko_ratioCutoff, 
                                                   P_cutoff = input$volcano_jeko_pCutoff, 
                                                   use_adjusted = input$volcano_jeko_pvalue,
                                                   labeling = input$volcano_jeko_labels) })
  
  # JEKO: Download Processed Data Table
  output$DL_jeko_processed <- downloadHandler(
    filename = function() {"processed_JEKO.csv"},
    content = function(file) {
      write.csv(jekoFIATO, file, row.names = FALSE)
    }
  )
  
  # JEKO: Download Raw Data Table
  output$DL_jeko_raw <- downloadHandler(
    filename = function() {"MaxQuant_proteinGroups_JEKO.csv"},
    content = function(file) {
      write.csv(jeko.raw, file, row.names = FALSE)
    }
  )
  
  # JEKO: Data Table
  output$dataTable_jeko <- DT::renderDataTable(
    DT::datatable(
      jekoFIATO2, 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # JEKO: Filtered Table Line Plot (in modal dialog)
  output$barPlot_jeko = renderPlot({
    barplot_popup(jekoFIATO, input$dataTable_jeko_rows_selected)
  })
  observeEvent(input$dataTable_jeko_rows_selected,
               {
                 showModal(modalDialog(
                   title = jekoFIATO[input$dataTable_jeko_rows_selected, "Name"],
                   plotOutput("barPlot_jeko"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
  
  
  # JURK PAGE --------------------------------------------------
  # JURK: Output summary table statistics
  output$summary_jurk <- renderTable({ summary_stat(jurkF, jurkFI, jurkFIA) })
  
  # JURK: Output pairwise plot
  output$pairwise_jurk <- renderPlot({ 
    if (input$pairwise_jurk_type == "WT + KO") {
      pattern = "^LOG2"
    } else if (input$pairwise_jurk_type == "WT") {
      pattern = "^LOG2.JURKAT.WT"
    } else (
      pattern = "^LOG2.JURKAT.KO"
    )
    if (input$pairwise_jurk_display == "Before imputation") {   # Restore NA to imputed data set
      df = select(jurkFIA, starts_with("LOG2"))
      for (i in 1:ncol(df)) df[jurkFIA[[sub("LOG2", "IMP", names(df)[i])]], i] = NA
      plot_pairs(df, pattern)
    } else plot_pairs(jurkFIA, pattern)
  })
  
  # JURK: Histogram of Imputed Value
  output$hist_jurk <- renderPlot({ hist_impute(jurkFIA, use_keep = TRUE) })
  
  # JURK: Volcano plot
  output$volcano_jurk <- renderPlot({ plot_volcano(jurkFIAT, 
                                                   ratio_cutoff = input$volcano_jurk_ratioCutoff, 
                                                   P_cutoff = input$volcano_jurk_pCutoff, 
                                                   use_adjusted = input$volcano_jurk_pvalue,
                                                   labeling = input$volcano_jurk_labels) })
  
  # JURK: Download Processed Data Table
  output$DL_jurk_processed <- downloadHandler(
    filename = function() {"processed_JURKAT.csv"},
    content = function(file) {
      write.csv(jurkFIATO, file, row.names = FALSE)
    }
  )
  
  # JURK: Download Raw Data Table
  output$DL_jurk_raw <- downloadHandler(
    filename = function() {"MaxQuant_proteinGroups_JURKAT.csv"},
    content = function(file) {
      write.csv(jurk.raw, file, row.names = FALSE)
    }
  )
  
  # JURK: Data Table
  output$dataTable_jurk <- DT::renderDataTable(
    DT::datatable(
      jurkFIATO2, 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # JURK: Filtered Table Line Plot (in modal dialog)
  output$barPlot_jurk = renderPlot({
    barplot_popup(jurkFIATO, input$dataTable_jurk_rows_selected)
  })
  observeEvent(input$dataTable_jurk_rows_selected,
               {
                 showModal(modalDialog(
                   title = jurkFIATO[input$dataTable_jurk_rows_selected, "Name"],
                   plotOutput("barPlot_jurk"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
}

shinyApp(ui = ui, server = server)


#############################################
# Step 10 - Deploy App!
#############################################
# Copy and run code in console
#library(rsconnect)
#setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Jae_surfaceomics/analysis/IFITM3/")
#deployApp()
