# IFITM3-surfaceomics
Surface proteomics analysis on IFITM3-deficient cell lines.

## Motivation
The aim of this experiment is to identify membrane proteins whose abundances are affected as a result of IFITM3 
(Interferon-induced transmembrane protein 3) deletion. IFITM3 is a single-pass transmembrane protein localized in 
the plasma membrane and endosomes. It has been shown to play a role in maintaining cholesterol distribution in the 
cell membrane for the formation of lipid rafts, which in turn serve as anchoring platforms for many membrane proteins. 
Using the biotin-hydrazide-based cell surfaceome method, we hope to characterize surface receptors that are lost when IFITM3 
is knocked out.

## Interactive Web App
An interactive web tool built using ShinyApp in R is available for data visualization. Visit the site [here](https://tony-lin.shinyapps.io/ifitm3/).

## File Description
+ *Jae_surfaceome_analysis.R*: R analysis script.
+ *Jae_workspace.RData*: Saved workspace after running the script above.
+ *IFITM3/app.R*: Source code behind interactive app.
+ *proteinGroups-jeko.txt*: MaxQuant output for JeKo-1.
+ *proteinGroups-jurkat.txt*: MaxQuant output for Jurkat.
+ *membrane-proteins-annotation.csv*: Annotation file for identifying membrane proteins.
