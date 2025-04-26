# Load the required libraries for analysis
library("xlsx")               # For reading Excel files
library("preprocessCore")     # For quantile normalization
library("msImpute")           # For missing data imputation in mass spectrometry-based proteomics
library("survival")           # For survival analysis
library("survminer")          # For visualization of survival analysis
library("RegParallel")        # For parallelized regression analysis
library("stringr")            # For advanced string manipulation
library("org.Hs.eg.db")       # For gene annotation (human genome database)
library("ReactomePA")         # For pathway enrichment analysis using Reactome
library("pathfindR")          # For pathway enrichment analysis

# Set the working directory and list files in the specified folder
setwd("C://Users/kptl632/OneDrive - AZCollaboration/Documents/personal/task_dfkz/")
list.files("Bioinformatic_task/") # Listing files in the 'Bioinformatic_task' directory

# Load the input data
prot_abun = read.csv("Bioinformatic_task/Protein_abundance.tsv", 
                     sep = "\t", 
                     stringsAsFactors = F, 
                     check.names = F, 
                     row.names = 1) # Protein abundance dataset, rows represent proteins

metadata = read.xlsx(file = "Bioinformatic_task/sampleAnnotation.xls", 
                     sheetIndex = 1) # Metadata for the samples

# Basic checks on the data
## Calculate the percentage of missing data in the protein abundance matrix
print(paste0("Percentage missing = ", 
             round(100 * (sum(is.na(prot_abun)) / (nrow(prot_abun) * ncol(prot_abun))), 2), " %"))

## Check dimensions of protein abundance and metadata
dim(prot_abun)
dim(metadata)

# Identify a sample present in metadata but missing in protein abundance
print("Metadata has 1 sample extra")
missing_sample = metadata$sample.ID[!metadata$sample.ID %in% colnames(prot_abun)]
metadata[metadata$sample.ID == missing_sample, ]

# Filter metadata to match the sample IDs in protein abundance data
print("Most of the metadata from the missing sample is not available")
metadata_filter = metadata[metadata$sample.ID %in% colnames(prot_abun), ]
row.names(metadata_filter) = metadata_filter$sample.ID
metadata_filter$sample.ID = NULL

# Reorder protein abundance data to match the metadata
prot_abun_order = prot_abun[, row.names(metadata_filter)]

# Data processing
## Apply quantile normalization
qn_pa = normalize.quantiles(as.matrix(prot_abun_order))
row.names(qn_pa) = row.names(prot_abun_order)
colnames(qn_pa) = colnames(prot_abun_order)

## Log-transform the normalized data
log_qn_pa = log2(qn_pa)

## Check rows with missing values and filter out features with excessive missing data
max(rowSums(is.na(log_qn_pa)))                 # Maximum number of missing values in any row
table(rowSums(is.na(log_qn_pa)))              # Distribution of missing values across rows
row.names(log_qn_pa[rowSums(is.na(log_qn_pa)) >= 46, ]) # Features with excessive missing data

filtered = log_qn_pa[rowSums(is.na(log_qn_pa)) <= 45, ] # Retain features with fewer than 45 missing values

# Impute missing values
imputed = msImpute(filtered, method = "v1") # General missing value imputation
imputed_mnar = msImpute(filtered, group = as.character(metadata_filter$died)) # Imputation considering missing not at random (MNAR)

# Check minimum values to assess signal strength
min(na.omit(prot_abun_order))
min(na.omit(filtered))
min(imputed_mnar)
min(imputed)

# Quality assessment
## Prepare data for Cox Proportional Hazards (Cox PH) model
coxdata = cbind(t(imputed), metadata_filter)
coxdata$time = as.numeric(coxdata$last.known.alive - coxdata$date.of.diagnosis) # Time to event
coxdata$total.protein.concentration = as.numeric(coxdata$total.protein.concentration) # Convert concentration to numeric
coxdata$freeThawCycle = as.factor(coxdata$freeThawCycle) # Convert freeze-thaw cycle to factor
coxdata$batch = as.factor(replace(coxdata$batch, coxdata$batch == "0", "7")) # Recode batch 0 as 7
coxdata$operator = as.factor(replace(coxdata$operator, coxdata$operator == "MG/FMA", "MG")) # Recode operator values
coxdata$event = ifelse(coxdata$died == TRUE, 2, 1) # Event indicator (1 = alive, 2 = died)

# Perform univariate Cox PH analysis for technical factors
batch_cox = coxph(Surv(time, event) ~ batch, data = coxdata)
op_cox = coxph(Surv(time, event) ~ operator, data = coxdata)
thaw_cox = coxph(Surv(time, event) ~ freeThawCycle, data = coxdata)
conc_cox = coxph(Surv(time, event) ~ total.protein.concentration, data = coxdata)

# Function to summarize Cox PH results
cox_summary = function(x) { 
  x = summary(x)
  p.value = signif(x$wald["pvalue"], digits = 2)
  wald.test = signif(x$wald["test"], digits = 2)
  beta = signif(x$coef[1], digits = 2) # Coefficient beta
  HR = signif(x$coef[2], digits = 2)   # Hazard ratio
  res = c(beta, HR, wald.test, p.value)
  names(res) = c("beta", "HR", "wald.test", "p.value")
  return(res)
}

# Summarize univariate Cox PH results
univariate_summary = data.frame(
  rbind(
    cox_summary(op_cox),
    cox_summary(thaw_cox),
    cox_summary(conc_cox),
    cox_summary(batch_cox)
  ), row.names = c("operator", "thaw", "protein concentration", "batch")
)

# Visualize survival curves
ggsurvplot(survfit(op_cox), data = coxdata, conf.int = FALSE)
new_df = data.frame("operator" = c("CB", "MG"), 
                    "freeThawCycle" = c("0"), 
                    "total.protein.concentration" = mean(coxdata$total.protein.concentration), 
                    "batch" = "3")
fit = survfit(op_cox, newdata = new_df)
op_survplot = ggsurvplot(fit, legend.labs = c("CB", "MG"), data = coxdata, ggtheme = theme_minimal(), conf.int = FALSE, pval = TRUE)
op_survplot

# PCA analysis to assess the effect of technical factors
t_imp = t(imputed)
pca_imp = prcomp(t_imp, scale. = TRUE, center = TRUE)
autoplot(pca_imp, data = coxdata, colour = 'batch') # Color by batch
autoplot(pca_imp, data = coxdata, colour = 'died') # Color by survival status
# No obvious cluster separation suggests technical factors do not drive clustering

# Identify prognostic protein markers using parallelized regression
protein_coxdata = coxdata[, !colnames(coxdata) %in% c("batch", "operator", "total.protein.concentration", "conc", "died", "last.known.alive", "date.of.diagnosis", "freeThawCycle")]
colnames(protein_coxdata) = make.names(colnames(protein_coxdata)) # Ensure column names are syntactically valid

prognosis = RegParallel(
  data = protein_coxdata,
  formula = 'Surv(time, event) ~ [*]',
  FUN = function(formula, data) coxph(formula = formula, data = data, ties = 'breslow', singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(protein_coxdata)[1:(ncol(protein_coxdata) - 2)],
  blocksize = 200,
  cores = 3,
  nestedParallel = FALSE,
  conflevel = 95
)

# Filter significant features
prognosis2 = prognosis[!is.na(prognosis$P), ]
prognosis2 = prognosis2[order(prognosis2$LogRank, decreasing = FALSE), ]
signif_feature = subset(prognosis2, LogRank < 0.01)
signif_feature2 = subset(prognosis2, LogRank < 0.05)

# Map significant features to genes and perform pathway enrichment
prognostic_proteins = str_split_fixed(signif_feature2$Variable, pattern = "\\.", 3)[, 2]
prognostic_gene = str_remove(str_split_fixed(signif_feature2$Variable, pattern = "\\.", 3)[, 3], "_HUMAN")

# Map genes to Entrez IDs for pathway analysis
hs = org.Hs.eg.db
my.symbols = prognostic_gene
entrez_map = select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
prognostic_gene_entrez = entrez_map[!is.na(entrez_map$ENTREZID), 2]

# Perform pathway enrichment analysis using ReactomePA
x = enrichPathway(gene = prognostic_gene_entrez, pvalueCutoff = 0.2, qvalueCutoff = 0.2, readable = TRUE)
reactome_pa_res = x@result
print(reactome_pa_res[1:15, c(2, 5, 6, 8)]) # Print top 15 enriched pathways

# Perform pathway enrichment analysis using pathfindR
pfr_input = signif_feature2[, c(1, 3, 9)]
pfr_input$Variable = str_remove(str_split_fixed(pfr_input$Variable, pattern = "\\.", 3)[, 3], "_HUMAN")
colnames(pfr_input) = c("Gene.symbol", "logFC", "adj.P.Val")
pfr_out = run_pathfindR(pfr_input)
