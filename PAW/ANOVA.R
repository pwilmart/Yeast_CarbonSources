library(tidyverse)
library(edgeR)
library(HybridMTest)

# load the data files into data frames. 
raw_data <- read_tsv(file = "edgeR_input.txt")
accession <- raw_data$Accession
tmt_raw <- raw_data[, -1]
# let's use more descriptive column names
colnames(tmt_raw) <- c("gal-1", "gal-2", "gal-3", 
                       "glu-1", "glu-2", "glu-3",
                       "raf-1", "raf-2", "raf-3")

# load data into edgeR structures and run TMM norm
group <- factor(rep(c("gal", "glu", "raf"), each = 3))
y <- DGEList(counts = tmt_raw, group = group, genes = accession)

# run TMM normalization
y <- calcNormFactors(y) 
y$samples

# Compute the normalized counts (start with tmt_data)
# sample loading adjusts each channel to the same average total
lib_facs <- mean(colSums(tmt_raw)) / colSums(tmt_raw)

# the TMM factors are library adjustment factors (so divide by them)
norm_facs <- lib_facs / y$samples$norm.factors

# compute the normalized data as a new data frame
tmt_tmm <- sweep(tmt_raw, 2, norm_facs, FUN = "*")
colnames(tmt_tmm) <- paste(colnames(tmt_raw), "tmm", sep = "_")

# set up design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# dispersion
y <- estimateDisp(y, design, robust = TRUE)

# run glm models
fit <- glmQLFit(y, design, robust = TRUE)

# make the contrast
con <- makeContrasts(gal.glu = gal - glu,
                     gal.raf = gal - raf,
                     glu.raf = glu - raf, levels = design)

# get the test results
anov <- glmQLFTest(fit, contrast = con)
topTags(anov)

# do multiple testing corrections
FDR3 <- p.adjust(anov$table$PValue, "BH")
length(FDR3[FDR3 < 0.01])

FDR4 <- p.adjust(anov$table$PValue, "bonferroni")
length(FDR4[FDR4 < 0.01])


###################################################################
# do old school one-way ANOVA
anova_test <- row.oneway.anova(log10(tmt_tmm), group)

# BH correction
FDR1 <- p.adjust(anova_test$pval, "BH")
length(FDR1[FDR1 < 0.01])

# Bonferroni correction
FDR2 <- p.adjust(anova_test$pval, "bonferroni")
length(FDR2[FDR2 < 0.01])

