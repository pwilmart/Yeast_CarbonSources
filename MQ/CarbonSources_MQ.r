
# load libraries first
library("psych")
library("tidyverse")
library("stringr")
library("gridExtra")
library("scales")
library("edgeR")
library("limma")

# load the data file into a data tibble 
raw_data <- read_tsv(file = "R-input.txt")

# separate row identifiers from the data
accession <- raw_data$Accession
tmt_raw <- raw_data %>% select(-c(Accession))

# let's use more descriptive column names
colnames(tmt_raw) <- c("gal-1", "gal-2", "gal-3", 
                        "glu-1", "glu-2", "glu-3",
                        "raf-1", "raf-2", "raf-3")
nrow(tmt_raw)

# let's see what the starting data look like
# set a 3 by 3 color vector
color_vector <- rep(c("red", "blue", "green"), each = 3)

# boxplots of RAW log intensities
boxplot(log10(tmt_raw), col = color_vector, 
        notch = TRUE, main = "RAW data: Gal (red), Glu (blue), and Raf(green)",
        xlab = "TMT Samples", ylab = "log of Intensity")

SL_Norm <- function(df) {
    # This makes each channel sum to the average grand total
        # df - data frame of TMT intensities
        # returns a new data frame with normalized values
    norm_facs <- mean(c(colSums(df))) / colSums(df)
    cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
    df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
}

# normalize the raw data
tmt_sl <- SL_Norm(tmt_raw)

# let's see what the SL normalized data look like
boxplot(log10(tmt_sl), col = color_vector, 
        notch = TRUE, main = "SL Norm data: Gal (red), Glu (blue), and Raf(green)",
        xlab = "TMT Samples", ylab = "log2 of Intensity")

# and check clustering
plotMDS(log10(tmt_sl), col = color_vector, main = "Clustering of SL Norm data")

# load data, study design, and row labels into edgeR object
group <- rep(c("gal", "glu", "raf"), each = 3)
y <- DGEList(counts = tmt_raw, group = group, genes = accession)

# run TMM normalization
y <- calcNormFactors(y) 
y$samples

apply_tmm_factors <- function(y) {
    # computes the tmm normalized data from the DGEList object
        # y - DGEList object
        # returns a dataframe with normalized intensities
    
    # compute grand total (library size) scalings
    lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size

    # the TMM factors are library adjustment factors (so divide by them)
    norm_facs <- lib_facs / y$samples$norm.factors

    # compute the normalized data as a new data frame
    tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
    colnames(tmt_tmm) <- str_c(colnames(tmt_raw), "_tmm")
    
    # return the data frame
    tmt_tmm
}
tmt_tmm <- apply_tmm_factors(y)

# look at intensity distributions across samples after TMM
boxplot(log10(tmt_tmm), 
        col = color_vector,
        xlab = "Samples", ylab = "Reporter Intensity", 
        main = "After TMM Normalization", notch = TRUE)

# check clustering after TMM with MDS plot
plotMDS(y, col = color_vector, main = "Clustering after TMM normalization")

# compare all samples to each other
pairs.panels(log10(tmt_tmm), method = "spearman", 
             lm = TRUE, main = "All versus all (normalized)")

# save the indexes of each condition
gal <- 1:3
glu <- 4:6
raf <- 7:9

CV <- function(df) {
    # Computes CVs of data frame rows
        # df - data frame, 
        # returns vector of CVs (%)
    ave <- rowMeans(df)    # compute averages
    sd <- apply(df, 1, sd) # compute standard deviations
    cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

labeled_boxplot <- function(df, ylim, title) {
    # Makes a box plot with the median value labeled
        # df - data frame with data to compute CVs of
        # ylim - upper limit for y-axis
        # title - plot title
    cv = CV(df)
    boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
    text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
}

# make CV distributions for each condition before and after normalization
par(mfrow = c(2, 3))

labeled_boxplot(tmt_raw[gal], 30, "Galactose RAW")
labeled_boxplot(tmt_raw[glu], 30, "Glucose RAW")
labeled_boxplot(tmt_raw[raf], 30, "Raffinose RAW")
labeled_boxplot(tmt_tmm[gal], 30, "Galactose TMM")
labeled_boxplot(tmt_tmm[glu], 30, "Glucose TMM")
labeled_boxplot(tmt_tmm[raf], 30, "Raffinose TMM")

par(mfrow = c(1, 1))

# print CV cutoffs that cature 90% of the proteins
cat("Gal - 90% of CVs less than:", round(quantile(CV(tmt_tmm[gal]), probs = 0.9), 2), "%\n")
cat("Glu - 90% of CVs less than:", round(quantile(CV(tmt_tmm[glu]), probs = 0.9), 2), "%\n")
cat("Raf - 90% of CVs less than:", round(quantile(CV(tmt_tmm[raf]), probs = 0.9), 2), "%\n")

# make a tidy data frame for plotting
cv_gal <- data.frame(cv = CV(tmt_tmm[gal]), sugar = "gal", stringsAsFactors = FALSE)
cv_glu <- data.frame(cv = CV(tmt_tmm[glu]), sugar = "glu", stringsAsFactors = FALSE)
cv_raf <- data.frame(cv = CV(tmt_tmm[raf]), sugar = "raf", stringsAsFactors = FALSE)
cv_long <- rbind(cv_gal, cv_glu, cv_raf)

# density plots
ggplot(cv_long, aes(x = cv, fill = sugar)) +
  geom_density(alpha = 0.3) +
  coord_cartesian(xlim = c(0, 25)) +
  ggtitle("CV distributions (TMM Norm)")

# compute the shared variance estimates and plot variance trends
y <- estimateDisp(y)
plotBCV(y, main = "Carbon Sources")

collect_results <- function(df, tt, x, xlab, y, ylab) {
    # Computes new columns and extracts some columns to make results frame
        # df - data in data.frame
        # tt - top tags from edgeR test
        # x - columns for first condition
        # xlab - label for x
        # y - columns for second condition
        # ylab - label for y
        # returns a new dataframe
    
    # condition average vectors
    ave_x <- rowMeans(df[x])
    ave_y <- rowMeans(df[y])
    
    # FC, direction, candidates
    fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
    direction <- ifelse(ave_y > ave_x, "up", "down")
    candidate = cut(tt$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                    labels = c("high", "med", "low", "no"))
    
    # make data frame
    temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                          PValue = tt$PValue, FDR = tt$FDR, 
                                          ave_x = ave_x, ave_y = ave_y, 
                                          direction = direction, candidate = candidate, 
                                          Acc = tt$genes)) 
    
    # fix column headers for averages
    names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
    
    temp # return the data frame
}

# compute the exact test models, p-values, FC, etc.
et <- exactTest(y, pair = c("gal", "glu"))

# see which proteins have the smallest p-values
topTags(et)$table

# make the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
gal_glu <- collect_results(tmt_tmm, tt, gal, "gal", glu, "glu")

pvalue_plots <- function(results, ylim, title) {
    # Makes p-value distribution plots
        # results - results data frame
        # ylim - ymax for expanded view
        # title - plot title
    p_plot <- ggplot(results, aes(PValue)) + 
        geom_histogram(bins = 100, fill = "white", color = "black") +
        geom_hline(yintercept = mean(hist(results$PValue, breaks = 100, 
                                     plot = FALSE)$counts[26:100]))

    # we will need an expanded plot
    p1 <- p_plot + ggtitle(str_c(title, " p-value distribution"))
    p2 <- p_plot + coord_cartesian(xlim = c(0, 1.0), ylim = c(0, 50)) + ggtitle("p-values expanded")
    grid.arrange(p1, p2, nrow = 2) # from gridExtra package
}

# check the p-value distribution
pvalue_plots(gal_glu, 50, "Gal vs Glu")

# see how many up and down candidates (10% FDR)
summary(decideTests(et, p.value = 0.10))

# use function from limma for MD plot
plotMD(et, main = "Gal vs Glu TMM Normalized", p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# see how many candidates are in each category
gal_glu %>% count(candidate)

log2FC_plots <- function(results, range, title) {
    # Makes faceted log2FC plots by candidate
        # results - results data frame
        # range - plus/minus log2 x-axis limits
        # title - plot title
    ggplot(results, aes(x = logFC, fill = candidate)) +
        geom_histogram(binwidth=0.1, color = "black") +
        facet_wrap(~candidate) +
        ggtitle(title) + 
        coord_cartesian(xlim = c(-range, range))
}

# can also look at log2FC distributions
log2FC_plots(gal_glu, 3, "LogFC by candidate for Gal vs Glu")

transform <- function(results, x, y) {
    # Make data frame with some transformed columns
        # results - results data frame
        # x - columns for x condition
        # y - columns for y condition
        # return new data frame
    df <- data.frame(log10((results[x] + results[y])/2), 
                     log2(results[y] / results[x]), 
                     results$candidate,
                     -log10(results$FDR))
    colnames(df) <- c("A", "M", "candidate", "P")
    
    df # return the data frame
}

MA_plots <- function(results, x, y, title) {
    # makes MA-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots 
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # 2-fold change lines
    ma_lines <- list(geom_hline(yintercept = 0.0, color = "black"),
                     geom_hline(yintercept = 1.0, color = "black", linetype = "dotted"),
                     geom_hline(yintercept = -1.0, color = "black", linetype = "dotted"))

    # make main MA plot
    ma <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("logFC (", y, "/", x, ")")) +
        scale_x_continuous("Ave_intensity") +
        ggtitle(title) + 
        ma_lines
    
    # make separate MA plots
    ma_facet <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("log2 FC (", y, "/", x, ")")) +
        scale_x_continuous("log10 Ave_intensity") +
        ma_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)"))

    # make the plots visible
    print(ma)
    print(ma_facet)
}    

scatter_plots <- function(results, x, y, title) {
    # makes scatter-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots
    
    # 2-fold change lines
    scatter_lines <- list(geom_abline(intercept = 0.0, slope = 1.0, color = "black"),
                          geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          scale_y_log10(),
                          scale_x_log10())

    # make main scatter plot
    scatter <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        ggtitle(title) + 
        scatter_lines

    # make separate scatter plots
    scatter_facet <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scatter_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)")) 

    # make the plots visible
    print(scatter)
    print(scatter_facet)
}

volcano_plot <- function(results, x, y, title) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # title - plot title string
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        ggtitle(str_c(title, " Volcano Plot"))
}

# do the gal vs glu MA visualizations
MA_plots(gal_glu, "ave_gal", "ave_glu", "Galactose vs Glucose")

# same comparison with scatter plots
scatter_plots(gal_glu, "ave_gal", "ave_glu", "Galactose vs Glucose")

# the ubiquitous volcano plot
volcano_plot(gal_glu, "ave_gal", "ave_glu", "Galactose vs Glucose")

# compare the conditions to each other
pairs.panels(log10(tmt_tmm[c(gal, glu)]), method = "spearman", 
             lm = TRUE, main = "Galactose versus Glucose")

# compute the exact test models, p-values, FC, etc.
et <- exactTest(y, pair = c("gal", "raf"))

# see which proteins have the smallest p-values
topTags(et)$table

# get the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
gal_raf <- collect_results(tmt_tmm, tt, gal, "gal", raf, "raf")

# check the p-value distribution
pvalue_plots(gal_raf, 50, "Gal vs Raf")

# see how many up and down candidates (10% FDR)
summary(decideTests(et, p.value = 0.10))

# use function from limma for MD plot
plotMD(et, main = "Gal vs Raf TMM Normalized", p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# see how many candidates are in each category
gal_raf %>% count(candidate)

# can also look at log2FC distributions
log2FC_plots(gal_raf, 3, "LogFC by candidate for Gal vs Raf")

# make the MA plots
MA_plots(gal_raf, "ave_gal", "ave_raf", "Galactose vs Raffinose")

# scatter plots
scatter_plots(gal_raf, "ave_gal", "ave_raf", "Galactose vs Raffinose")

volcano_plot(gal_raf, "ave_gal", "ave_raf", "Galactose vs Raffinose")

# compare the conditions to each other
pairs.panels(log10(tmt_tmm[c(gal, raf)]), method = "spearman", 
             lm = TRUE, main = "Galactose versus Raffinose")

# compute the exact test models, p-values, FC, etc.
et <- exactTest(y, pair = c("glu", "raf"))

# see which proteins have the smallest p-values
topTags(et)$table

# get the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
glu_raf <- collect_results(tmt_tmm, tt, glu, "glu", raf, "raf")

# check the p-value distribution
pvalue_plots(glu_raf, 50, "Glu vs Raf")

# see how many up and down candidates (10% FDR)
summary(decideTests(et, p.value = 0.10))

# use function from limma for MD plot
plotMD(et, main = "Glu vs Raf TMM Normalized", p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# see how many candidates are in each category
glu_raf %>% count(candidate)

# can also look at log2FC distributions
log2FC_plots(glu_raf, 3, "LogFC by candidate for Glu vs Raf")

# MA plots for glu vs raf
MA_plots(glu_raf, "ave_glu", "ave_raf", "Glucose vs Raffinose")

# scatter plots
scatter_plots(glu_raf, "ave_glu", "ave_raf", "Glucose vs Raffinose")

volcano_plot(glu_raf, "ave_glu", "ave_raf", "Glucose vs Raffinose")

# compare the conditions to each other
pairs.panels(log10(tmt_tmm[c(glu, raf)]), method = "spearman", 
             lm = TRUE, main = "Glucose versus Raffinose")

# save the results file to add back to the main spreadsheet
results <- cbind(gal_glu, gal_raf, glu_raf)
write.table(results, "CarbonSources_results.txt", sep = "\t", row.names = FALSE, na = " ")

# log the session information
sessionInfo()


