# Load necessary packages
library(limma)
library(edgeR)
library(Amelia)
library(MASS)

# Load the raw expression count data and the clinical metadata
dat.raw <- read.delim("PRAD.uncv2.mRNAseq_raw_counts.txt", na.strings=c("", " ", NA))
dat.clin <- read.delim("prad_tcga_clinical_data.tsv", na.strings=c("", " ", NA))

# Number of factors, when each numerical Gleason score is taken as an independent category
table(dat.clin$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer) # 5 cagetories

# Replace the "." with "-" within the sample ID names so they match the clinical metadata sample IDs
temp_ID <- gsub("\\.", "-", colnames(dat.raw)[-1])
colnames(dat.raw)[-1] <- temp_ID

# Consider expression data from primary cancer biopsies only (e.g. Sample IDs that end with "-01")
dat.raw2 <- dat.raw[,colnames(dat.raw) %in% dat.clin$Sample.ID]
rownames(dat.raw2) <- dat.raw$HYBRIDIZATION.R # add gene names as the row names for the expression matrix
dim(dat.raw2) # There are 498 columns (sample IDs) ~ 1 sample does not seem to have expression data even though it has clinical metadata
dat.clin$Sample.ID[dat.clin$Sample.ID %in% colnames(dat.raw) == FALSE] # The Sample "TCGA-KC-A4BO-01" does not seem to have expression data
"TCGA-KC-A4BO-01" %in% colnames(dat.raw) # Gives FALSE

# Remove TCGA-KC-A4BO-01" from dat.clin
dat.clin1 <- dat.clin[-which(dat.clin$Sample.ID == "TCGA-KC-A4BO-01"),]

# Create a desing object for normalization and comparisons
class(dat.clin1$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer) # This is integer. I will treat each Gleason score as a categorical variable
G <- as.factor(dat.clin1$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer)
design <- model.matrix(~0 + G)

# Create a DGEList object using the expression data
dge <- DGEList(counts=dat.raw2)

# Filter out genes with low counts
keep <- filterByExpr(dge, design)
dge <- dge[keep,, keep.lib.sizes = FALSE]

# Scale normalization and the TMM normalization
dge <- calcNormFactors(dge)

# Perform voom normalization
v <- voom(dge, design, plot=TRUE)

# Differential expression between different gleason score groups
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(G10-G6, G9-G6, G8-G6, G7-G6, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
topTable(fit2, coef=2, adjust="BH")
topTable(fit2, coef=3, adjust="BH")
topTable(fit2, coef=4, adjust="BH")

results <- decideTests(fit2)

vennDiagram(results)

topTableF(fit2, adjust.method = "BH", number=30) # Top 30 genes associated with increasing gleason score

# Visualize expression of FA2H and CD163 genes ~ gleason score
temp <- dat.clin1[,c("Sample.ID", "Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer")]
dge.counts.t <- t(dge$counts)
fam72b_exp <- dge.counts.t[,"FAM72B|653820"]
slc7a4_exp <- dge.counts.t[,"SLC7A4|6545"]
dat.temp <- cbind(temp, fam72b_exp, slc7a4_exp)

par(mfrow=c(1,2))
boxplot(fam72b_exp ~ Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, 
        data = dat.temp, 
        main = "FAM72B gene ~ Gleason Score",
        xlab = "Gleason Score",
        ylab = "Normalized Expression")
boxplot(slc7a4_exp ~ Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, 
        data = dat.temp, 
        main = "SLC7A4 gene ~ Gleason Score",
        xlab = "Gleason Score",
        ylab = "Normalized Expression")



# Covariate analysis

# Make a heatmap of missing values to remove variables with a lot of missing values
missmap(dat.clin1, 
        main = "Missing values vs observed", 
        x.cex = 0.4, 
        y.cex = 0.4,
        y.labels = dat.clin1$Sample.ID)

# Summarize missing values in the clinical metadata
summary.missing <- sapply(dat.clin1,function(x) sum(is.na(x)))

# Remove clinical variables where there are more than 10 missing values
dat.clin2 <- dat.clin1[,!summary.missing > 10] 
missmap(dat.clin2, 
        main = "Missing values vs observed", 
        x.cex = 0.4, 
        y.cex = 0.4,
        y.labels = dat.clin2$Sample.ID)

# Remove samples in row 381 and 443 because they are missing data for multiple clinical variables
dat.clin3 <- dat.clin2[-c(381, 443),]
missmap(dat.clin3, 
        main = "Missing values vs observed", 
        x.cex = 0.4, 
        y.cex = 0.4,
        y.labels = dat.clin3$Sample.ID)

# Remove the variables with only one value again
summary.unique <- sapply(dat.clin3, function(x) length(unique(x))) # 8 columns have only 1 variable
dat.clin3 <- dat.clin3[,summary.unique > 1]
missmap(dat.clin3, 
        main = "Missing values vs observed", 
        x.cex = 0.4, 
        y.cex = 0.4,
        y.labels = dat.clin3$Sample.ID)

# Do not inlcude ID related columns in the model - so subset the metadata for molecular/genetic/clinical variables
dat.clin4 <- dat.clin3[,c(2:7, 9, 12:17, 20:21, 24:32)]
missmap(dat.clin4, 
        main = "Missing values vs observed", 
        x.cex = 0.4, 
        y.cex = 0.4,
        y.labels = dat.clin4$Sample.ID)

# Impute the missing values for the remaining columns with the missing values (using the mode for factors and mean for numeric/integer)
## Investigate the data type of the variables that have missing values
class(dat.clin4$Primary.Tumor.Laterality) # factor
class(dat.clin4$Days.to.Last.Followup) # integer
class(dat.clin4$Primary.Lymph.Node.Presentation.Assessment.Ind.3) # factor
class(dat.clin4$Fraction.Genome.Altered) # numeric
class(dat.clin4$American.Joint.Committee.on.Cancer.Tumor.Stage.Code) # factor
class(dat.clin4$Initial.pathologic.diagnosis.method) # factor
class(dat.clin4$Disease.Free..Months.) # numeric
class(dat.clin4$Disease.Free.Status) # factor
class(dat.clin4$Tissue.Prospective.Collection.Indicator) # factor
class(dat.clin4$Tissue.Retrospective.Collection.Indicator) # factor

dat.clin4$Days.to.Last.Followup <- as.numeric(dat.clin4$Days.to.Last.Followup)
## Create a function to get mode
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

dat.clin4$Primary.Tumor.Laterality[is.na(dat.clin4$Primary.Tumor.Laterality)] <- getmode(dat.clin4$Primary.Tumor.Laterality)
dat.clin4$Primary.Lymph.Node.Presentation.Assessment.Ind.3[is.na(dat.clin4$Primary.Lymph.Node.Presentation.Assessment.Ind.3)] <- getmode(dat.clin4$Primary.Lymph.Node.Presentation.Assessment.Ind.3)
dat.clin4$American.Joint.Committee.on.Cancer.Tumor.Stage.Code[is.na(dat.clin4$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)] <- getmode(dat.clin4$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
dat.clin4$Initial.pathologic.diagnosis.method[is.na(dat.clin4$Initial.pathologic.diagnosis.method)] <- getmode(dat.clin4$Initial.pathologic.diagnosis.method)
dat.clin4$Disease.Free.Status[is.na(dat.clin4$Disease.Free.Status)] <- getmode(dat.clin4$Disease.Free.Status)
dat.clin4$Tissue.Prospective.Collection.Indicator[is.na(dat.clin4$Tissue.Prospective.Collection.Indicator)] <- getmode(dat.clin4$Tissue.Prospective.Collection.Indicator)
dat.clin4$Tissue.Retrospective.Collection.Indicator[is.na(dat.clin4$Tissue.Retrospective.Collection.Indicator)] <- getmode(dat.clin4$Tissue.Retrospective.Collection.Indicator)

dat.clin4$Fraction.Genome.Altered[is.na(dat.clin4$Fraction.Genome.Altered)] <- mean(dat.clin4$Fraction.Genome.Altered,na.rm=T)
dat.clin4$Days.to.Last.Followup[is.na(dat.clin4$Days.to.Last.Followup)] <- mean(dat.clin4$Days.to.Last.Followup,na.rm=T)
dat.clin4$Disease.Free..Months.[is.na(dat.clin4$Disease.Free..Months.)] <- mean(dat.clin4$Disease.Free..Months.,na.rm=T)

missmap(dat.clin4, 
         main = "Missing values vs observed", 
         x.cex = 0.4, 
         y.cex = 0.4,
         y.labels = dat.clin4$Sample.ID)

# Then pick the model 
dat.clin4$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer <- as.factor(dat.clin4$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer)

m <- polr(Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer ~ ., 
          data=dat.clin4[,-1],
          Hess = TRUE,
          method = "logistic")
summary(m)

## store table
ctable <- coef(summary(m))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
ctable <- cbind(ctable, "p value" = p)
ctable <- as.data.frame(ctable)

## calculate odds ratios
OR <- exp(coef(m))
ctable <- cbind(ctable[-(60:63),], OR)

ctable[ctable$`p value` <= 0.05, c(4,5)]

# Run voom with the covariates selected above
G <- as.factor(dat.clin4$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer)
A <- cut(dat.clin4$Diagnosis.Age, c(40,60,80), labels = c("Y", "O")) # Age below 60 as Y (young) and above 60 as O (old)
S <- as.factor(dat.clin4$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
FGA <- dat.clin4$Fraction.Genome.Altered
TSS <- as.factor(dat.clin4$Tissue.Source.Site)

design2 <- model.matrix(~ 0 + G + A + S + FGA + TSS, data=dat.clin4)
colnames(design2)

#resubset the raw data for the sample ids that we have in dat.clin4
dat.raw3 <- dat.raw2[,colnames(dat.raw2) %in% dat.clin4$Sample.ID]

# Create a DGEList object using the expression data
dge2 <- DGEList(counts=dat.raw3)

# Filter out genes with low counts
keep2 <- filterByExpr(dge2, design2)
dge2 <- dge2[keep2,, keep.lib.sizes = FALSE]

# Scale normalization and the TMM normalization
dge2 <- calcNormFactors(dge2)

v2 <- voom(dge2, design2, plot=FALSE)

fit_new <- lmFit(v2, design2)

contrast.matrix2 <- makeContrasts(G10-G6, G9-G6, G8-G6, G7-G6, levels=design2)
fit_new2 <- contrasts.fit(fit_new, contrast.matrix2)
fit_new2 <- eBayes(fit_new2)

results2 <- decideTests(fit_new2)

vennDiagram(results2)

topTable(fit_new2, coef=1, number=30) # G10-G6 difference
topTableF(fit_new2, number=30) 

# Visualize expression of FA2H and CD163 genes ~ gleason score
temp <- dat.clin4[,c("Sample.ID", "Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer")]
dge2.counts.t <- t(dge2$counts)
fa2h_exp <- dge2.counts.t[,"FA2H|79152"]
cd163_exp <- dge2.counts.t[,"CD163|9332"]
dat.temp <- cbind(temp, fa2h_exp, cd163_exp)

par(mfrow=c(1,2))
boxplot(fa2h_exp ~ Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, 
        data = dat.temp, 
        main = "FA2H gene ~ Gleason Score",
        xlab = "Gleason Score",
        ylab = "Normalized Expression")
boxplot(cd163_exp ~ Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer, 
        data = dat.temp, 
        main = "CD163 gene ~ Gleason Score",
        xlab = "Gleason Score",
        ylab = "Normalized Expression")

### ADD SESSION INFO
sessionInfo()
