library(data.table)
library(Seurat)
library(stringr)
library(caret)
library(ggplot2)

set.seed(42)

setwd("/stornext/Bioinf/data/lab_phipson/givanna/ozsinglecell_2022/")

dat <- readRDS("aml_dataset/figshare/Healthy.rds")

Cells(dat)
# 49507 cells
ncol(dat)

# 568 RNA + ADT
nrow(dat)

rownames(dat)
colnames(dat)

# Antibodies seem to have AB prefix
ab_names <- rownames(dat)[str_which(rownames(dat), ".*-AB")]

# Strange shit has happened where cell type for Cell_Type until Cell_Type3 are only for samples in BM1
table(dat$Cell_Type, dat$Batch)
table(dat$Cell_Type2, dat$Batch)
table(dat$Cell_Type3, dat$Batch)
table(dat$Cell_Type4, dat$Batch)

# Maybe the higher level annotation is applicable to all samples? 
# If the sum of BM1 for Cell_Type is the same as the sum of all samples in Cell_Type4, then we have the answer.
ct1 <- data.table(table(dat$Cell_Type, dat$Batch))
ct2 <- data.table(table(dat$Cell_Type4, dat$Batch))
sum(ct1$N)
sum(ct2[ct2$V2 == 'BM1']$N)
# Nope, the annotation for Cell_Type until Cell_Type3 are only applicable for BM1.

## Let's pick BM1 sample for now and use the highest level annotation
# Extract the annotation and the antibodies into data table

ab_dat <- data.table(FetchData(dat, vars = c(ab_names, "Cell_Type", "Batch")))

# Keep just BM1 and remove cells that are not assigned any annotation
# 9751 cells
ab_dat <- ab_dat[ab_dat$Batch == 'BM1' & !is.na(ab_dat$Cell_Type)]

# Remove ADT which are not expressed at all by any cells (i.e. sum == 0)
which(colSums(ab_dat[, 1:105]) == 0, TRUE) # CD138 and CD99

ab_dat <- ab_dat[, c(which(colSums(ab_dat[, 1:105]) > 0), 106, 107), with=FALSE]
ncol(ab_dat)

cell_type <- as.factor(as.character(ab_dat$Cell_Type))

# Remove the cell type and batch info (the last 2 columns)
ab_dat <- ab_dat[, 1:103]

# use 5 fold to do parameter tuning
seeds <- vector(mode = "list", length = 26)
for(i in 1:51) seeds[[i]] <- rep(42, 3)

tr_ctrl <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        seeds = seeds)

# random subsample just to test
# rand_sub <- sample(x = c(1:nrow(ab_dat)), size=1000)
# ab_dat_sub <- ab_dat[rand_sub,]
# cell_type_sub <- as.factor(as.character(cell_type[rand_sub]))

## Use random forest
# TODO have to work out the pre-processing. 
# Center and scale have issue with features with small variance

# Subsetted data just for testing.
# class_model_sub <- train(x=ab_dat_sub, y=cell_type_sub, 
#                      method='rf',
#                      trControl = tr_ctrl,
#                      allowParallel=TRUE,
#                      seeds = seeds)

class_model <- train(x=ab_dat, y=cell_type, 
                     method='rf',
                     trControl = tr_ctrl,
                     allowParallel=TRUE,
                     seeds = seeds)

# feat_imp_sub <- data.table(varImp(class_model_sub)$importance, keep.rownames = TRUE)
feat_imp <- data.table(varImp(class_model)$importance, keep.rownames = TRUE)
# 
# feat_imp <- merge.data.table(feat_imp, feat_imp_sub, 
#                              suffixes = c("_all_cells", "_subsampled_1000"),
#                              by = 'rn')
# 
# # Interestingly, there seem to be some correspondence between subsampled and full data.
# ggplot(feat_imp, aes(x=feat_imp$Overall_all_cells, y=feat_imp$Overall_subsampled_1000)) +
#   geom_point() +
#   geom_text(label=feat_imp$rn)

# Turn scientific notation off
feat_imp <- feat_imp[order(-rank(Overall))]
feat_imp$Rank_all_cells <- c(1:nrow(feat_imp))
feat_imp$Overall <- format(feat_imp$Overall, scientific = FALSE)


# Some plot
plot(varImp(class_model))

conf_mat <- data.table(class_model$finalModel$confusion, keep.rownames = TRUE)
fwrite(conf_mat, "confusion_matrix_rf.csv")

# Rough look on the top 10 ADTs found by RF to be important
dat_sub <- subset(dat, subset = Batch == 'BM1')

# Assuming we can use a panel of 10. What would it be?
FeaturePlot(dat_sub, features = c('CD11b-AB', 'CD8-AB', 'CD133-AB', 'IgG-AB', 'CD45RA-AB', 'CD3-AB',
                                  'CD19-AB', 'CD4-AB', 'CD38-AB', 'CD34-AB'))
Idents(dat_sub) <- dat_sub$Cell_Type 
DimPlot(dat_sub, reduction = "UMAPni", label = TRUE, pt.size = 0.5) + NoLegend()

rmarkdown::render("ch1.R", "pdf_document")

