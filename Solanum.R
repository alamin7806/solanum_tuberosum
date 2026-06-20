library(readr)
library(magrittr)
library(dplyr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)

View(read_csv("SraRunTable.csv"))
table_826 = read_csv("SraRunTable.csv") %>%
  select(Run, BioProject, BioSample, Cultivar, Experiment,`Sample Name`,
         `SRA Study`, treatment)
View(table_826)
colnames(table_826)[c(6,7)] = c("sample_name", "SRA_study")

Sample_file_826 = paste0(pull(table_826, Run), "/quant.sf")
ref_gene = read_csv("gene_id.txt", col_names = c("transcript", "gene"))

## from transcript to gene level count

count_data_826 = tximport(files = Sample_file_826,
                          type = "salmon",
                          tx2gene = ref_gene,
                          ignoreTxVersion = TRUE)

table_826 = as.data.frame(table_826)

condition = factor(table_826$treatment)

repeat_data = c(rep("control", times = 4),
                rep("drought_stress", times = 3),
                rep("control", times = 3),
                rep("drought_stress", times = 4),
                rep("control", times = 2),
                rep("drought_stress", times = 3),
                rep("control", times = 4),
                rep("drought_stress", times = 1),
                rep("control", times = 1),
                rep("drought_stress", times = 4),
                rep("control", times = 4),
                rep("drought_stress", times = 3),
                rep("control", times = 3),
                rep("drought_stress", times = 4),
                rep("control", times = 3),
                rep("drought_stress", times = 2))
condition = factor(repeat_data)                

table_826$condition = condition
View(table_826)

## construct a DESeq DataSet from count data

deseq_data826 = DESeqDataSetFromTximport(txi = count_data_826,
                                         colData = table_826,
                                         design = ~condition)

## pre-filtering data

deseq_data826_norm = rowSums(counts(deseq_data826)) >=15
deseq_data826_norm = deseq_data826[deseq_data826_norm,]

## Differential expression analysis

dds826 = DESeq(deseq_data826_norm)
res826 = results(dds826)
vst_826 = vst(dds826)
exp_vst = assay(vst_826)

plotPCA(vst_826, intgroup = "condition") + ggtitle("GSE77826") + geom_point(size= 4)
dds826_res_df = as.data.frame(res826)

dds826_f1 = dplyr::filter(dds826_res_df, complete.cases(dds826_res_df))
dds826_f2 = dplyr::filter(dds826_f1, padj <0.05)
dds826_f3 = dplyr::filter(dds826_f2, abs(log2FoldChange) > 1)

## Prepare Machine Learning Dataset

sig_genes_826 = rownames(dds826_f3)
expr_sig826 = exp_vst[sig_genes_826,]
rf_data826 = as.data.frame(t(expr_sig826))
rf_data826$condition = table_826$condition


## Train-Test Split
set.seed(123)
train_index = createDataPartition(rf_data826$condition, p= 0.7,
                                  list = FALSE)

train_data = rf_data826[train_index,]
test_data =  rf_data826[-train_index,]

## Random Forest Model

rf_model826 = randomForest(condition ~ .,
                           data = train_data,
                           ntree = 500,
                           importance = TRUE)

print(rf_model826)

## Prediction
pred = predict(rf_model826, test_data)

## Confusion Matrix
confusionMatrix(pred, test_data$condition)

## Accuracy
Accuracy = mean(pred == test_data$condition)
Accuracy

## ROC Curve & AUC
prob = predict(rf_model826,
               test_data,
               type = "prob")

roc_obj = roc(test_data$condition, prob[,2])
plot(roc_obj)
auc(roc_obj)

## Feature Importance
importance_gene = importance(rf_model826)
top_rf826_genes = rownames(importance_gene)[
  order(importance_gene[,1], decreasing = TRUE)
][1:20]
top_rf826_genes

## Importance Plot
varImpPlot(rf_model826)

pheatmap(expr_sig826[top_rf826_genes,])

pheatmap(expr_sig826[top_rf826_genes,],
         annotation_col = table_826)


## Gene Annotation
ensembl_plants = useEnsemblGenomes(
  biomart = "plants_mart",
  dataset = "stuberosum_eg_gene"
)

annotation <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = top_rf826_genes,
  mart = ensembl_plants
)

write.csv(annotation,"top_RF_biomarkers.csv",row.names = FALSE)

annotated_826full = left_join(dds322_f1, annotation_322,
                            by = c("ensgene" = "ensembl_gene_id"))



