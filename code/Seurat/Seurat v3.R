library(Seurat)
library(anndata)
library(Matrix)


# Set path
current_working_directory <- getwd()
datasets_dir <- file.path(current_working_directory, "datasets/")
datasets_dir <- normalizePath(datasets_dir)
outputs_dir <- file.path(current_working_directory, "outputs/")
outputs_dir <- normalizePath(outputs_dir)
if (!dir.exists(outputs_dir)) {
  dir.create(outputs_dir, recursive = TRUE)
}

save_dir <- file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/Seurat/")
save_dir <- normalizePath(save_dir)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


# Load h5ad data
data_train <- read_h5ad(file.path(datasets_dir, "different samples/CITE-SLN111-Gayoso/Mouse1.h5ad"))
data_test <- read_h5ad(file.path(datasets_dir, "different samples/CITE-SLN111-Gayoso/Mouse2.h5ad"))


# Convert h5ad data to SeuratObject
rna_train <- Matrix::t(data_train$X)
rna_test <- Matrix::t(data_test$X)

adt_train <- Matrix::t(data_train$obsm$protein_expression)
row.names(adt_train) <- data_train$uns$protein_name
colnames(adt_train) <- colnames(rna_train)

adt_test <- Matrix::t(data_test$obsm$protein_expression)
row.names(adt_test) <- data_test$uns$protein_name
colnames(adt_test) <- colnames(rna_test)

train_seurat <- CreateSeuratObject(counts=rna_train, assay="RNA", meta.data=data_train$obs, row.names=row.names(rna_train))
protein_train_assay <- CreateAssayObject(counts=adt_train)
train_seurat[["ADT"]] <- protein_train_assay

test_seurat <- CreateSeuratObject(counts=rna_test, assay="RNA", meta.data=data_test$obs, row.names=row.names(rna_test))
protein_test_assay <- CreateAssayObject(counts=adt_test)
test_seurat[["ADT"]] <- protein_test_assay


# Preprocess gene expression data
train_seurat <- NormalizeData(train_seurat)
test_seurat <- NormalizeData(test_seurat)

train_seurat <- FindVariableFeatures(train_seurat, nfeatures=2000)
test_seurat <- FindVariableFeatures(test_seurat, nfeatures=2000)

train_seurat <- ScaleData(train_seurat)
test_seurat <- ScaleData(test_seurat)


# Preprocess protein expression data
train_seurat <- NormalizeData(train_seurat, normalization.method="CLR", margin=2, assay="ADT")
test_seurat <- NormalizeData(test_seurat, normalization.method="CLR", margin=2, assay="ADT")


# Save preprocessed protein expression data for evaluation
test_protein_clr <- test_seurat@assays$ADT@data
write.table(test_protein_clr, file.path(save_dir, "test_protein_clr.txt"))


# Mask protein expression data in test set
test_seurat_protein <- test_seurat[["ADT"]]
test_seurat[["ADT"]] <- NULL


# Find anchors
transfer_anchors <- FindTransferAnchors(reference=train_seurat, query=test_seurat, reduction="cca", dims=1:50)   # reduction = "cca" for Seurat v3 (CCA) and reduction = "pcaproject" for Seurat v3 (PCA)


# Transfer protein data to test set
refdata <- GetAssayData(object=train_seurat, assay='ADT', slot='data')
test_protein_imputation <- TransferData(anchorset=transfer_anchors, refdata=refdata, query=test_seurat, dims=1:50, k.weight=50, 
                                        weight.reduction="pca")


# Save prediction
test_protein_prediction <- as.data.frame(test_protein_imputation@assays[["id"]]@data)
write.table(test_protein_prediction, file.path(save_dir, "test_protein_prediction (Seurat v3 (CCA)).txt"))   # CCA or PCA 







