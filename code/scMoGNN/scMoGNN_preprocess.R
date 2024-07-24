library(scran)
library(SingleCellExperiment)
library(anndata)
library(Seurat)


# Set path
current_working_directory <- getwd()
datasets_dir <- file.path(current_working_directory, "datasets/")
datasets_dir <- normalizePath(datasets_dir)
outputs_dir <- file.path(current_working_directory, "outputs/")
outputs_dir <- normalizePath(outputs_dir)


# Load gene expression data
data_train <- read.table(file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/train_raw_gene_expression.txt"), 
                         sep=",", header=1, row.names=1)
data_test <- read.table(file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/test_raw_gene_expression.txt"), 
                        sep=",", header=1, row.names=1)


# Load clustering results
meta_train <- read.table(file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/train_cell_metadata.txt"), 
                         sep=",", header=1, row.names=1)
meta_test <- read.table(file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/test_cell_metadata.txt"), 
                        sep=",", header=1, row.names=1)


# Convert gene expression data to SingleCellExperiment object
sce_train <- SingleCellExperiment(list(counts=data_train))
sce_test <- SingleCellExperiment(list(counts=data_test))


# Calculate the size factor for normalization
clusters_train <- meta_train$groups
clusters_test <- meta_test$groups

sum_fac_train <- computeSumFactors(sce_train, clusters=clusters_train, min.mean=0.1)
sum_fac_test <- computeSumFactors(sce_test, clusters=clusters_test, min.mean=0.1)

size_fac_train <- sizeFactors(sum_fac_train)
size_fac_test <- sizeFactors(sum_fac_test)


# Preprocess gene expression data
data_train_normalized <- t(data_train) / size_fac_train
data_test_normalized <- t(data_test) / size_fac_test

data_train_normalized <- log(data_train_normalized+1)
data_test_normalized <- log(data_test_normalized+1)


# Save preprocessed gene expression data
write.table(data_train_normalized, 
            file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/train_normalized_gene_expression.txt"))
write.table(data_test_normalized,
            file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/test_normalized_gene_expression.txt"))


# Preprocess protein gene expression data
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


# Preprocess protein expression data
train_seurat <- NormalizeData(train_seurat, normalization.method="CLR", margin=2, assay="ADT")
test_seurat <- NormalizeData(test_seurat, normalization.method="CLR", margin=2, assay="ADT")


# Save preprocessed protein expression data
train_protein_clr <- train_seurat@assays$ADT@data
test_protein_clr <- test_seurat@assays$ADT@data
write.table(train_protein_clr, file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/train_protein_clr.txt"))
write.table(test_protein_clr, file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN/test_protein_clr.txt"))







