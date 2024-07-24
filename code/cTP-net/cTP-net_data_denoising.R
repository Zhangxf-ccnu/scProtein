library(anndata)
library(SAVERX)


# Set path
current_working_directory <- getwd()
outputs_dir <- file.path(current_working_directory, "outputs/")
outputs_dir <- normalizePath(outputs_dir)
if (!dir.exists(outputs_dir)) {
  dir.create(outputs_dir, recursive = TRUE)
}

save_dir <- file.path(outputs_dir, "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/cTP-net")
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


# load data
data <- read_h5ad(file.path(datasets_dir, "different samples/CITE-SLN111-Gayoso/Mouse1.h5ad"))
rna_expression <- t(as.matrix(data$X))


# Run SAVER-X
file <- saverx(rna_expression)


# Save denoised gene expression data
denoised.data <- readRDS(file)
denoised.data <- denoised.data$estimate
write.table(denoised.data, file.path(save_dir, "Mouse1_denoised_gene_expression_data.txt"))






