## Installation
For instructions on installing Seurat package, please refer to the official documentation for more details <https://satijalab.org/seurat/articles/install_v5>.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Run model
You can run Seurat v3 and Seurat v4 by following the instructions in `Seurat v3.R` and `Seurat v4.R`, respectively. You can choose to use CCA or PCA for dimensionality reduction by specifying the `reduction` parameter in the `FindTransferAnchors` function.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `Seurat_evaluation.ipynb` for detailed instructions.