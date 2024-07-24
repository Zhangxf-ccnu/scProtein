## Installation
For instructions on installing cTP-net package, please refer to its official GitHub repository <https://github.com/zhouzilu/ctpnetpy>. cTP-net uses SAVER-X to denoise the raw gene expression data. For installation instructions, please refer to its official GitHub repository <https://github.com/jingshuw/SAVERX>.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Denoise gene expression data
Before running the model, cTP-net uses SAVER-X to denoise the raw gene expression data. Please refer to `cTP-net_data_denoising.R` for detailed instructions.

## Run model
You can run cTP-net following the detailed instructions provided in `cTP-net.ipynb`. 

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `cTP-net_evaluation.ipynb` for detailed instructions.



