## Installation
For instructions on installing moETM package, please refer to its official GitHub repository <https://github.com/manqizhou/moETM>.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Prepare data
To accommodate the input requirements of moETM, we need to convert gene expression data and protein expression data to anndata format separately. Please refer to `moETM_prepare_data.ipynb` for detailed instructions.

## Run model
You can run moETM following the detailed instructions provided in `moETM.ipynb`. 

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `moETM_evaluation.ipynb` for detailed instructions.




