## Installation
For instructions on installing TotalVI package, please refer to the official documentation for more details <https://scvi-tools.org/>.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example. Additionally, we demonstrate how scMoGNN handles data containing batch information using the example of CITE-PBMC-Li: Group1→Group2 experiment in different samples scenario.

## Run model
You can run TotalVI following the detailed instructions provided in `TotalVI.ipynb`. Additionally, when the dataset contains batch information, you can follow the instructions in `TotalVI_batch.ipynb` to run TotalVI.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `TotalVI_evaluation.ipynb` for detailed instructions.

