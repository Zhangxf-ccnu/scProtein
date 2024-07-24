## Installation
scMoGNN is based on the dance package for Python. Please refer to the official GitHub repository <https://github.com/OmicsML/dance> for instructions on how to install the dance package.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example. Additionally, we demonstrate how scMoGNN handles data containing batch information using the example of CITE-PBMC-Li: Group1→Group2 experiment in different samples scenario.

## Cluster gene expression data
First, please run `scMoGNN_cluster.ipynb` to obtain the clustering results required for preprocessing gene expression data.

## Preprocess data
Please run `scMoGNN_preprocess.R` to preprocess the gene expression data and protein expression data.

## Prepare data
To accommodate the input requirements of scMoGNN, you need to convert gene expression data and protein expression data to anndata format separately. Please refer to `scMoGNN_prepare_data.ipynb` for detailed instructions.

## Run model
Babel is run using the command line. Please make sure to set your working directory to `your customized path/code/scMoGNN` first. Use the command line to train the model:
```python
python scMoGNN.py --prefix "CITE-SLN111-Gayoso-Mouse1toMouse2" --subtask "openproblems_bmmc_cite_phase2_rna" --log_folder "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN" --model_folder "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN" --result_folder "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMoGNN" --pathway --no_batch_features -inres -sb -sf -ws                                                                      
```
The `prefix` parameter specifies the name of the experiment and is used for naming files when saving results. The `log_folder` parameter specifies the path for saving log files. The `model_folder` parameter specifies the path for saving trained model. The `result_folder` parameter specifies the path for saving imputed protein expression for the test set. Ensure the `no_batch_features` parameter is included in the command line to indicate that `no_batch_features` is set to True.
Additionally, use the following command line to handle datasets with batch information:
```python
python scMoGNN_batch.py --prefix "CITE-PBMC-Li-Group1toGroup2" --subtask "openproblems_bmmc_cite_phase2_rna" --log_folder "your customized path/outputs/different samples/CITE-PBMC-Li-Group1toGroup2/scMoGNN" --model_folder "your customized path/outputs/different samples/CITE-PBMC-Li-Group1toGroup2/scMoGNN" --result_folder "your customized path/outputs/different samples/CITE-PBMC-Li-Group1toGroup2/scMoGNN" --pathway -inres -sb -sf -ws --batch_information "donor"                                                                      
```
Ensure the `no_batch_features` parameter is not included in the command line to indicate that `no_batch_features` is set to False. The `batch_information` parameter is the column name in the dataset's `obs` that contains sample batch information.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `scMoGNN_evaluation.ipynb` for detailed instructions.


