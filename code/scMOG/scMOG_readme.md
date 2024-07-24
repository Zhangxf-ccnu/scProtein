## Installation
For instructions on installing scMOG package, please refer to its official GitHub repository <https://github.com/GaoLabXDU/scMOG>.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Prepare data
To accommodate the input requirements of scMOG, we need to convert both the gene expression data and protein expression data from the h5ad file to CSV format. Please refer to `scMOG_prepare_data.ipynb` for detailed instructions.

## Run model
scMOG is run using the command line. To train the model with custom datasets, you need to make the following modifications to the `train_protein.py` file:
 ```python
rna_path = dataset_path + '/CITE_train_10k_RNA'
modal_path = dataset_path + '/CITE_train_10k_protein'
rna_test_path = dataset_path + '/CITE_test_5k_RNA'
modal_test_path = dataset_path + '/CITE_test_5k_protein'
```
replace with
```python
rna_path = dataset_path + '/train_gene_expression.csv'
modal_path = dataset_path + '/train_protein_expression.csv'
rna_test_path = dataset_path + '/test_gene_expression.csv'
modal_test_path = dataset_path + '/test_protein_expression.csv'
```
Please make sure to set your working directory to `your customized path/code/scMOG/scMOG-main/scMOG_code/bin` first, then use the command line to train the model:
```python
python train_protein.py --outdir "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMOG" --dataset_path "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMOG"
```
The `outdir` parameter specifies the path where the trained model will be saved. The `dataset_path` parameter is the path to the input data.
Once training is complete, you can perform inference on the test set by running the following command:
```python
python "your customized path/code/scMOG/scMOG-main/scMOG_code/bin/predict-protein.py" --outdir "your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMOG"
```
Note that to correctly load the trained model, you need to set the working directory to `your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMOG`. The `outdir` parameter specifies the path where the imputed protein expression data for the test set will be saved.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `scMOG_evaluation.ipynb` for detailed instructions.



