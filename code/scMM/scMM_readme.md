## Installation
We use scMM based on the dance package for Python. Please refer to the official GitHub repository <https://github.com/OmicsML/dance> for instructions on how to install the dance package.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Prepare data
To accommodate the input requirements of scMM, we need to convert gene expression data and protein expression data to anndata format separately. Please refer to `scMM_prepare_data.ipynb` for detailed instructions.

## Run model
scMM is run using the command line. Please make sure to set your working directory to `your customized path/code/scMM` first. Use the command line to train the model:
```python
python scMM.py --subtask openproblems_bmmc_cite_phase2_rna --device cuda --output_path "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMM"
```
The `output_path` parameter specifies the relative path within the `outputs` folder of the parent directory where the prepared dataset from the previous step is saved. The results from this step will also be stored in this folder. For example, if the prepared dataset is saved in `your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMM/openproblems_bmmc_cite_phase2_rna`, then the `output_path` parameter should be `different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/scMM`. After running model, the imputed protein expression data for the test set, the ground truth protein expression data for the test set after preprocessing, and the trained model will be saved.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `scMM_evaluation.ipynb` for detailed instructions.







