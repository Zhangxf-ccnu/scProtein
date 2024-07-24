## Installation
We use Babel based on the dance package for Python. Please refer to the official GitHub repository <https://github.com/OmicsML/dance> for instructions on how to install the dance package.

After installation, we will demonstrate how to run this method using the CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example.

## Preprocess data
We first preprocess gene expression data and protein expression data. Additionally, to accommodate the input requirements of Babel, we need to convert gene expression data and protein expression data to anndata format separately. Please refer to `Babel_preprocess.ipynb` for detailed instructions.

## Run model
Babel is run using the command line. Please make sure to set your working directory to `your customized path/code/Babel` first.
Use the command line to train the model:
```python
python Babel.py --subtask openproblems_bmmc_cite_phase2_rna --device cuda --outdir "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/Babel" --model_folder "different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/Babel"
```
The `outdir` parameter specifies the relative path within the `outputs` folder of the parent directory where the preprocessed dataset from the previous step is saved. For example, if the preprocessed dataset is saved in `your customized path/outputs/different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/Babel/openproblems_bmmc_cite_phase2_rna`, then the `outdir` parameter should be `different samples/CITE-SLN111-Gayoso-Mouse1toMouse2/Babel`.
After running model, the imputed protein expression data for the test set, the ground truth protein expression data for the test set, and the log file will be saved in this folder. The `model_folder` parameter is also a relative path within the `outputs` folder, used to store the trained models.

## Evaluation
The prediction can be evaluated using PCC and RMSE at the protein and cell levels. Please refer to `Babel_evaluation.ipynb` for detailed instructions.