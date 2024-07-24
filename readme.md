# Benchmarking single-cell cross-omics imputation methods for surface protein expression
## Overview
In this study, we present an extensive benchmark of twelve state-of-the-art imputation methods using eleven CITE-seq and REAP-seq datasets across six distinct scenarios. We employ various accuracy measures to quantitatively evaluate the predicted values at both the protein and cell levels. Additionally, we assess the methods’ sensitivity to training data size, robustness across experiments, and efficiency in terms of time and memory. We also compare their popularity based on the number of stars on their official GitHub repositories and evaluate the quality of their documentation, code, and tutorials.

![image](Overview.png 'pipeline')

## Datasets
The all datasets used in the manuscript have been curated and uploaded to a public repository in Zenodo with a DOI assignment (DOI: <https://doi.org/10.5281/zenodo.12725699>). The datasets are stored according to the experiments. Each experiment folder contains datasets that have undergone quality control and an intersection of genes and proteins, making them directly usable as input data for the experiments. If you find our datasets useful in your research, please consider citing:
```bibtex
@article{li2024benchmark,
  title={Benchmarking single-cell cross-omics imputation methods for surface protein expression},
  author={Li Chen-yang, Hong Yong-jia, Li Bo and Zhang Xiao-fei},
  url={https://doi.org/10.5281/zenodo.12725699},
  year={2024},
}
```

## Implementation
After downloading the dataset and code, please place them according to the following directory structure:
```
root directory
│
├── datasets
│   ├── different clinical states
│   ├── different protocols
│   ├── ...
│   └── random holdout&different training data sizes
│
├── code
│   ├── Babel
│   ├── cTP-net
│   ├── ...
│   └── TotalVI
│
└── outputs
    ├── different clinical states
    ├── different protocols
    ├── ...
    └── random holdout&different training data sizes
```
In `code` folder, we provide detailed tutorials and documentation for 12 methods. The implementation pipeline of each method is illustrated using CITE-SLN111-Gayoso: Mouse1→Mouse2 experiment in different samples scenario as example. For methods that can handle batch information within the dataset, we also provide CITE-PBMC-Li: Group1→Group2 experiment in different samples scenario as example. Note that, unless otherwise specified, the working directory for each program is `root directory`. The `outputs` folder is used to store the results of experiments.

## Contact
Please do not hesitate to contact Mr. Chen-yang Li (<lichy@mails.ccnu.edu.cn>) or Dr. Xiao-Fei Zhang (<zhangxf@ccnu.edu.cn>) to seek any clarifications regarding any contents or operation of the archive.






