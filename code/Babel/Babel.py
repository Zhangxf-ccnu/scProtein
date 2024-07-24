# -*- coding: utf-8 -*-
import argparse
import logging
import os

import pandas as pd
import torch

from dance import logger
from dance.datasets.multimodality import ModalityPredictionDataset
from dance.modules.multi_modality.predict_modality.babel import BabelWrapper
from dance.utils import set_seed
import anndata as ad


if __name__ == "__main__":
    OPTIMIZER_DICT = {
        "adam": torch.optim.Adam,
        "rmsprop": torch.optim.RMSprop,
    }
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--subtask", default="openproblems_bmmc_cite_phase2_rna")
    parser.add_argument("-device", "--device", default="cuda")
    parser.add_argument("-cpu", "--cpus", default=1, type=int)
    parser.add_argument("-seed", "--seed", default=1, type=int)
    parser.add_argument("--runs", type=int, default=1, help="Number of repetitions")
    parser.add_argument("-m", "--model_folder", default="./models")
    parser.add_argument("--outdir", "-o", default="./logs", help="Directory to output to")
    parser.add_argument("--lossweight", type=float, default=1., help="Relative loss weight")
    parser.add_argument("--lr", "-l", type=float, default=0.01, help="Learning rate")
    parser.add_argument("--batchsize", "-b", type=int, default=64, help="Batch size")
    parser.add_argument("--hidden", type=int, default=64, help="Hidden dimensions")
    parser.add_argument("--earlystop", type=int, default=2, help="Early stopping after N epochs")
    parser.add_argument("--naive", "-n", action="store_true", help="Use a naive model instead of lego model")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max_epochs", type=int, default=40)
    args = parser.parse_args()
    args.resume = True
    
    args.model_folder = os.path.join(os.path.abspath(os.path.join(os.getcwd(), '../../outputs')), args.model_folder)
    args.outdir = os.path.join(os.path.abspath(os.path.join(os.getcwd(), '../../outputs')), args.outdir)
    
    torch.set_num_threads(args.cpus)
    rndseed = args.seed
    set_seed(rndseed)
    
    dataset = ModalityPredictionDataset(args.subtask, root=args.outdir)
    dataset._maybe_download()
    modalities = []
    for mod_path in dataset.data_paths:
        logger.info(f"Loading {mod_path}")
        modalities.append(ad.read_h5ad(mod_path))
    raw_data = modalities
    data = dataset._raw_to_dance(raw_data)
    
    data.set_config(feature_mod="mod1", label_mod="mod2", feature_channel_type="X", label_channel_type="X")
    
    device = args.device
    args.outdir = os.path.abspath(args.outdir)
    os.makedirs(args.model_folder, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)
    
    fh = logging.FileHandler(f"{args.outdir}/training_{args.subtask}_{args.seed}.log", "w")
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)
    
    for arg in vars(args):
        logger.info(f"Parameter {arg}: {getattr(args, arg)}")
    
    x_train, y_train = data.get_train_data(return_type="torch")
    x_test, y_test = data.get_test_data(return_type="torch")
    x_train, y_train, x_test, y_test = x_train.float(), y_train.float(), x_test.float(), y_test.float()
    
    set_seed(args.seed)
    model = BabelWrapper(args, dim_in=x_train.shape[1], dim_out=y_train.shape[1])
    model.fit(x_train, y_train, val_ratio=0.15)
    
    predicted_protein_expression = pd.DataFrame(model.predict(x_test).cpu().numpy(), columns=raw_data[3].var.index)
    predicted_protein_expression.to_csv(os.path.join(args.outdir, "test_protein_prediction.txt"), sep="\t")
    
    true_protein_expression = pd.DataFrame(y_test, columns=raw_data[3].var.index)
    true_protein_expression.to_csv(os.path.join(args.outdir, "test_protein_groundtruth.txt"), sep="\t")
    
    torch.save(model, os.path.join(args.model_folder, "model.pth"))     







