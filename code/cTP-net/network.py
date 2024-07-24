# -*- coding: utf-8 -*-
import torch.nn as nn


class Shared_Generator(nn.Module):
    def __init__(self, input_dim, latent_dim1, latent_dim2):
        super(Shared_Generator, self).__init__()
        
        self.fc1 = nn.Sequential(
                  nn.Linear(input_dim, latent_dim1),
                  nn.ReLU(), 
                  nn.Linear(latent_dim1, latent_dim2),
                  nn.ReLU()
                  )

    def forward(self, x):
        out01 = self.fc1(x)
        return out01


class Protein_Predictor(nn.Module):
    def __init__(self, latent_dim2, protein_latent_dim, protein_name):
        super(Protein_Predictor, self).__init__()
        
        self.fc2 = nn.ModuleDict({})
        for p in protein_name:
           self.fc2[p] = nn.Sequential(
               nn.Linear(latent_dim2, protein_latent_dim),
               nn.ReLU()
               )
           
        self.fc3 = nn.ModuleDict({})
        for p in protein_name:
           self.fc3[p] = nn.Sequential(
               nn.Linear(protein_latent_dim, 1), 
               )

        self.protein_name = protein_name
        
    def forward(self, out01):
        out02 = {}
        for p in self.protein_name:
            out02[p] = self.fc3[p](self.fc2[p](out01))
        return out02











