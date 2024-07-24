# -*- coding: utf-8 -*-
import torch.nn as nn
import torch
from torch.utils.data import DataLoader
import torch.utils.data as Data
from network import Shared_Generator, Protein_Predictor
import para


def _init_weight(m):
    if type(m) == nn.Linear:
        torch.nn.init.xavier_uniform_(m.weight)
        m.bias.data.fill_(0.01)
    elif type(m) == nn.BatchNorm1d:
        torch.nn.init.normal_(m.weight.data, 1.0, 0.02)
        torch.nn.init.constant_(m.bias.data, 0.0)
        

def training(source_rna, source_protein, protein_name):
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    n_genes = source_rna.shape[1]
    
    
    shared_generator = Shared_Generator(input_dim = n_genes, 
                                latent_dim1 = para.latent_dims[0], 
                                latent_dim2 = para.latent_dims[1]).to(device)

    protein_predictor = Protein_Predictor(latent_dim2  = para.latent_dims[1], 
                                protein_latent_dim = para.latent_dims[2], 
                                protein_name = protein_name).to(device)


    shared_generator.apply(_init_weight)
    protein_predictor.apply(_init_weight)
    
    
    optimizer_sg = torch.optim.Adam(shared_generator.parameters(), lr= para.lr)
    optimizer_pp = torch.optim.Adam(protein_predictor.parameters(), lr= para.lr)

       
    xx = torch.tensor(source_rna).to(device)
    yy = torch.tensor(source_protein).to(device)
    tor_dataset = Data.TensorDataset(xx,yy)
    dataloader = DataLoader(tor_dataset, batch_size = para.batch_size)
    
    
    for epoch in range(para.n_epochs):
    
        for i, (batch_rna, batch_protein) in enumerate(dataloader):
    
            optimizer_sg.zero_grad()
            optimizer_pp.zero_grad()
            
            latent_dim_rep = shared_generator(batch_rna.float().to(device))
            prediction = protein_predictor(latent_dim_rep.float().to(device))
            
            protein_prediction = None
            for p in protein_name:
                temp = prediction[p]
                if protein_prediction is None:
                    protein_prediction = temp
                else:
                    protein_prediction = torch.cat((protein_prediction, temp), 1)
    
            
            #L1 loss
            loss_fun = torch.nn.L1Loss()
            prediction_loss = loss_fun(batch_protein.float(), protein_prediction)
            
            total_loss = prediction_loss
           
            
            total_loss.backward()
            
            optimizer_sg.step()
            optimizer_pp.step()
    
            print(
                "[gener Epoch %d/%d] [gener Batch %d/%d] [gener L1 loss: %f] " 
                % (epoch, para.n_epochs, i, len(dataloader), prediction_loss.item())
            )
            
    return shared_generator, protein_predictor


def predicting(shared_generator, protein_predictor, target_rna, protein_name):
    
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    shared_generator.eval()
    protein_predictor.eval()
    
    xx = torch.tensor(target_rna)
    
    
    with torch.no_grad():
        latent_dim_rep = shared_generator(xx.float().to(device))
        prediction = protein_predictor(latent_dim_rep.float().to(device))
        
        protein_prediction = None
        for p in protein_name:
            temp = prediction[p]
            if protein_prediction is None:
                protein_prediction = temp
            else:
                protein_prediction = torch.cat((protein_prediction, temp), 1)
        
        protein_prediction = protein_prediction.data.detach().cpu().numpy()
    
    return protein_prediction



