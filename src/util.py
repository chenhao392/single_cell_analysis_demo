#!/usr/bin/python

import numpy as np

#functions copied from scVI official tutorial

def get_score(normalized_adata, gene_set):
    score = np.zeros(normalized_adata.n_obs)
    for gene in gene_set['positive']:
        expression = np.array(normalized_adata[:, gene].X)
        score += expression.flatten()
    for gene in gene_set['negative']:
        expression = np.array(normalized_adata[:, gene].X)
        score -= expression.flatten()
    return score

def get_cell_mask(normalized_adata, gene_set):
    score = get_score(normalized_adata, gene_set)
    cell_idx = score.argsort()[-50:]
    mask = np.zeros(normalized_adata.n_obs)
    mask[cell_idx] = 1
    return mask.astype(bool)
