# -*- coding: utf-8 -*-

## Stand-alone version of DEICODE (https://github.com/biocore/DEICODE)

from biom import load_table
from deicode.rpca import rpca
from deicode.rpca import auto_rpca
import numpy as np
import pandas as pd
import os
from skbio import DistanceMatrix
from deicode.preprocessing import rclr

from os import chdir, getcwd
wd=getcwd()
chdir(wd)

os.chdir('/Users/johndoe/TransplantLines/mock_data/decoide_input/')

# This is the abundance table in .biom format
biom = load_table('Mocks_Taxa_Filtered.biom')
# print(biom.shape)
# This is the abundance table in .csv format
table = pd.read_csv('Mocks_Taxa_Filtered.csv', index_col=0)
# print(table.shape)
# print(table)

# We use the auto_rpca function to run DECOIDE; it finds the number of dimensions that most parsimoniously represents the data
# No filtering is applied here as this has already been done on the input data
ordination, distance = auto_rpca(biom, 
								 min_feature_count = 0, 
                                 min_sample_count = 0, 
                                 min_feature_frequency = 0)

# Output ordination
# the sample loadings
spca_df = ordination.samples
# print(spca_df)
pd.DataFrame(spca_df).to_csv("/Users/johndoe/TransplantLines/mock_data/decoide_output/Mocks_Taxa_Filtered_sample_loadings.csv")
# the feature loadings
fpca_df = ordination.features
# print(fpca_df)
pd.DataFrame(fpca_df).to_csv("/Users/johndoe/TransplantLines/mock_data/decoide_output/Mocks_Taxa_Filtered_feature_loadings.csv")


# To extract the RCLR matrix, we simply run the rclr function on the table
table = table.as_matrix()
# print(table)
# print(table)
# This runs the rclr function (clr on non-zero values)
table_rclr = rclr(table)
# print(table_rclr)
# print(table_rclr.shape)
pd.DataFrame(table_rclr).to_csv("/Users/johndoe/TransplantLines/mock_data/decoide_output/Mocks_Taxa_Filtered_rclr.csv")

# Below are two outputs that we don't use, but that DECOIDE can produce. 

# Output the Aitchison distance matrix
aitchison = distance.to_data_frame()
# print(aitchison.shape)
pd.DataFrame(aitchison).to_csv("/Users/johndoe/TransplantLines/mock_data/decoide_output/Mocks_Taxa_Filtered_aitchison.csv")

# USV^T - imputed CLR
USVT = ordination.samples.values @ np.diag(ordination.eigvals.values) @ ordination.features.values.T
USVT = USVT - USVT.mean(axis=0)
USVT = USVT - USVT.mean(axis=1).reshape(-1, 1)
# print(USVT.T.shape)
pd.DataFrame(USVT.T).to_csv("/Users/johndoe/TransplantLines/mock_data/decoide_output/Mocks_Taxa_Filtered_usvtt.csv")


