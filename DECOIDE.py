# -*- coding: utf-8 -*-

## Stand-alone version of DECOIDE

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

os.chdir('/Users/johndoe/decoide/')

# This is the abundance table in .biom format
biom = load_table('feature_tab.biom')
# print(biom.shape)
# This is the abundance table in .csv format
table = pd.read_csv('feature_tab.csv', index_col=0)
# print(table.shape)
# print(table)

# This runs DECOIDE
# no filtering is applied as this has already been done on the feature_tab
ordination, distance = auto_rpca(biom, 
								 min_feature_count = 0, 
                                 min_sample_count = 0, 
                                 min_feature_frequency = 0)
# Output the Aitchison distance matrix
# print(distance.shape)
aitchison = distance.to_data_frame()
# print(aitchison.shape)
pd.DataFrame(aitchison).to_csv("/Users/johndoe/decoide/aitchison.csv")

# Output ordination
# the sample loadings
spca_df = ordination.samples
# print(spca_df)
pd.DataFrame(spca_df).to_csv("/Users/johndoe/decoide/sample_loadings.csv")
# the feature loadings
fpca_df = ordination.features
# print(fpca_df)
pd.DataFrame(fpca_df).to_csv("/Users/johndoe/decoide/feature_loadings.csv")


# To extract the RCLR matrix, we simply run the rclr function on the table
table = table.as_matrix()
# print(table)
# print(table)
# This runs the rclr function (clr on non-zero values)
table_rclr = rclr(table)
# print(table_rclr)
# print(table_rclr.shape)
pd.DataFrame(table_rclr).to_csv("/Users/johndoe/decoide/decoide/rclr.csv")

# USV^T - imputed CLR
print(ordination.samples.shape)
print(ordination.features.shape)
print(np.diag(ordination.eigvals).shape)

USVT = ordination.samples.values @ np.diag(ordination.eigvals.values) @ ordination.features.values.T
USVT = USVT - USVT.mean(axis=0)
USVT = USVT - USVT.mean(axis=1).reshape(-1, 1)
# print(USVT.T.shape)
pd.DataFrame(USVT.T).to_csv("/Users/userX/decoide/usvtt.csv")


