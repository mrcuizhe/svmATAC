#!/usr/bin/env python
# coding: utf-8

# ####  Generate Barcode

# In[1]:


import itertools
import random
import pandas as pd
import numpy as np
import os


list_barcodes = [''.join(x) for x in itertools.product('ATCG', repeat=8)]

metadata = pd.read_csv('../input/metadata_corces2016_sorted.tsv',sep='\t',index_col = 0)
metadata.insert(0,'ID',metadata.index)

random.seed(2019)
barcodes = random.sample(list_barcodes,metadata.shape[0])

# In[8]:
metadata.insert(0,'barcode', barcodes)

metadata.to_csv('../tmp/corces2016_barcode_metadata.tsv',sep='\t',header=True,index=False)


# #### Prepend barcode to each cell

# `./prepend_barcodes.sh`

# #### Merge all the bam files

path = '../tmp/sc-bams_barcodes/'
files = [os.path.join(path,f) for f in os.listdir(path) if f.endswith(".bam")]

pd.DataFrame(files).to_csv('../tmp/list_bamfiles.txt',sep='\t',header=False,index=False)




