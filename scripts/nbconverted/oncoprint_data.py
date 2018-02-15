
# coding: utf-8

# # Generating OncoPrint Data Files
# 
# The script will process all variant files and output files in an ingestible format for the R OncoPrint function.
# 
# It will output oncoprint data for both replicate files and the merged variant callsets.

# In[1]:

import os
import pandas as pd


# In[2]:

# Load all cosmic variants called in this dataset
# This file was generated in filter_variants.ipynb
cosmic_all_file = os.path.join('results', 'all_cosmic_variants.tsv')
cosmic_all_df = pd.read_table(cosmic_all_file)
cosmic_all_df.head(3)


# In[3]:

# What are the 50 most commonly altered COSMIC genes?
top_n = 50
paad_genes = cosmic_all_df['Gene.refGene'].value_counts().head(top_n).index.tolist()
cosmic_all_df['Gene.refGene'].value_counts().head(20)


# ## Generate OncoPrint Data
# 
# ### For All Replicates

# In[4]:

get_ipython().run_cell_magic('time', '', "\n# Process each replicate by observed COSMIC mutation\nvariant_file_path = os.path.join('results', 'processed_vcfs')\nvariant_assign = []\ncase_id = []\nfor variant_file in os.listdir(variant_file_path):\n    # Load and subset file to only variants in the COSMIC db\n    variant_df = pd.read_table(os.path.join(variant_file_path, variant_file), index_col=0)\n    variant_sub_df = variant_df[variant_df['cosmic70'] != '.']\n    \n    # Define mutated genes if they exist for the given variant\n    variant_class = ['MUT;' if x in variant_sub_df['Gene.refGene'].tolist() else ''\n                     for x in paad_genes]\n    \n    # Store results\n    variant_assign.append(variant_class)\n    case_id.append(variant_file.replace('_001_processed_variants.tsv.bz2', ''))")


# In[5]:

# Generate and save oncoprint data for all replicates
oncoprint_file = os.path.join('results', 'oncoprint_replicates.tsv')
oncoprint_df = pd.DataFrame(variant_assign, index=case_id, columns=paad_genes)
oncoprint_df.index.name = 'Case.ID'
oncoprint_df.to_csv(oncoprint_file, sep='\t')


# In[6]:

oncoprint_df.head(3)


# ### Merged Samples

# In[7]:

variant_assign_consensus = []
case_id_consensus = []
for sample_id in set(cosmic_all_df['sample_name']):
    # Subset file to given sample ID
    variant_sub_df = cosmic_all_df.query('sample_name == @sample_id')
    
    # Define mutated genes if they exist for the given variant
    variant_class = ['MUT;' if x in variant_sub_df['Gene.refGene'].tolist() else ''
                     for x in paad_genes]
    
    # Store results
    variant_assign_consensus.append(variant_class)
    case_id_consensus.append(sample_id)


# In[8]:

# Generate and save oncoprint data for consensus samples
oncoprint_consensus_file = os.path.join('results', 'oncoprint_merged.tsv')

oncoprint_consensus_df = (
    pd.DataFrame(variant_assign_consensus,
                 index=case_id_consensus,
                 columns=paad_genes)
    )
oncoprint_consensus_df.index.name = 'Case.ID'
oncoprint_consensus_df.to_csv(oncoprint_consensus_file, sep='\t')


# In[9]:

oncoprint_consensus_df.head(3)


# ## COSMIC Mutational Similarity
# 
# Output mutational similarity data for all replicates and consensus samples. The COSMIC mutational similarity is built from a (0,1) sample by COSMIC mutation matrix.
# 
# ### For All Replicates

# In[10]:

# How many COSMIC mutation IDs are in the entire set and how many are unique?
print('All COSMIC mutations: {}'.format(cosmic_all_df.shape[0]))
unique_cosmic_ids = set(cosmic_all_df['cosmic70'])
print('Unique COSMIC mutations: {}'.format(len(unique_cosmic_ids)))


# In[11]:

case_id = []
cosmic_similarity_list = []
for variant_file in os.listdir(variant_file_path):
    # Load and subset file to only variants in the COSMIC db
    variant_df = pd.read_table(os.path.join(variant_file_path, variant_file), index_col=0)
    variant_sub_df = variant_df[variant_df['cosmic70'] != '.']
    
    # Define membership in COSMIC IDs
    cosmic_class = [1 if x in variant_sub_df['cosmic70'].tolist() else 0
                    for x in unique_cosmic_ids]
    
    # Store results
    cosmic_similarity_list.append(cosmic_class)
    case_id.append(variant_file.replace('_001_processed_variants.tsv.bz2', ''))


# In[12]:

# Generate COSMIC id membership data (for downstream similarity matrix)
cosmic_common_file = os.path.join('results', 'cosmic_similarity_replicates.tsv')
cosmic_common_df = (
    pd.DataFrame(cosmic_similarity_list,
                 index=case_id,
                 columns=unique_cosmic_ids)
    )
cosmic_common_df.index.name = 'Case.ID'
cosmic_common_df.to_csv(cosmic_common_file, sep='\t')


# ### Consensus samples

# In[13]:

case_id_consensus = []
cosmic_similarity_consensus_list = []
for sample_id in set(cosmic_all_df['sample_name']):
    # Subset file to given sample ID
    variant_sub_df = cosmic_all_df.query('sample_name == @sample_id')
    
    # Define membership in COSMIC IDs
    cosmic_class = [1 if x in variant_sub_df['cosmic70'].tolist() else 0
                    for x in unique_cosmic_ids]
    
    # Store results
    cosmic_similarity_consensus_list.append(cosmic_class)
    case_id_consensus.append(sample_id)


# In[14]:

common_cosmic_consensus_file = os.path.join('results', 'cosmic_similarity_merged.tsv')
cosmic_common_consensus_df = pd.DataFrame(cosmic_similarity_consensus_list,
                                          index=case_id_consensus,
                                          columns=unique_cosmic_ids)
cosmic_common_consensus_df.index.name = 'Case.ID'
cosmic_common_consensus_df.to_csv(common_cosmic_consensus_file, sep='\t')

