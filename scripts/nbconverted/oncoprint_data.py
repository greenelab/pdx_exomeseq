
# coding: utf-8

# # Generating OncoPrint Data Files
# 
# The script will process all variant files and output files in an ingestible format for the R OncoPrint function.

# In[1]:

import os
import pandas as pd


# In[2]:

cosmic_file = os.path.join('results', 'all_common_replicate_COSMIC_variants.tsv')
cosmic_df = pd.read_table(cosmic_file, index_col=0)
print(cosmic_df.shape)
cosmic_df.head()


# In[3]:

# What are the 50 most commonly altered genes?
top_n = 50
paad_genes = cosmic_df['Gene.refGene'].value_counts().head(top_n).index.tolist()
print(paad_genes)


# ## Generate OncoPrint Data
# 
# ### For All Replicates

# In[4]:

get_ipython().run_cell_magic('time', '', "variant_file_path = os.path.join('results', 'processed_vcfs')\nvariant_assign = []\ncase_id = []\nall_cosmic_ids = []\nfor variant_file in os.listdir(variant_file_path):\n    # Load and subset file to only variants in the COSMIC db\n    variant_df = pd.read_table(os.path.join(variant_file_path, variant_file), index_col=0)\n    variant_sub_df = variant_df[variant_df['cosmic70'] != '.']\n    \n    # Build a list of all COSMIC IDs for separate R visualization (similarity heatmaps)\n    all_cosmic_ids += variant_sub_df['cosmic70'].tolist()\n    \n    # Define mutated genes if they exist for the given variant\n    variant_class = ['MUT;' if x in variant_sub_df['Gene.refGene'].tolist() else ''\n                     for x in paad_genes]\n    \n    # Store results\n    variant_assign.append(variant_class)\n    case_id.append(variant_file.replace('_001_processed_variants.tsv.bz2', ''))")


# In[5]:

# Generate and save oncoprint data for all replicates
oncoprint_file = os.path.join('results', 'oncoprint_replicates.tsv')
oncoprint_df = pd.DataFrame(variant_assign, index=case_id, columns=paad_genes)
oncoprint_df.index.name = 'Case.ID'
oncoprint_df.to_csv(oncoprint_file, sep='\t')


# ### Consensus samples

# In[6]:

# Generate oncoprint data for consensus samples (COSMIC variant exists in all replicates)
full_variant_file = os.path.join('results', 'all_common_replicate_COSMIC_variants.tsv')
full_variant_df = pd.read_table(full_variant_file, index_col=0)
full_variant_df.head(2)


# In[7]:

variant_assign_consensus = []
case_id_consensus = []
for sample_id in set(full_variant_df['sample_id']):
    # Subset file to given sample ID
    variant_sub_df = full_variant_df[full_variant_df['sample_id'] == sample_id]
    
    # Define mutated genes if they exist for the given variant
    variant_class = ['MUT;' if x in variant_sub_df['Gene.refGene'].tolist() else ''
                     for x in paad_genes]
    
    # Store results
    variant_assign_consensus.append(variant_class)
    case_id_consensus.append(sample_id)


# In[8]:

# Generate and save oncoprint data for consensus samples
oncoprint_consensus_file = os.path.join('results', 'oncoprint_consensus.tsv')
oncoprint_consensus_df = pd.DataFrame(variant_assign_consensus, index=case_id_consensus,
                                      columns=paad_genes)
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
print('All COSMIC mutations: {}'.format(len(all_cosmic_ids)))
unique_cosmic_ids = set(all_cosmic_ids)
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
cosmic_common_df = pd.DataFrame(cosmic_similarity_list, index=case_id, columns=unique_cosmic_ids)
cosmic_common_df.index.name = 'Case.ID'
cosmic_common_df.to_csv(cosmic_common_file, sep='\t')


# ### Consensus samples

# In[13]:

case_id_consensus = []
cosmic_similarity_consensus_list = []
for sample_id in set(full_variant_df['sample_id']):
    # Subset file to given sample ID
    variant_sub_df = full_variant_df[full_variant_df['sample_id'] == sample_id]
    
    # Define membership in COSMIC IDs
    cosmic_class = [1 if x in variant_sub_df['cosmic70'].tolist() else 0
                    for x in unique_cosmic_ids]
    
    # Store results
    cosmic_similarity_consensus_list.append(cosmic_class)
    case_id_consensus.append(sample_id)


# In[14]:

common_cosmic_consensus_file = os.path.join('results', 'cosmic_similarity_consensus.tsv')
cosmic_common_consensus_df = pd.DataFrame(cosmic_similarity_consensus_list,
                                          index=case_id_consensus,
                                          columns=unique_cosmic_ids)
cosmic_common_consensus_df.index.name = 'Case.ID'
cosmic_common_consensus_df.to_csv(common_cosmic_consensus_file, sep='\t')

