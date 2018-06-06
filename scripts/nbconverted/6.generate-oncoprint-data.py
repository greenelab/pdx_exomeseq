
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

# What are the 50 most commonly altered COSMIC genes?
top_n = 50
paad_genes = cosmic_all_df['Gene.refGene'].value_counts().head(top_n).index.tolist()
cosmic_all_df['Gene.refGene'].value_counts().head(20)


# ### Define Functions for Oncoprint Data Processing

# In[3]:


def process_variants(variant_dir, focus_variants, strip_text, process_cosmic=False,
                     id_updater=None):
    """
    Retrieve VCF files from an input directory and determine membership
    
    Arguments:
    variant_dir - the directory to search for variant files to load
    focus_variants - a list of genes or variants to search for in samples
    strip_text - a string of text to strip from variant files
    process_cosmic - boolean to determine if cosmic variants are to be processed
    id_updater - a dictionary of sample ID mappings (defaults to None)
    
    Output:
    A dataframe that is ready for input into oncoprint function
    """

    variant_assign = []
    case_ids = []
    for variant_file in os.listdir(variant_dir):
        # Load and subset file to only variants in the COSMIC db
        variant_df = pd.read_table(os.path.join(variant_dir, variant_file), index_col=0)
        variant_sub_df = variant_df[variant_df['cosmic70'] != '.']

        # Define mutated genes or variants if they exist for the given sample
        if process_cosmic:
            variant_class = [1 if x in variant_sub_df['cosmic70'].tolist() else 0
                             for x in focus_variants]
        else:
            variant_class = ['MUT;' if x in variant_sub_df['Gene.refGene'].tolist() else ''
                             for x in focus_variants]

        # Store results
        sample_id = variant_file.replace(strip_text, '')
        variant_assign.append(variant_class)
        if id_updater is not None:
            sample_id = variant_file.replace(variant_file.split('_')[0],
                                             id_updater[variant_file.split('_')[0]])
        case_ids.append(sample_id)
    
    # Process output variants
    output_df = pd.DataFrame(variant_assign,
                             index=case_ids,
                             columns=focus_variants).sort_index()
    output_df.index.name = 'Case.ID'

    return output_df


# ## Generate OncoPrint Data
# 
# ### For All Replicates

# In[4]:


# Process each replicate by observed COSMIC mutation
replicate_file_path = os.path.join('results', 'processed_vcfs')
replicate_strip_text = '_001_processed_variants.tsv.bz2'

replicate_oncoprint_df = process_variants(variant_dir=replicate_file_path,
                                          focus_variants=paad_genes,
                                          strip_text=replicate_strip_text,
                                          process_cosmic=False,
                                          id_updater=None)

# Output file
replicate_output_file = os.path.join('results', 'oncoprint_replicates.tsv')
replicate_oncoprint_df.to_csv(replicate_output_file, sep='\t')


# ### For Merged Samples

# In[5]:


# Process each replicate by observed COSMIC mutation
merged_file_path = os.path.join('results', 'processed_merged_vcfs')
merged_strip_text = '_processed_variants.tsv.bz2'

merged_oncoprint_df = process_variants(variant_dir=merged_file_path,
                                       focus_variants=paad_genes,
                                       strip_text=merged_strip_text,
                                       process_cosmic=False,
                                       id_updater=None)

# Output file
merged_output_file = os.path.join('results', 'oncoprint_merged.tsv')
merged_oncoprint_df.to_csv(merged_output_file, sep='\t')


# ## COSMIC Mutational Similarity
# 
# Output mutational similarity data for all replicates and consensus samples.
# The COSMIC mutational similarity is built from a (0,1) sample by COSMIC mutation matrix.

# In[6]:


# How many COSMIC mutation IDs are in the entire set and how many are unique?
print('All COSMIC mutations: {}'.format(cosmic_all_df.shape[0]))
unique_cosmic_ids = set(cosmic_all_df['cosmic70'])
print('Unique COSMIC mutations: {}'.format(len(unique_cosmic_ids)))


# ### For All Replicates

# In[7]:


# Obtain replicate cosmic similarity matrix
replicate_cosmic_df = process_variants(variant_dir=replicate_file_path,
                                       focus_variants=unique_cosmic_ids,
                                       strip_text=replicate_strip_text,
                                       process_cosmic=True,
                                       id_updater=None)

replicate_common_file = os.path.join('results', 'cosmic_similarity_replicates.tsv')
replicate_cosmic_df.to_csv(replicate_common_file, sep='\t')


# ### Consensus samples

# In[8]:


# Obtain consensus cosmic similarity matrix
merged_cosmic_df = process_variants(variant_dir=merged_file_path,
                                       focus_variants=unique_cosmic_ids,
                                       strip_text=merged_strip_text,
                                       process_cosmic=True,
                                       id_updater=None)

merged_common_file = os.path.join('results', 'cosmic_similarity_merged.tsv')
merged_cosmic_df.to_csv(merged_common_file, sep='\t')


# ## What about prefiltered variants (i.e. before COSMIC filtering)
# 
# Observed merged samples with cosmic similarity only

# In[11]:


# Load all prefiltered cosmic variants called in this dataset
# This file was generated in filter_variants.ipynb
file = os.path.join('results', 'all_cosmic_prefiltered_variants.tsv')
cosmic_prefiltered_df = pd.read_table(file)
cosmic_prefiltered_df.head(3)


# In[12]:


paad_prefiltered_genes = cosmic_prefiltered_df['Gene.refGene'].value_counts().head(top_n).index.tolist()
cosmic_prefiltered_df['Gene.refGene'].value_counts().head(20)


# In[13]:


# Only consider mutations that change amino acid sequence
cosmic_prefiltered_df = cosmic_prefiltered_df[cosmic_prefiltered_df['AAChange.refGene'] != '.']


# In[14]:


# How many prefiltered COSMIC mutation IDs are in the entire set and how many are unique?
print('All prefiltered COSMIC mutations: {}'.format(cosmic_prefiltered_df.shape[0]))
unique_mutations = set(cosmic_prefiltered_df['AAChange.refGene'])
print('Unique nonsilent mutations: {}'.format(len(unique_mutations)))


# In[15]:


case_id_consensus = []
cosmic_similarity_consensus_list = []
for sample_id in set(cosmic_prefiltered_df['final_id']):
    # Subset file to given sample ID
    variant_sub_df = cosmic_prefiltered_df.query('final_id == @sample_id')
    
    # Define membership in COSMIC IDs
    cosmic_class = [1 if x in variant_sub_df['AAChange.refGene'].tolist() else 0
                    for x in unique_mutations]
    
    # Store results
    cosmic_similarity_consensus_list.append(cosmic_class)
    case_id_consensus.append(sample_id)


# In[16]:


common_cosmic_consensus_file = os.path.join('results', 'cosmic_prefiltered_similarity_merged.tsv')
cosmic_common_consensus_df = pd.DataFrame(cosmic_similarity_consensus_list,
                                          index=case_id_consensus,
                                          columns=unique_mutations)
cosmic_common_consensus_df.index.name = 'Case.ID'
cosmic_common_consensus_df.to_csv(common_cosmic_consensus_file, sep='\t')

