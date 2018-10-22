
# coding: utf-8

# # Determining the difference in variant calling in human-only samples `004` and `005`
# 
# **Gregory Way 2018**
# 
# Samples `004` and `005` are human tumors.
# They were previously included in the entire `disambiguate` pipeline, where the WES reads were aligned to both human and mouse genomes.
# 
# In the pipeline, all WES reads are aligned to both genomes.
# Reads are then "disambiguated" to either species in which the species with the highest alignment score per read is assigned that given read.
# There was some interesting variants observed in the human only samples that appeared to belong only to mouse.
# We hypothesized that the reason we are observing these variants is because of an error in the disambiguation step.
# 
# To determine if this was the case, I also aligned the two samples to the human genome only, and called variants in these files.
# 
# The following script compares the variants called in `004` and `005` between the `disambiguate` pipeline and the `human-only` pipeline.
# 
# ## Note - we use the human-only pipeline reads for all downstream analyses for these two samples
# 
# Specifically, this means that the entire analysis pipeline is performed using variants called from the human-only pipeline for sampels `004` and `005`.

# In[1]:


import os
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


human_dir = os.path.join('results', 'annotated_vcfs_humanonly')
disambig_dir = os.path.join('results', 'annotated_merged_vcfs')

human_samples = ['004-primary', '005-primary']


# In[4]:


human_sample_dict = {}

for pdx_dir, pdx_type in [[human_dir, 'human'], [disambig_dir, 'disambig']]:
    for human_sample in human_samples:
        file = os.path.join(pdx_dir,
                            '{}.annotated.hg19_multianno.csv'.format(human_sample))
        pdx_df = pd.read_csv(file)
        pdx_key = '{}_{}'.format(pdx_type, human_sample)
        human_sample_dict[pdx_key] = pdx_df


# In[5]:


human_004 = set(human_sample_dict['human_004-primary'].cosmic70)
disambig_004 = set(human_sample_dict['disambig_004-primary'].cosmic70)

human_005 = set(human_sample_dict['human_005-primary'].cosmic70)
disambig_005 = set(human_sample_dict['disambig_005-primary'].cosmic70)


# ## Visualize the difference in called variants

# In[6]:


venn2([human_004, disambig_004], set_labels = ('human', 'disambig'))
plt.title("COSMIC Variants in Sample 004")

file = os.path.join('figures', 'human_only_venn_004.png')
plt.savefig(file)


# In[7]:


venn2([human_005, disambig_005], set_labels = ('human', 'disambig'))
plt.title("COSMIC Variants in Sample 005")

file = os.path.join('figures', 'human_only_venn_005.png')
plt.savefig(file)


# ## What are the variants themselves?

# In[8]:


mouse_only_004 = disambig_004 - human_004
human_sample_dict['disambig_004-primary'].query('cosmic70 in @mouse_only_004')


# In[9]:


human_only_004 = human_004 - disambig_004
human_sample_dict['human_004-primary'].query('cosmic70 in @human_only_004')


# In[10]:


mouse_only_005 = disambig_005 - human_005
human_sample_dict['disambig_005-primary'].query('cosmic70 in @mouse_only_005')


# In[11]:


human_only_005 = human_005 - disambig_005
human_sample_dict['human_005-primary'].query('cosmic70 in @human_only_005')

