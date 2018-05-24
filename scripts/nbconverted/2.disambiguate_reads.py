
# coding: utf-8

# # Disambiguate Reads Summary
# 
# Visualizes the proportion of reads filtered to human, mouse, and ambiguous by the `disambiguate` tool. The notebook will output a figure describing these proportions across all replicates from all 30 samples.

# In[1]:


import os
import pandas as pd

import plotnine as gg


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[3]:


# Load Phenotype Data 
file = 'pdx_phenotype.csv'
pheno_df = pd.read_table(file, sep=',')

# Create variable for ID updating
id_updater = dict(zip([x[0] for x in pheno_df.read_id.str.split('_')],
                      pheno_df.final_id))
id_updater


# In[4]:


summary_dir = os.path.join('results', 'disambiguate_summary')
summary_files = [os.path.join(summary_dir, x) for x in os.listdir(summary_dir)]


# In[5]:


summary_list = []
for summary_file in summary_files:
    summary_list.append(pd.read_table(summary_file))


# In[6]:


summary_df = pd.concat(summary_list)
summary_df.columns = ['sample', 'human', 'mouse', 'ambiguous']
summary_df = summary_df.assign(base_sample = [id_updater[x[0]] for x in summary_df['sample'].str.split('_')])
summary_df['sample'] = summary_df['sample'].replace(id_updater, regex=True)
summary_df = summary_df.assign(lane = [x[2] for x in summary_df['sample'].str.split('_')])

summary_df.head()


# In[7]:


total_reads = summary_df['human'] + summary_df['mouse'] + summary_df['ambiguous']
human_percent = (summary_df['human'] / total_reads) * 100
mouse_percent = (summary_df['mouse'] / total_reads) * 100
ambig_percent = (summary_df['ambiguous'] / total_reads) * 100

summary_df = summary_df.assign(human_percent = human_percent.round(1))
summary_df = summary_df.assign(mouse_percent = mouse_percent.round(1))
summary_df = summary_df.assign(ambig_percent = ambig_percent.round(1))


# In[8]:


summary_df.to_csv(os.path.join('results', 'full_disambiguate_summary.tsv'), sep='\t')
summary_df.head()


# In[9]:


summary_melt_df = summary_df.melt(id_vars=['base_sample', 'lane', 'sample', 'human_percent',
                                           'mouse_percent', 'ambig_percent'],
                                  value_vars=['human', 'mouse', 'ambiguous'],
                                  var_name='species', value_name='pairs')
summary_melt_df.head()


# In[10]:


summary_melt_df.loc[summary_melt_df['species'] != 'human', 'human_percent'] = ''
summary_melt_df.loc[summary_melt_df['species'] != 'mouse', 'mouse_percent'] = ''
summary_melt_df.loc[summary_melt_df['species'] != 'ambiguous', 'ambig_percent'] = ''


# In[11]:


p = (
    gg.ggplot(summary_melt_df, gg.aes(x='lane', y='pairs', fill='species')) +
    gg.geom_bar(stat='identity', position='stack') +
    gg.geom_text(gg.aes(y=3.5e7, label='ambig_percent'), size=4.5) +
    gg.geom_text(gg.aes(y=2.7e7, label='human_percent'), size=4.5) +
    gg.geom_text(gg.aes(y=1.9e7, label='mouse_percent'), size=4.5) +
    gg.facet_wrap('~ base_sample') +
    gg.theme_bw() +
    gg.theme(axis_text_x=gg.element_text(angle='90'),
             axis_text=gg.element_text(size=8),
             axis_title=gg.element_text(size=14))
    )
p


# In[12]:


figure_file = os.path.join('figures', 'disambiguate_results.pdf')
gg.ggsave(p, figure_file, height=5.5, width=6.5, dpi=500)

