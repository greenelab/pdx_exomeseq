
# coding: utf-8

# # Visualizing Read Depth Statistics
# 
# The depth and number of mapped reads are visualized in a series of plots. The data come from `mosdepth` and `samtools flagstat` commands.

# In[1]:

import os
import csv
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:

get_ipython().magic('matplotlib inline')
plt.style.use('seaborn-notebook')


# ## Mosdepth calculations of Coverage and Depth

# In[3]:

depth_dir = os.path.join('results', 'wes_stats')
depth_files = os.listdir(depth_dir)

plt.xlabel('Depth')
plt.ylabel('Proportion of Genome')
full_depth_list = []
for depth_file in depth_files:
    # Load and process file
    depth_file = os.path.join(depth_dir, depth_file)
    depth_df = pd.read_table(depth_file, names=['chrom', 'depth', 'prop_cov'])
    depth_df = depth_df.assign(sample_id = depth_file)
    
    # Set a reasonable cutoff to view bulk of distribution
    depth_sub_df = depth_df[depth_df['depth'] < 300]
    
    # There are chromosome specific estimates, here, view total
    ax = plt.plot('depth', 'prop_cov', data=depth_sub_df[depth_sub_df['chrom'] == 'total'])
    
    # Save processed files to list for later extraction
    full_depth_list.append(depth_df)

depth_fig_file = os.path.join('figures', 'mosdepth_estimation.pdf')
plt.savefig(depth_fig_file)


# In[4]:

# Get full depth matrix
full_depth_df = pd.concat(full_depth_list)

# 50% of the exome is covered at what depth?
half_depth_df = full_depth_df[full_depth_df['prop_cov'] == 0.5]
half_depth_df = half_depth_df[half_depth_df['chrom'] == 'total']


# In[5]:

half_depth_df.describe()


# In[6]:

# 75% of the exome is covered at what depth?
threeq_depth_df = full_depth_df[full_depth_df['prop_cov'] == 0.75]
threeq_depth_df = threeq_depth_df[threeq_depth_df['chrom'] == 'total']


# In[7]:

threeq_depth_df.describe()


# ## General flagstat read mapping distributions

# In[8]:

read_stat_dir = os.path.join('results', 'read_counts')
read_stat_files = os.listdir(read_stat_dir)


# In[9]:

read_stat_file = os.path.join(read_stat_dir, read_stat_files[0])


# In[10]:

total_read_counts = []
mapped_read_counts = []

for read_stat_file in read_stat_files:
    read_stat_file = os.path.join(read_stat_dir, read_stat_file)

    with open(read_stat_file) as csvfile:
        readstat_reader = csv.reader(csvfile, delimiter=' ')
        line_idx = 0
        for row in readstat_reader:
            if line_idx == 0:
                total_read_counts.append(row[0])
            if line_idx == 2:
                mapped_read_counts.append(row[0])
            line_idx +=1


# In[11]:

# Get flagstat dataframe
read_depth_df = pd.DataFrame([total_read_counts, mapped_read_counts],
                             index=['total', 'mapped'], dtype='float64').T


# In[12]:

read_depth_df.plot(x='total', y='mapped', kind='scatter')


# In[13]:

# Output distribution of total and mapped reads
read_depth_df.describe()

