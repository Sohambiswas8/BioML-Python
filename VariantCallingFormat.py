#!/usr/bin/env python
# coding: utf-8

# # Operations with VCF files

# Python script for basic operations on VCF files. 
# 
# Downloaded file: https://browser.genomeasia100k.org/service/web/download_files/22.substitutions.annot.cont_withmaf.vcf.gz from GenomeAsia 1000k. 

# # Read VCF file

# Import libraries

# In[5]:


# import libraries
# either, 'pip install scikit-allel' or run here in Jupyter

get_ipython().system('pip install scikit-allel')
import allel
import h5py
import zarr
import numcodecs


# In[6]:


# Read VCF file 
# simple

with open('22.substitutions.annot.cont_withmaf.vcf', mode='r') as vcf:
    print(vcf.read())


# In[7]:


# Read VCF file
# with 'read_vcf()' of allel package
# once check version

callset = allel.read_vcf('22.substitutions.annot.cont_withmaf.vcf')

# callset object will give array of numpy keys
print(callset)


# In[8]:


# We see there are several VCF columns as keys:
# CHROM, POS, ID, REF, ALT, QUAL. Other keys can be there like FILTER, INFO etc.
# The keys will be return alphabetically not by column in actual VCF file
# We can call each key and see its values i.e. values in each column

sorted(callset.keys())


# In[12]:


# We see only 7 columns are there
# call each key

print(callset['variants/CHROM'])
len(callset['variants/CHROM'])


# In[14]:


# We see 860546 SNPs are there
# call IDs
callset['variants/ID']


# In[15]:


# Similarly we can call other keys
# Next how to convert the VCF dataset into dataframe


# In[16]:


# import pandas
import pandas as pd
df = allel.vcf_to_dataframe('22.substitutions.annot.cont_withmaf.vcf')
df


# In[17]:


# check first 100 rows

df.head(100)


# In[19]:


# We see ALT column has been broken into ALT_1, ALT_2 and ALT_3.
# This means if the VCF file has multiple alternative alles, then those allele values will be broken into seperate columns
# Check if we have multiple alternative alleles
# Only print the values for ALT_2 and ALT_3

df['ALT_1']


# In[20]:


print(df.ALT_2)


# In[21]:


print(df.ALT_3)


# In[22]:


# Still ALT_2 and ALT_3 can have values in alleles that are not visible here
# Find the better way to look into


# But wait !!!!!!! I just realised. 
# The data that I downloaded only include substituions not InDels (insertions/deletions). Therefore the data will only include Ref and single ALT allele.

# The data looks like quality checked. Still I need to be sure. 
# As of now save it in CSV format.

# In[23]:


# Save in CSV

allel.vcf_to_csv('22.substitutions.annot.cont_withmaf.vcf', 'genomeasia.chr22.subst.cont_withmaf.csv',
                fields = ['CHROM', 'POS', 'ID','REF','ALT','QUAL','FILTER_PASS'])


# In[24]:


with open('genomeasia.chr22.subst.cont_withmaf.csv', mode='r') as f:
    print(f.read())


# # Use with 1000 Genomes Project data

# Download data from from 1000 Genomes project. Data in VCF format with chromosome wise SNP information 

# Idea: 
# 1. First save the file path from FTP into an object (vcf_path)
# 2. Then look into the file size and the number of lines in the file

# In[25]:


vcf_path = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
get_ipython().system('ls -lh {vcf_path}')


# In[ ]:




