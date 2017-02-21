#!/bin/bash

# @Author: Dileep Kishore <dileep>
# @Date:   December 7, 2016 9:24:25 AM
# @Filename: rna_env.sh
# @Last modified by:   dileep
# @Last modified time: February 7, 2017 8:22:42 PM


#If using a module system such as SCC: module load a preinstalled version of anaconda
# module load anaconda2/4.2.0

## General way to create a conda environment
# Prepending the path just in case
export PATH=$HOME/anaconda3/bin:$PATH
#Create virtual environment
#FIXME: Python2 vs. Python3
conda create -n rna_env python=2 anaconda --yes

#Activate virtual environment
source activate rna_env

# Add the bioconda channel
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# Install additional python packages
conda install -n rna_env fastqc cutadapt trim-galore star rseqc stringtie multiqc --yes
# ## install nextflow as well (might not be the lastest version)
# conda install -n rna_env nextflow --yes

## Install non-conda packages
# pip install [package]

# Export environment to yaml files
conda env export > rna_env.yml

# # Deactivate virtual environment
# deactivate rna_env

# #Following Montilab's example: Create an environment folder that stores all the files needed to run the pipeline
# #channel represents bioconda channel we can create with specific packages needed
# conda create -p ./rna_env -c [channel] --yes pipeliner
