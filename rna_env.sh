#!/bin/bash

#If using a module system such as SCC: module load a preinstalled version of anaconda
module load anaconda2/4.2.0

#General way to create a conda environment
#Create virtual environment
conda create -n rna_env python=2.7 anaconda --yes

#Activate virtual environment
source active rna_env

#Install additional python packages
conda install --name rna_env [package]

#Install non-conda packages
pip install [package]

#Install packages from Anaconda.org i.e. Install package from a channel
conda install -c [channel] [package]


#Deactivate virtual environment
deactivate rna_env

#Following Montilab's example: Create an environment folder that stores all the files needed to run the pipeline
#channel represents bioconda channel we can create with specific packages needed

conda create \
    -p ./rna_env \
    -c [channel] \
    --yes \
    pipeliner
