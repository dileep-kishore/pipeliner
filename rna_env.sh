# @Author: Dileep Kishore <dileep>
# @Date:   December 7, 2016 9:24:25 AM
# @Filename: rna_env.sh
# @Last modified by:   dileep
# @Last modified time: December 7, 2016 10:30:23 AM

#!/bin/bash

#If using a module system such as SCC: module load a preinstalled version of anaconda
# module load anaconda2/4.2.0

#General way to create a conda environment
#Create virtual environment
conda create -n rna_env python=3 anaconda --yes

#Activate virtual environment
source activate rna_env

# Add the bioconda channel
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

#Install additional python packages
conda install --name rna_env [package]

#Install non-conda packages
pip install [package]

#Deactivate virtual environment
deactivate rna_env

#Following Montilab's example: Create an environment folder that stores all the files needed to run the pipeline
#channel represents bioconda channel we can create with specific packages needed

conda create \
    -p ./rna_env \
    -c [channel] \
    --yes \
    pipeliner
