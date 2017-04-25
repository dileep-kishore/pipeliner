
==================
PIPELINER: A Versatile RNA-Seq Pipeline
==================

Features
========

* Modular directory structure: It is designed to generate automated result directory based on the names of the samples and tools used to process them
* Platform independent: It is bundled with an anaconda repository which contains pre-compiled tools as well as pre-built environments that can use used directly.
* Modular architecture: It allows the expert users to customize, modify processes, or add additional tools based on their needs.
* Automated job parallelization, job recovery, and reproducibility

Running the pipeline
====================


**Pre-requisites**

Pipeliner requires Nextflow and Java 7 or greater for implementation. All other tools for implementation are wrapped in an environment described below. 

*1. Download Nextflow*

Make sure Java 7  or 7 is installed and install Nextflow in directory you plan to work in::

  java -version
  curl -fsSL get.nextflow.io | bash

*2. Download Conda*

`Conda` is available through 'Anaconda <https://www.continuum.io/downloads>`_ and  `Miniconda <http://conda.pydata.org/miniconda.html>`_.

If you're using a module system such as on the shared computing cluster (SCC) at Boston University you can just load a preinstalled version::

   module purge
   module load anaconda2/4.3.0


*2. Create Conda Environment (rna_env)*

To install a basic development environment, download rna_env from Pipeliner repository in anaconda cloud::

  conda env create pipeliner/rna_env
  source activate rna_env
  
There is also an option to create a './rna_env folder and store all the files needed to run the pipeline as such::

  conda create \
    -p ./rna_env \
    -c https://anaconda.org/Pipeliner \
    --yes \


Running the pipeline
====================

Activate the environment (follow instructions above to create environment)::
 
  source activate rna_env

OR::

  source activate ./rna_env
  
Pipeliner consists of a main nextflow script parametrized using a configuration file. The configuration file includes all parameters necessary to run the pipeline including  parameters to direct the path of files and results, as well as selecting specific tools and processes to run. The example gallus gallus (chicken) dataset is all within the ggal folder. All files necessary to run this example is in the folder "Gallus_Example".


Once the appropriate paths and tools have been set in the configuration file, you can run the pipeline with the configuration file::

  nextflow main.nf -c config




