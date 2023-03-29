# scASAP: single-cell ATAC Sequencing Analysis Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[Pipeline status](https://github.com/cfusterot/scASAP/commits/master)

**scASAP** is a Snakemake pipeline that performs a comprehensive single-cell ATAC-seq data analysis, covering from the basic steps (QC, alignment, quantification) to more advanced downstream analyses.

**scASAP** makes extensive use of Snakemake's integration with the conda package manager, to automatically take care of software requirements and dependencies. Furthermore, it executes `cellranger count` on FASTQ files produced using 10X Genomics kits. The pipeline's flexibility allows users to adjust it to the distinct experimental characteristics. 

## Authors

* Coral Fustero-Torre

## Setup

For setting up the pipeline, three configuration files need to be modified. A general description of these files follows. See the *Usage* section for more details.

### Configuration files

* **config.yaml** contains all pipeline parameters.
* **units.tsv** contains information on the samples to be analysed and their paths.
* **samples.tsv**: contains information about the experimental conditions.

### Input files

* raw data in gzip compressed FASTQ files

## Usage 

### 1. Set up the environment 

**scASAP** requires the conda package manager in order to work. Please install conda by following the [bioconda installation instructions](http://bioconda.github.io/user/install.html#install-conda). In addition, of course, it is essential to install Snakemake; following the steps in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

To run the pipeline, the user needs to create the conda environments first, which will take some minutes.
This step is done automatically using this command:

    snakemake --use-conda --conda-create-envs-only --conda-frontend mamba


### 2. Download **scASAP** repository from Gitlab.
Use git clone command to create a local copy. 

    git clone https://github.com/cfusterot/scASAP.git

### 3. Configure the pipeline.

Before executing the pipeline, the users must configure it according to their samples. To do this, they must fill these files:

> TIP: different analysis can be run using just one cloned repository. This is achieved by changing the outdir and logdir in the configuration file. Also different parameters values can be used in the different analysis.

#### **a. config.yaml**

This is the pipeline configuration file, where you can tune all the available parameters to customise your scATAC-seq analysis. 

#### **b. units.tsv**

This file is used to configure the FASTQ input files.

An example file ([template_units.tsv](https://github.com/cfusterot/scASAP/master/template_units.tsv) is included in the repository.

Rename it to `units.tsv` and edit its contents according to the following table:

| **Field name** 	| **Description**                  |
|------------	|-----------------------------------------------------	|
| **sample**     	| Sample name (must match the sample name specified in *samples.tsv*).         	|
| **unit**       	| Lane identifier from which the samples will be concatenated.|
| **fq1**        	| Path to FASTQ file for read 1.  	|
| **fq2**        	| Path to FASTQ file for read 2.    | 
| **fq3**        	| Path to FASTQ file for read 3.    | 

#### **c. samples.tsv**

This table contains the name of each sample and the experimental condition it belongs to. 

An example file ([template_samples.tsv)](https://github.com/cfusterot/scASAP/master/template_samples.tsv)) is included in the repository. Rename it to `samples.tsv` and edit its contents. 

### 4. Run the pipeline.

Once the pipeline is configured and conda environments are created, the user just needs to run scASAP.

    snakemake --use-conda --use-envmodules --jobs 32 

The mandatory arguments are:
* **--use-conda**: to install and use the conda environemnts.
* **--use-envmodules**: to install and use the module environments.
* **-j**: number of threads/jobs provided to snakemake.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

