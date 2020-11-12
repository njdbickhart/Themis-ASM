# Themis-ASM
---

A snakemake pipeline to compare assemblies with comparative alignment data. The output of this pipeline consists of many popular assembly metrics in an easy to digest report which highlights large and small-scale differences in assembly content.

## Table of Contents:
* [Quick Start guide](#quick)
* [What is Themis-ASM?](#what)
* [What do I need to run Themis-ASM?](#who)
	* [Input file requirements](#input)
	* [Software installation requirements](#software)
	* [Hardware requirements and recommendations](#hardware)
* [Running Themis-ASM for comparative assembly analysis](#running)
	* [Single server instructions](#single)
	* [Distributed server instructions](#cluster)
	* [Cluster format table](#format)
* [What is the output of Themis-ASM?](#output)

<a name="quick"></a>
## Quick Start guide

Ensure that you have python 3.6+ and snakemake installed on your path along with any other listed prerequisites. Run the pipeline by invoking the wrapper script with the following example arguments:

```bash
python ThemisASM.py  \
	# Enter the full path to all assembly files
	-a /path/to/assembly1.fa  \
	-a /path/to/assembly2.fa  \

	# In the order in which they were entered, assign a name to the assembly fastas
	-n Assembly1	 \
	-n Assembly2     \

	# Select a busco database for comparison
	-b mammalia_odb10 \

	# Enter a comma delimited list of short paired-end reads to align to the assemblies
	-f /path/to/short_reads_R1.fq,/path/to/short_reads_R2.fq \

	# For each "s" entry, give a name to that sample
	-s reference_animal  \

	# If you have a limited number of CPU threads, specify a number of concurrent jobs appropriately
	-j 20

	### Only add the following argument if you are running on a cluster! 
	# If you are submitting on a distributed system, just enter a cluster submission prefix string similar to the following example:
	-c 'sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}'
```

<a name="what"></a>
## What is Themis-ASM?

Named after the Greek Titaness of natural law, this workflow attempts to embody the aspect of [Themis](https://en.wikipedia.org/wiki/Themis) who was known as the "lady of good counsel." In practice, determining if your assembly is suitable for use by outside research groups often requires comparative analysis to identify potential flaws (or misassemblies). Themis-ASM is a workflow that combines several useful quality assessment tools for comparative analysis of assemblies against a reference version or even prior versions of the same assembly. 

<a name="who"></a>
## What do I need to run Themis-ASM?

Themis-ASM has several requirements for use. Let's first discuss the input datasets that are required prior to setting up the pipeline. We assume that you have at least two assemblies for your comparison and that you also have sequenced your reference/target sample with short-read sequence data. You will also need to select a [BUSCO](https://busco.ezlab.org/) database to align against your assembly for gene content statistics.

<a name="input"></a>
#### Input file requirements:

* At least two assembly fasta files 
* At least one paired-end short-read sequence dataset (fastq format) to align against your assemblies
* A BUSCO database label (ie. mammalia_odb10)

If you do not have short-read data for your assemblies, it is possible to use publicly available data from online repositories such as the [SRA](https://www.ncbi.nlm.nih.gov/sra). The files must be paired-end at this time (there is no support for single-end reads... yet!).

Now that we've detailed the dataset inputs for Themis-ASM, let's talk about software requirements. As we will touch on this in the next section, here is a brief listing of the requirements to run the pipeline:

<a name="software"></a>
#### Software installation requirements:

* Python 3.6+ installation
	* Please make sure that Snakemake is installed and on your path!
	* Also, please make sure that you have installed the following python packages:
		* itertools
		* pysam
* [Merqury](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02134-9) must be on your path
* [BWA v 0.7+](https://github.com/lh3/bwa) on your path
* [Samtools v 1.9+](http://www.htslib.org/) on your path
* [Miniconda or an equivalent "conda-based" package manager](https://docs.conda.io/en/latest/miniconda.html) on your path

#### NOTE: Themis-ASM has only been tested on Linux-based systems! If you are trying to run the package on Windows, Mac or your smart toaster, we cannot provide support. 

Themis-ASM can run on single computers or on distributed systems. Running tasks on clusters is accomplished through built-in job submission code within the Snakemake workflow system. Here are our required and recommended computational statistics for using the pipeline:

<a name="hardware"></a>
#### Hardware requirements and recommendations:

| Component  | Required | Recommended |
| :--------- | :------- | :---------- |
| RAM | 48 gigabytes | 300+ gigabytes recommended|
| CPU threads | 8 threads |70+ recommended |
| Storage per assembly| ~ 100 gigabytes | 200+ gigabytes |
| Storage per short-read dataset | ~ 100 gigabytes | 200+ gigabytes |

Hard disk footprints are estimated for genomes that are 1 gigabases in size and for short-read datasets that are 100 million reads of 150 x 150 paired end Illumina data, so please adjust these estimates accordingly based on your datasets! 

<a name="running"></a>
## Running Themis-ASM for comparative assembly analysis

We've created a script that should automate the tedious process of setting up your configuration files for the workflow. The name of the configuration file will always be **default.json** and it will reside in the directory in which you run the **ThemisASM.py** workflow script. The script is designed to set up your configuration and queue up the main workflow snakefile automatically, with an option to use cluster-based job submission instructions. Here are brief instructions on how to run Themis-ASM on single nodes and distributed systems:

<a name="single"></a>
#### Single server instructions

First, you must start in your designated working directory. To run the pipeline on a single server, run a modification of the following toy command:

```bash
python ThemisASM.py  \
	# Enter the full path to all assembly files
	-a /path/to/assembly1.fa  \
	-a /path/to/assembly2.fa  \

	# In the order in which they were entered, assign a name to the assembly fastas
	-n Assembly1	 \
	-n Assembly2     \

	# Select a busco database for comparison
	-b mammalia_odb10 \

	# Enter a comma delimited list of short paired-end reads to align to the assemblies
	-f /path/to/short_reads_R1.fq,/path/to/short_reads_R2.fq \

	# For each "s" entry, give a name to that sample
	-s reference_animal  \

	# If you have a limited number of CPU threads, specify a number of concurrent jobs appropriately
	-j 20

```

The script will then create the necessary configuration files and attempt to run the snakemake workflow itself. If the workflow is interrupted, or terminates prematurely, snakemake "unlock" functionality is wrapped within the script:

```bash
python ThemisASM.py \
	# Enter all of your previous arguments

	# Use this to unlock your directory if snakemake exited prematurely: 
	--unlock

	# Then, after the directory is unlocked, use the following flag to resume your job
	--resume
```

<a name="cluster"></a>
#### Distributed server instructions

In order to queue tasks on clusters, snakemake must be fed a job submission prefix that is specific to your cluster system. Themis-ASM already has suggested settings for specific rules in the following file in the base level of this repository:

> cluster.json

If you queue the wrapper script using the "--cluster_string" flag, these settings will be automatically applied to your jobs. More details are in the table below. An example of submitting distributed jobs on the Slurm job management system is provided below:

```bash
python3 ThemisASM.py  \
	...  # All previous arguments found in the single server example

	# Enter a formatted prefix string for running jobs on your cluster:
	-c 'sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q {cluster.qos} -o {cluster.stdout}'
```

If you want to ignore Themis-ASM's default settings (or they don't apply to your cluster setup) you can set these values to flat instances or ignore them like so:

```bash
# Same as above

	-c 'sbatch --nodes=1 --ntasks-per-node={cluster.ntasks-per-node} --mem=50000 --partition=my_partition'
```

<a name="format"></a>
#### Cluster format table

Here is a brief listing of all of the cluster settings and format strings used by Themis. Please adjust them according to your needs and cluster environment:

| Resource  | Format String  | Default Value |
| :---  | :--- | :--- |
| Node count | {cluster.nodes} | 1 |
| Cpu threads | {cluster.ntasks-per-node} | 1 |
| Memory (Mb) | {cluster.mem} | 5000 |
| Partition (or account)   | {cluster.partition} | "priority", but you should probably change! |
| Job STDOUT | {cluster.stdout} | "logs/{rule}.{wildcards.sample}.stdout" |
| Job Name  | {cluster.jobname} | "{rule} [{wildcards.sample}]" |

#### NOTE: Themis-ASM has only been tested on Slurm HPCs. We welcome input on how to adjust the pipeline for SGE or other systems, with the notable exception of smart refrigerators

<a name="output"></a>
## What is the output of Themis-ASM?

Assuming everything went well (congratulations!), you should notice a whole bunch of files in your working directory! Here is a very brief listing of the folders you should expect:


```bash
workingdir
|__assembly_qc.zip   <- This is your packaged, final report for your assembly comparisons!
|__busco
|  |___asm_name
|      |___busco_summary.txt  <- These are the busco scores for the assembly
|__calls			 <- this directory contains intermediate files for your mapped sequence read calls
|__fastas			<- This contains assembly alignment indicies
|__final			 <- This is the folder that contains all of the plots and tables. It's packaged in the zip file above
|__logs			  <- Check the log files for the progress on your pipeline!
|__mapped			<- This contains intermediate alignment files for the assemblies
|__merqury			
|  |___asm_name	  <- This contains intermediate merqury files for each assembly
|__summary_page.html		<- this is a html report file that is packaged in the zip file listed above

```

The most important files are within the "final" folder and the zip archive produced by the pipeline. You can share the entire archive with colleagues, but please ask them to unzip the folder before viewing the html summary page!

