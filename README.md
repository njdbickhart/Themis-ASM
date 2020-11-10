# Themis-ASM
---

A snakemake pipeline to compare assemblies with comparative alignment data. 

## Table of Contents:

## What is Themis-ASM?

Named after the Greek Titaness of natural law, this workflow attempts to embody the aspect of [Themis](https://en.wikipedia.org/wiki/Themis) who was known as the "lady of good counsel." In practice, determining if your assembly is suitable for use by outside research groups often requires comparative analysis to identify potential flaws (or misassemblies). Themis-ASM is a workflow that combines several useful quality assessment tools for comparative analysis of assemblies against a reference version or even prior versions of the same assembly. 

## What do I need to run Themis-ASM?

Themis-ASM has several requirements for use. Let's first discuss the input datasets that are required prior to setting up the pipeline. We assume that you have at least two assemblies for your comparison and that you also have sequenced your reference/target sample with short-read sequence data. You will also need to select a [BUSCO](https://busco.ezlab.org/) database to align against your assembly for gene content statistics.

#### Input file requirements:

* At least two assembly fasta files 
* At least one paired-end short-read sequence dataset (fastq format) to align against your assemblies
* A BUSCO database label (ie. mammalia_odb10)

If you do not have short-read data for your assemblies, it is possible to use publicly available data from online repositories such as the [SRA](https://www.ncbi.nlm.nih.gov/sra). The files must be paired-end at this time (there is no support for single-end reads... yet!).

Now that we've detailed the dataset inputs for Themis-ASM, let's talk about software requirements. As we will touch on this in the next section, here is a brief listing of the requirements to run the pipeline:

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

#### Hardware requirements and recommendations:

| Component  | Required | Recommended |
| :--------- | :------- | :---------- |
| RAM | 48 gigabytes | 300+ gigabytes recommended|
| CPU threads | 8 threads |70+ recommended |
| Storage per assembly| ~ 100 gigabytes | 200+ gigabytes |
| Storage per short-read dataset | ~ 100 gigabytes | 200+ gigabytes |

Hard disk footprints are estimated for genomes that are 1 gigabases in size and for short-read datasets that are 100 million reads of 150 x 150 paired end Illumina data, so please adjust these estimates accordingly based on your datasets! 

## Running Themis-ASM for comparative assembly analysis