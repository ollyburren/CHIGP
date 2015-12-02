#CHiGP

CHiGP stands for **C**apture **Hi**-C **G**ene **P**rioritisation. It's purpose is to take a bunch of p-vals from a genome-wide association study and integrate these with capture Hi-C data to prioritise genes. The computational engine is written in R and whilst it should be possible to process on a laptop, things will go considerably faster if you have access to a compute cluster.

#Prerequisites

  * Linux or Mac OS X (preferably the former)
  * R (developed with 3.2.2)
    * data.table
    * snpStats (bioconductor)
    * GenomicRanges (bioconductor)
  * make
  * (optional) PERL
    * Macd queue software and dependencies.
  * [htslib](https://github.com/samtools/htslib) for tabix


#Data files

  * capture Hi-C data (we support a few tab delimited formats)
    * PeakMatrix format
    * WashU format
  * 1000 Genomes VCF files.
  * coding SNP bed file
  
#Initial setup

This assumes that you have R installed and have the ability to install packages. Firstly head to [bioconductor](https://www.bioconductor.org/install/) and read how to install snpStats and GenomicRanges dependencies. After you have these install data.table

```
install.packages("data.tables")
```
Next grab the test file bundle. This contains a cut down set of real data taken from Misfud et al. and GWAS from Okada et al.
```
cd CHIGP
curl -s http://www.immunobase.org/Downloads/CHIGP/test_data.tgz
tar -xvzf test_data.tgz
```
This should have created a dir 'data' within the current CHIGP dir 
open up ~/.Rprofile in a text editor and add the following lines
```
CHIGP_DIR<-<PATH TO YOUR CHIGCP DIR>
CHIGP_OUTDIR<-<PATH TO WHERE YOU WANT CHIGP TO WRITE DATA>
```



#Test run
Once you have all of the above dependencies and data installed you should be able to take CHIGP for a spin. Running CHICGP consists of two steps
  1. Convert p-vals to posterior probabilities using a reference set of genotypes.
  2. Integrate CHiC interaction and other functional data to prioritise genes.
For a given GWAS it should only be neccessary to run step 1 only once (unless of course you want to fiddle with parameters). You might want to run step 2 multiple times depending on what capture HiC datasets are available to you.

## Converting p-vals to posterior probabilities.

For a detailed discussion on how we do this please see Wakefield(2009), Giombartelli(2014) and Pickrell(2014). In fact the code used here is adapted from Giombartelli. Briefly, we split the genome into blocks based on a recombination frequency of 0.1cM, for each region we then estimate the minor allele frequency (MAF) of each SNP within that region available in the reference genotyping dataset. **Note that if a SNP in your GWAS is not in your reference genotyping set then it will be ommited from future steps**. We use this to estimate an approximate Bayes factor for that SNP to be causal given the data, we can then for a given region estimate the posterior probability (PPi) that a SNP is causal. 

How to do this on the command line.
```
Rscript --vanilla PMI.R \\
 --region_file=./0.1cM_shuffled_regions/0.1cMcs \\
 --out_dir=/home/uname/captureHIC_support/GWAS/PMIPV/ --gwas_tbx=./GWAS/QC_TABIX/Platelet_Volume.bed.gz \\
 --gwas_type=QUANT --n_samples=18600 --prop_cases=1 --kg_compress_dir=./1KGenome/VCF/EUR/by.chr.phase3/ALL. \\
 --kg_compress_suffix=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz \\
 --tabix_bin=/bin/tabix --pi_1=0.0001 
