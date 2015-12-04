#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
TABIX_BIN=`which tabix`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN ../R/computePPi.R --kg_compress_dir=../DATA/1kgenome/VCF/EUR/by.chr.phase3/ALL. --region_file=../DATA/support/0.1cM_chr22 --gwas_tbx=../DATA/gwas/OKADA_RA_22.bed.gz --gwas_type=QUANT --n_samples=58293 --prop_cases=0.25 --kg_compress_suffix=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz --tabix_bin=$TABIX_BIN --pi_i=0.0001 --out_dir=../DATA/out/
