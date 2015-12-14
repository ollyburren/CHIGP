#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN ../R/computeGeneScore.R --pmi_file=../DATA/out/0.1cM_chr22.pmi --out_file=../DATA/out/0.1cM_chr22.geneScores.pmi.tab --csnps=../DATA/RDATA/test_csnps.by.ld.RData --int=../DATA/RDATA/test_interactions.RData --frags=../DATA/RDATA/test_frags.by.ld.RData --target.gene.cSNPs.only=1
