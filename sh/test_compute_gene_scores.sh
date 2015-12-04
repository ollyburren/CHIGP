#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN ../R/computeGeneScore.R --pmi_file=../DATA/out/0.1cM_chr22.imp --out_file=../DATA/0.1cM_chr22.geneScores.tab --csnps=../DATA/RDATA/test_cnps.by.ld.RData --int=../DATA/RDATA/test_interactions.RData --frags=../DATA/RDATA/test_frags.by.ld.RData --target.gene.cSNPs.only=1
