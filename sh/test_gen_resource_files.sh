#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN ../R/generateResourceFiles.R --prefix="test_" --cSNP_file="../DATA/support/cSNPs_w_ENSG.e75_chr22.bed" --interaction_file="../DATA/chic/mifsud_et_al.pm.chr22.tab" --pchic.thresh=5 --res_frag_file='../DATA/support/Digest_Human_HindIII_chr22.bed' --region_bed_file='../DATA/support/0.1cM.regions.b37_chr22.bed' out_dir='../DATA/RDATA/'
