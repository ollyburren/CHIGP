#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN ../R/blockshifter.R --contacts_file=../DATA/chic/misfud_et_al.pm.chr22.tab --pmi_file=../DATA/out/0.1cM_chr22.pmi --perm_no=100 --test_tissue=GM12878 --control_tissue=CD34 --output_file=../DATA/out/mifsud_GM12878vsCD34_0.1cM_chr22.txt --metric=ppi
