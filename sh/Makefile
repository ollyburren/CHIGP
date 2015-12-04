## makefile for running test analysis

GRPATH = ../../
DATADIR = ../DATA/
RDATADIR = $(DATADIR)RDATA/

#all:../DATA/out/0.1cM_chr22.geneScores.tab

$(DATADIR)/out/0.1cM_chr22.geneScores.tab:$(RDATADIR)/test_interactions.RData $(RDATADIR)/test_frags.by.ld.RData $(RDATADIR)/test_cnps.by.ld.RData
	./test_compute_gene_scores.sh
	
$(RDATADIR)/test_interactions.RData:$(DATADIR)/out/0.1cM_chr22.imp
	./test_gen_resource_files.sh
	
$(RDATADIR)/test_frags.by.ld.RData:$(DATADIR)/out/0.1cM_chr22.imp
	./test_gen_resource_files.sh
	
$(RDATADIR)/test_cnps.by.ld.RData:$(DATADIR)/out/0.1cM_chr22.imp
	./test_gen_resource_files.sh
	
$(DATADIR)/out/0.1cM_chr22.imp: 
	./test_ppi.sh
	
clean:
	rm $(DATADIR)/out/0.1cM_chr22.imp $(DATADIR)/out/0.1cM_chr22.geneScores.tab $(RDATADIR)/test*.RData
