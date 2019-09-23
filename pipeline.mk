SHELL=/bin/bash

PIPELINE_DIR := $(shell realpath --relative-to=. $(dir $(lastword $(MAKEFILE_LIST))))

default:: all

all::

# distclean::
# 	rm -rf data

# ------------------------------------------------------------------------

all:: data/genome_genbank/genome.fna

data/assembly_genbank :
	${PIPELINE_DIR}/download-assembly \
		${GENBANK_ACCESSION} data/assembly_genbank
data/genome_genbank/genome.fna : data/assembly_genbank
	${PIPELINE_DIR}/demux-genome ${MOLECULE_NAMES} \
		data/assembly_genbank data/genome_genbank

# ------------------------------------------------------------------------

all:: data/genome_refseq/genome.fna

data/assembly_refseq :
	${PIPELINE_DIR}/download-assembly \
		${REFSEQ_ACCESSION} data/assembly_refseq
data/genome_refseq/genome.fna : data/assembly_refseq
	${PIPELINE_DIR}/demux-genome ${MOLECULE_NAMES} \
		data/assembly_refseq data/genome_refseq

# ------------------------------------------------------------------------

all:: data/genome_prokka/annotation.gff
data/genome_prokka/annotation.gff : \
			data/genome_refseq/genome.fna \
			${PIPELINE_DIR}/run-prokka
	@rm -rf data/genome_prokka
	time ${PIPELINE_DIR}/run-prokka \
		-G ${PROKKA_GRAM} -g ${PROKKA_GENUS} -s ${PROKKA_SPECIES} \
		-S ${PROKKA_STRAIN} -p prokka \
		-o data/genome_prokka -t data/genome_prokka/temp \
		data/genome_refseq/genome.fna
	@echo ''

# ------------------------------------------------------------------------

RUN_PGAP_ARGS=
ifdef PGAP_DIR
RUN_PGAP_ARGS+= -p ${PGAP_DIR} # use existing PGAP instance
endif
#RUN_PGAP_ARGS+= -k # keep temporary files
#RUN_PGAP_ARGS+= -n # fake the run

PGAP_ARGS=
PGAP_ARGS+= --quiet
PGAP_ARGS+= --report-usage-false
# PGAP_ARGS+= --report-usage-true
# PGAP_ARGS+= --verbose

all:: data/genome_pgap/annotation.gff
data/genome_pgap/annotation.gff : \
			data/genome_refseq/genome.fna \
			${PIPELINE_DIR}/run-pgap
	time ${PIPELINE_DIR}/run-pgap ${RUN_PGAP_ARGS} \
			   -t ${PGAP_TAXID} -S ${PGAP_STRAIN} \
			   -o data/pgap data/genome_refseq/genome.fna \
			   -- ${PGAP_ARGS}
	@mkdir data/genome_pgap
	cp data/pgap/output/annot.fna data/genome_pgap/genome.fna
	cp data/pgap/output/annot.faa data/genome_pgap/proteome.faa
	cp data/pgap/output/annot.gff data/genome_pgap/annotation.gff
	cat data/pgap/output/annot.gbk \
	   | ${PIPELINE_DIR}/split-gbk -d data/genome_pgap -n

# ------------------------------------------------------------------------

all:: data/reconciliation.gff

data/reconciliation.gff : \
		${SEARCH_RESULTS} \
		data/genome_genbank/annotation.gff \
		data/genome_refseq/annotation.gff \
		data/genome_prokka/annotation.gff \
		data/genome_pgap/annotation.gff \
		${PIPELINE_DIR}/reconcile-annotations.pl
	cat ${SEARCH_RESULTS} \
	    | egrep -v '^(crap|contam)'$$'\t' \
	    | fgrep -v $$'\t'complex_peptide$$'\t' \
		    > data/pg-search-results-genome.gff
	time ${PIPELINE_DIR}/reconcile-annotations.pl \
	    data/pg-search-results-genome.gff \
	    genbank:data/genome_genbank/annotation.gff \
	    refseq:data/genome_refseq/annotation.gff \
	    prokka:data/genome_prokka/annotation.gff \
	    pgap:data/genome_pgap/annotation.gff \
	    > data/reconciliation.gff

