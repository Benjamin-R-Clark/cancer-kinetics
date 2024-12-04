#!/bin/bash

nextflow run nf-core/scrnaseq -resume  -profile docker -c nextflow.config  --outdir gs://bc-scseq-data/nextflow-scrnaseq-final \
--input samplesheet.csv \
--star_index gs://bc-scseq-data/genomes/index \
--genome GRCh38 \
--gtf gs://bc-scseq-data/genomes/hg38.knownGene.gtf \
--protocol 10XV2 \
--aligner star \

