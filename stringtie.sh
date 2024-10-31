#!/bin/bash
pacbio_gtf=GBC_filter_rules.filtered.gtf
stringtie ${samplename}_Aligned.out.sort.bam -p 60 -G ${pacbio_gtf} -e -A ${samplename}_gene_abund.tab -B -o ${samplename}.gtf -l ${samplename} -C ${samplename}_cov_refs.gtf
