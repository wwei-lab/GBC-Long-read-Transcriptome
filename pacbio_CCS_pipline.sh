#!/bin/bash
###Using cell and normal tissue pacbio ccs data as examples

# 01 ccs 
# Consensus calling 
cd 2_ccs

ccs GBC_SD_merge.subreads.bam ./GBC_SD_merge.ccs.bam --min-rq 0.9
ccs L2F7_merge.subread.bam ./L2F7_merge.ccs.bam --min-rq 0.9
ccs NOZ_merge.subreads.bam ./NOZ_merge.ccs.bam --min-rq 0.9
ccs ZJU0430_merge.subreads.bam ./ZJU0430_merge.ccs.bam --min-rq 0.9
ccs OCUG.subreads.bam ./OCUG.ccs.bam --min-rq 0.9

# 02 lima 
# Primer removal and demultiplexing
cd 3_lima
lima ./2_ccs/GBC_SD_merge.ccs.bam primer.fasta GBC_SD_merge.demux.bam --isoseq --peek-guess
lima ./2_ccs/L2F7_merge.ccs.bam primer.fasta L2F7_merge.demux.bam --isoseq --peek-guess
lima ./2_ccs/NOZ_merge.ccs.bam primer.fasta NOZ_merge.demux.bam --isoseq --peek-guess
lima ./2_ccs/ZJU0430_merge.ccs.bam primer.fasta ZJU0430_merge.demux.bam --isoseq --peek-guess
lima ./2_ccs/OCUG.ccs.bam primer.fasta OCUG.demux.bam --isoseq --peek-guess

lima ./2_ccs/normalEpi.ccs.bam   normalEpi.demux.bam --isoseq --peek-guess

path="/picb/sysgenomics3/project/GBC/PBRNA/ccs"
lima ${path}/p1/p1.ccs.bam  ./lyb7-9-10-27-PBRNA-barcode-1.fasta  p1.demux.bam --isoseq --peek-guess

path="/picb/sysgenomics3/project/GBC/PBRNA/ccs"
lima ${path}/p2/p2.ccs.bam  ./lyb7-9-10-27-PBRNA-barcode-2.fasta  p2.demux.bam --isoseq --peek-guess

# 03 refine
# Trimming of poly(A) tails
# Rapid concatemer identification and removal
# If your sample has poly(A) tails, use --require-polya. 
# This filters for FL reads that have a poly(A) tail with at least 20 base pairs (--min-polya-length) and removes identified tail

cd 4_refine
##for cell data
isoseq3 refine --require-polya 3_lima/GBC_SD_merge.demux.GBC_5p--GBC_3p.bam primer.fasta GBC_SD_merge.flnc.bam
isoseq3 refine --require-polya 3_lima/L2F7_merge.demux.L2F7_5p--L2F7_3p.bam primer.fasta L2F7_merge.flnc.bam
isoseq3 refine --require-polya 3_lima/NOZ_merge.demux.NOZ_5p--NOZ_3p.bam primer.fasta NOZ_merge.flnc.bam
isoseq3 refine --require-polya 3_lima/ZJU0430_merge.demux.ZJU0430_5p--ZJU0430_3p.bam primer.fasta ZJU0430_merge.flnc.bam
isoseq3 refine --require-polya 3_lima/OCUG.demux.OCUG_5p--OCUG_3p.bam primer.fasta OCUG.flnc.bam 
##for normal tissue data
isoseq3 refine --require-polya 3_lima/normalEpi.demux.Clontech_5p--NGE_31_3p.bam NGE_31.flnc.bam
isoseq3 refine --require-polya 3_lima/normalEpi.demux.Clontech_5p--NGE_32_3p.bam NGE_32.flnc.bam
isoseq3 refine --require-polya 3_lima/normalEpi.demux.Clontech_5p--NGE_34_3p.bam NGE_34.flnc.bam
isoseq3 refine --require-polya 3_lima/normalEpi.demux.Clontech_5p--OD_GB_E_3p.bam NGE_34.flnc.bam

## 04 dataset 
## merge all  flnc.bam
dataset=/software/smrt/smrtlink/smrtcmds/bin/dataset 
${dataset} create --type ConsensusReadSet combined_demux.consensusreadset_rmPT.xml \
  GBC_SD_merge.flnc.bam  L2F7_merge.flnc.bam NOZ_merge.flnc.bam ZJU0430_merge.flnc.bam OCUG.flnc.bam \
  NGE_31.flnc.bam NGE_32.flnc.bam NGE_34.flnc.bam OD_GB_E.flnc.bam \
  T10.flnc.bam T22.flnc.bam T25.flnc.bam T26.flnc.bam \
  T27.flnc.bam T33.flnc.bam T7.flnc.bam T9.flnc.bam

## 05 cluster
##Cluster FLNC reads and generate transcripts (FLNC to TRANSCRIPTS)
isoseq3 cluster combined_demux.consensusreadset_rmPT.xml polish_rmPT.bam --verbose --use-qvs

## 06 pbmm2
## Map the hq to ref
ref_fa=GRCh38.primary_assembly.genome.fa
pbmm2 align --preset ISOSEQ --sort polish_rmPT.hq.bam ${ref_fa} GBC_8s_rmPT.mapped.hq.bam

## 07 Collapse into unique isoforms
isoseq3 collapse GBC_8s_rmPT.mapped.hq.bam GBC_8s_rmPT.collapse.gff --do-not-collapse-extra-5exons

## 08 Cupcakes pipeline

get_abundance_post_collapse.py GBC_8s_rmPT.collapse polish_rmPT.cluster_report.csv
##note:polish_rmPT.cluster_report.csv from isoseq3 cluster step
filter_by_count.py --min_count 3 --dun_use_group_count GBC_8s_rmPT.collapse
filter_away_subset.py GBC_8s_rmPT.collapse.min_fl_3

## 09 SQANTI3
sqanti3="/software/SQANTI3-5.1.1/sqanti3_qc.py"

file_path="/picb/sysgenomics/gaoli/cell_data/cell_pb_trans_merge/8_isoseq_collapse_rmPT"
ref_gtf=gencode.v39.annotation.gtf
ref_fa=GRCh38.primary_assembly.genome.fa

polyA_motif_list=/software/SQANTI3-5.1.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt
cage_peak=/software/SQANTI3-5.1.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
out_dir=9_squanti_rmPT
aboundance_file=GBC_8s_rmPT.collapse.min_fl_3.abundance.txt

python ${sqanti3} ${file_path}/GBC_8s_rmPT.collapse.min_fl_3.filtered.gff \
               ${ref_gtf} ${ref_fa} \
              --CAGE_peak ${cage_peak} \
              --polyA_motif_list ${polyA_motif_list} \
              -o GBC \
              -d ${out_dir} \
              -fl ${aboundance_file} \
              --cpus 10 \
              --chunks 6 \
              --report both \
              --isoAnnotLite \
              --saturation \
              --short_reads GBC.RNA.fq.fofn
