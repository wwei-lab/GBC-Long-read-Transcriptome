# GBC-Long-read-Transcriptome
## A Novel ERBB2 Alternative Splicing Variant Induces Drug Resistance in Gallbladder Cancer Revealed by Long-read Transcriptome Sequencing

we employed the advanced long-read full-length transcriptome sequencing technology on GBC, normal tissues, and cell lines to establish a comprehensive transcriptomic atlas. One isoform with a novel exon was identified designated as ERRB2 i14e in multiple independent GBC datasets involving 179 GBC cases and the expression of ERRB2 i14e was associated with disease worse prognosis.

![image](https://github.com/user-attachments/assets/2998c0ad-3540-43b1-90a1-efcbc40c88a7)

## Pacbio CCS data processing

Refer to the official website and follow the step-by-step process
### isoseq3
https://isoseq.how/clustering/cli-workflow.html

### Followed by Cupcake：
https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step

### Followed by SQANTI3：
https://github.com/ConesaLab/SQANTI3

### For Interchromosomal fusion transcripts
We used FusionSeeker (v1.0.1) to identify interchromosomal fusion transcripts. Next, we selected the confident_genefusion_transcript_sequence.fa files and used BLAST to manually verify the correct fusion transcripts.


