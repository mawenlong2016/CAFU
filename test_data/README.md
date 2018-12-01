### SE RNA-Seq
- Eight single-end RNA-Seq datasets (SRR2144382.fastq.zip, SRR2144383.fastq.zip, SRR2144384.fastq.zip, SRR2144385.fastq.zip, SRR2144410.fastq.zip, SRR2144411.fastq.zip, SRR2144442.fastq.zip, SRR2144443.fastq.zip)
- DE_Info: A semicolon seperated matrix containing used for identifying differential expression transcripts (DE_Info).
- mapping_info: RNA-Seq information for aligning reads to reference genome.

### genomes
This folder includes three demo of reference genomes (maize.fa, Oryza.fa, Sorghum.fa) and an annotation of maize genomes (maize.gff3).

### Transcripts
This folder includes two well-annotated transcripts ("ref_trans.fasta", "ref_trans_1.fasta").

### SAT
This folder includes four files of sequences (negative_cds.fa, negative_cds.fa_BG.fa, positive_cds.fa, positive_cds.fa_BG.fa, test_cds.fa).

### others
This folder includes seven files (Assembled_transcript.fasta, Assembled_transcript_expression, AS_test.gtf, Cont_seq.fasta, Differentially_expressed_transcript_expression, RNA-Seq_sample_information).


### For each function:
- (1) Quality Control

  This function needs a file or a collection of RNA-Seq. User can use RNA-Seq data in the folder "RNA-Seq".

- (2) Trim Raw Reads

  The test data for this function include a collection of RNA-Seq from the folder "RNA-Seq". (upload as a collection, SEE USER MANNUAL)

- (3) Extract Unmapped Reads

  The test data for this function include "maize.fa" from the folder "Genomes", and trimmed RNA-Seq data generated from the function "Trim Raw Reads".

- (4) Remove Contamination

  This function needs these test data:
  
	  1) Unmapped reads generated from the function "Extract Unmapped Reads";
	  2) Contamination sequence: "Cont_seq.fasta" from the folder "Others";

- (5) Assemble Unmapped Reads

  The test data for this function is unmapped reads from the function either "Extract Unmapped Reads" or "Remove Contamination".

- (6) Expression-level Evidence

  This function needs the assembled transcripts generated from the function "Assemble Unmapped Reads".

- (7) Genome-level Evidence

  This function needs reference genomes from the folder "Genomes" (upload as a collection), and "Assembled_transcript.fasta" from the folder "Others".

- (8) Transcript-level Evidence

  This function needs well-annotated transcripts from the folder "Transcripts" (upload as a collection), and "Assembled_transcript.fasta" from the folder "Others".

- (9) Protein-level Evidence

  This function needs "Assembled_transcript.fasta" from the folder "Others".

- (10) Species Assignment of Transcripts
  
  The test data for this function are in the folder "SAT".

- (11) Characterize Nucleic-acid-based Features
  The test data for this function are "Assembled_transcript.fasta" from the folder "Others", and a well-annotated transcripts (such as "ref_trans.fasta") from the folder "Transcripts".

- (12) Characterize Amino-acid-based Features
  The test data for this function are the same as the function "Characterize Nucleic-acid Feature"

- (13) Detect Alternative Splicing Events
  
  This function needs the file "AS_test.gtf" in the folder "Others".

- (14) Analysis of Condition Specificity

  This function needs the file "Assembled_transcript_expression", and "RNA-Seq_sample_information" in the folder "Others".

- (15) Heterogeneous Analysis

  This function needs the file "Assembled_transcript_expression" in the folder "Others".

- (16) Differential Expression Analysis

  This function needs the file "Assembled_transcript.fasta" in the folder "Others" and RNA-Seq data from the folder "RNA-Seq" (upload as a collection)

- (17) Co-expression and Gene Ontology Enrichment Analysis

  This function needs the file "Differentially_expressed_transcript_expression".

- (18) Extract Sequences

  This function needs the file "maize.fa", "maize.gff3" and "Assembled_transcript.fasta" in the folder "Others" and a transcript ID list generated from "Analyze Differential Expression".

- (19) Remove Batch Effect 

  This function needs the file "Assembled_transcript_expression", and "Batch_information" in the folder "Others".
