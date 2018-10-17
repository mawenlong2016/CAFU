[![Docker Repository on Quay](https://quay.io/repository/bgruening/galaxy-rna-workbench/status "Docker Repository on Quay")](https://hub.docker.com/r/malab/cafu/)


## CAFU
- CAFU is a Galaxy-based bioinformatics framework for comprehensive assembly and functional annotation of unmapped RNA-seq data from single- and mixed-species samples which integrates plenty of existing NGS analytical tools and our developed programs, and features an easy-to-use interface to manage, manipulate and most importantly, explore large-scale unmapped reads. Besides the common process of reads cleansing, reads mapping, unmapped reads generation and novel transcription assembly, CAFU optionally offers the multiple-level evidence analysis of assembled transcripts, the sequence and expression characteristics of assembled transcripts, and the functional exploration of assembled transcripts through gene co-expression analysis and genome-wide association analysis. Taking the advantages of machine learning (ML) technologies, CAFU also effectively addresses the challenge of classifying species-specific transcript assembled using unmapped reads from mixed-species samples. The CAFU project is hosted on GitHub(https://github.com/cma2015/CAFU) and can be accessed from http://bioinfo.nwafu.edu.cn:4001. In addition, in order to enable large-scale analysis, we also provided a standardized Docker image: [CAFU Docker image](https://hub.docker.com/r/malab/cafu/).

![CAFU](https://github.com/cma2015/CAFU/blob/master/CAFU_images/Overview%20of%20CAFU.png)

## Overview of functional modules in CAFU
- [**Extraction of unmapped reads**](https://github.com/cma2015/CAFU/blob/master/tutorial/Extraction_mapped_reads.md)
- [***De novo* transcript assembly of unmapped reads**](https://github.com/cma2015/CAFU/blob/master/tutorial/De_novo_transcript_assembly_of_unmapped_reads.md)
- [**Evidence support of assembled transcripts**](https://github.com/cma2015/CAFU/blob/master/tutorial/Evidence_support_of_assembled_transcripts.md)
- [**Species assignment of assembled transcripts**](https://github.com/cma2015/CAFU/blob/master/tutorial/SAT.md)
- [**Sequence characterization of assembled transcripts**](https://github.com/cma2015/CAFU/blob/master/tutorial/Sequence%20characterization%20of%20assembled%20transcripts.md)
- [**Expression profiles of assembled transcripts**](https://github.com/cma2015/CAFU/blob/master/tutorial/Expression%20profiles%20of%20assembled%20transcripts.md)
- [**Function annotation of assembled transcripts**](https://github.com/cma2015/CAFU/blob/master/tutorial/Function%20annotation%20of%20assembled%20transcripts.md)

## CAFU Docker image installation
- Step 1: [Docker installation](https://github.com/cma2015/CAFU/blob/master/tutorial/Docker_installation.md)
- Step 2: CAFU installation from Docker Hub
--------------------------------
  For Mac OS X and Linux operating systems, open the terminal, for Windows operating system, open CMD. Typing the following command:
```bash
# Pull CAFU from Docker Hub
$ docker pull malab/cafu
```

- Step 3: Qucikly start
```bash
$ docker run -it --net=host malab/cafu bash
$ cd /home/galaxy
$ bash run.sh
```
Then you can access CAFU instance via http://localhost:8080 


## References
1. Andrews S. FastQC: a quality control tool for high throughput sequence data 2010.
2. Pertea G. Fqtrim: v0. 9.4 release. Zenodo, 2015.
3. Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 2014;30:2114-2120.
4. Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods 2015;12:357-360.
5. Li H, Handsaker B, Wysoker A et al. The sequence alignment/map format and SAMtools. Bioinformatics 2009;25:2078-2079.
6. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 2010;26:841-842.
7. Grabherr MG, Haas BJ, Yassour M et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol 2011;29:644-652.
8. Li W, Godzik A. Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics 2006;22:1658-1659.
9. Huang X, Madan A. CAP3: A DNA sequence assembly program. Genome Res 1999;9:868-877.
10. Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 2011;12:323-338.
11. Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods 2012;9:357-359.
12. Wu TD, Watanabe CK. GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 2005;21:1859-1875.
13. Kang YJ, Yang DC, Kong L et al. CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features. Nucleic Acids Res 2017;45:W12-W16.
14. Finn RD, Bateman A, Clements J et al. Pfam: the protein families database. Nucleic Acids Res 2014;42:D222-D230.
15. Leng N, Dawson JA, Thomson JA et al. EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics 2013;29:1035-1043.
16. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008;9:559-571.
17. Alexa A, Rahnenfuhrer J. topGO: enrichment analysis for gene ontology. R package version 2010;2.
