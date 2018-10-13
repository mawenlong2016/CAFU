## CAFU
CAFU is a Galaxy-based bioinformatics framework for comprehensive assembly and functional annotation of unmapped RNA-seq data from single- and mixed-species samples which integrates plenty of existing NGS analytical tools and our developed programs, and features an easy-to-use interface to manage, manipulate and most importantly, explore large-scale unmapped reads. Besides the common process of reads cleansing, reads mapping, unmapped reads generation and novel transcription assembly, CAFU optionally offers the multiple-level evidence analysis of assembled transcripts, the sequence and expression characteristics of assembled transcripts, and the functional exploration of assembled transcripts through gene co-expression analysis and genome-wide association analysis. Taking the advantages of machine learning (ML) technologies, CAFU also effectively addresses the challenge of classifying species-specific transcript assembled using unmapped reads from mixed-species samples. The CAFU project is hosted on GitHub (https://github.com/cma2015/CAFU), publically available via an Amazon Elastic Cloud disk image and a standardized Docker container.

## Overview of functional modules in CAFU
### Extraction of unmapped reads

<table>
    <tr>
        <td>Functions</td>
        <td>Applications</td>
        <td>Input files</td>
        <td>Main output files</td>
        <td>Programs</td>
        <td>References</td>
   </tr>
    <tr>
        <td>Quality Control</td>
        <td>Examine the quality of RNA-Seq data</td>
        <td>Quality examination reports (html)</td>
        <td>Main output files</td>
        <td>FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)</td>
        <td>[1]</td>
   </tr>
    <tr>
        <td rowspan="3">Trim Raw Reads</td>
        <td rowspan="3">Trim poly-A/T and low-quality reads</td>
        <td rowspan="3">Raw RNA-seq data</td>
        <td rowspan="3">High-quality RNA-seq data (Fastq)</td>
   </tr>
    <tr>
        <td>Fqtrim (version 0.9.7; https://ccb.jhu.edu/software/fqtrim/)</td>
        <td>[2]</td>
    </tr>
    <tr>
        <td>Trimmomatic (version 0.36; http://www.usadellab.org/cms/?page=trimmomatic)</td>
        <td>[3]</td>
    </tr>
    <tr>
        <td rowspan="4">Generate Unmapped Reads</td>
        <td rowspan="4">Align trimmed reads and obtain unmapped reads</td>
        <td rowspan="4">High-quality RNA-seq data; Reference sequences of corresponding species</td>
        <td rowspan="4">Alignment results (BAM); Unmapped reads (Fastq)</td>
    </tr>
    <tr>
        <td>HISAT2 (version 2.1.0; https://ccb.jhu.edu/software/hisat2/index.shtml)</td>
        <td>[4]</td>
    </tr>
    <tr>
        <td>SAMTools (version 1.8; http://samtools.sourceforge.net/)</td>
        <td>[5]</td>
    </tr>
    <tr>
        <td>BEDTools (version 2.27.0; http://bedtools.readthedocs.io/en/latest/)</td>
        <td>[6]</td>
    </tr>
 </table>
