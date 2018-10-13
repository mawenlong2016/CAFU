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
        <td rowspan="2">Trim Raw Reads</td>
        <td rowspan="2">Trim poly-A/T and low-quality reads</td>
        <td rowspan="2">Raw RNA-seq data</td>
        <td rowspan="2">High-quality RNA-seq data (Fastq)</td>
        <td>Fqtrim (version 0.9.7; https://ccb.jhu.edu/software/fqtrim/)</td>
        <td>[2]</td>
   </tr>
    <tr>
        <td>Trimmomatic (version 0.36; http://www.usadellab.org/cms/?page=trimmomatic)</td>
        <td>[3]</td>
    </tr>
    <tr>
        <td rowspan="3">Generate Unmapped Reads</td>
        <td rowspan="3">Align trimmed reads and obtain unmapped reads</td>
        <td rowspan="3">High-quality RNA-seq data; Reference sequences of corresponding species</td>
        <td rowspan="3">Alignment results (BAM); Unmapped reads (Fastq)</td>
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
 
### De novo transcript assembly of unmapped reads
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
    <td rowspan="3">Assemble Unmapped Reads</td>
    <td rowspan="3">De novo assemble unmapped reads, remove redundancy of transcript fragments and re-assemble transcript fragments</td>
    <td rowspan="3">Unmapped reads</td>
    <td rowspan="3">Assembled transcript sequences (Fasta)</td>
    <td>Trinity (version 2.2.0; https://github.com/trinityrnaseq/trinityrnaseq/wiki)</td>
    <td>[7]</td>
</tr>
<tr>
    <td>CD-HIT-EST (version 4.6.8; http://weizhongli-lab.org/cd-hit/)</td>
    <td>[8]</td>
</tr>
<tr>
    <td>CAP3 (version 12/21/07; http://seq.cs.iastate.edu/cap3.html)</td>
    <td>[9]</td>
</tr>
</table>



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
    <td rowspan="3">Expression-level Evidence</td>
    <td rowspan="3">Obtain confident transcripts according to read coverage and expression abundance</td>
    <td rowspan="3">All transcripts (combine reference transcript and assembled transcripts); High-quality RNA-seq data </td>
    <td rowspan="3">Expression matrix; Confident transcript ID</td>
    <td>RSEM (version 1.3.0; https://deweylab.github.io/RSEM/) ; Bowtie2 (version 2.3.4.1; http://bowtie-bio.sourceforge.net/index.shtml)</td>
    <td>[10，11]</td>
</tr>
<tr>
    <td>BEDTools (version 2.27.0; http://bedtools.readthedocs.io/en/latest/)</td>
    <td>[6]</td>
</tr>
<tr>
    <td>In-house scripts</td>
    <td>Our study</td>
</tr>
<tr>
    <td rowspan="2">Genome-level Evidence</td>
    <td rowspan="2">Obtain confident transcripts according to transcript-genome alignments</td>
    <td rowspan="2">Assembled transcripts; Reference sequences of corresonding and closely related speces</td>
    <td rowspan="2">Alignment results (PSL); Confident transcript ID</td>
    <td>GMAP (version 2015-09-29; https://github.com/juliangehring/GMAP-GSNAP)</td>
    <td>[12]</td>
</tr>
<tr>
    <td>In-house scripts</td>
    <td>Our study</td>
</tr>

<tr>
    <td rowspan="2">Transcript-level Evidence</td>
    <td rowspan="2">Obtain confident transcripts according to transcript-transcript alignments</td>
    <td rowspan="2">Assembled transcripts; Reference transcripts of corresonding and closely related speces</td>
    <td rowspan="2">Alignment results (PSL); Confident transcript ID</td>
    <td>GMAP (version 2015-09-29; https://github.com/juliangehring/GMAP-GSNAP)</td>
    <td>[12]</td>
</tr>
<tr>
    <td>In-house scripts</td>
    <td>Our study</td>
</tr>


<tr>
    <td rowspan="2">Protein-level Evidence</td>
    <td rowspan="2">Obtain confident transcripts according to the protein potential</td>
    <td rowspan="2">Assembled transcripts</td>
    <td rowspan="2">Confident transcripts ID; Protein potential assessment results; Domain/family results</td>
    <td>CPC2 (version 0.1; http://cpc2.cbi.pku.edu.cn/)</td>
    <td>[13]</td>
</tr>
<tr>
    <td>Pfam (version 31.0; https://pfam.xfam.org/)</td>
    <td>[14]</td>
</tr>

</table>
