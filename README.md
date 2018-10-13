[![Docker Repository on Quay](https://quay.io/repository/bgruening/galaxy-rna-workbench/status "Docker Repository on Quay")](https://hub.docker.com/r/malab/cafu/)

## CAFU
- CAFU is a Galaxy-based bioinformatics framework for comprehensive assembly and functional annotation of unmapped RNA-seq data from single- and mixed-species samples which integrates plenty of existing NGS analytical tools and our developed programs, and features an easy-to-use interface to manage, manipulate and most importantly, explore large-scale unmapped reads. Besides the common process of reads cleansing, reads mapping, unmapped reads generation and novel transcription assembly, CAFU optionally offers the multiple-level evidence analysis of assembled transcripts, the sequence and expression characteristics of assembled transcripts, and the functional exploration of assembled transcripts through gene co-expression analysis and genome-wide association analysis. Taking the advantages of machine learning (ML) technologies, CAFU also effectively addresses the challenge of classifying species-specific transcript assembled using unmapped reads from mixed-species samples. The CAFU project is hosted on GitHub(https://github.com/cma2015/CAFU) and can be accessed from 210.27.96.13:4001. In addition, in order to enable larger-scale analysis, we also provided a standardized Docker image: [CAFU Docker image](https://hub.docker.com/r/malab/cafu/).

## Overview of functional modules in CAFU
### Extraction of unmapped reads

<table>
    <tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
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
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
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


### Evidence support of assembled transcripts
<table>
<tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
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


### Species assignment of assembled transcripts
<table>
<tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
</tr>
<tr>
        <td>SAT</td>
        <td>Machine learning-based prediction of the original species of assembled transcripts</td>
        <td>CDSs of two species; Assembled transcripts</td>
        <td>SAT score</td>
        <td>In-house scripts</td>
        <td>Our study</td>
</tr>
</table>



### Sequence characterization of assembled transcripts
<table>
<tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
</tr>
<tr>
        <td>Characterize Nucleic-acid Feature</td>
        <td>Explore the similarity between assembled and annotated transcripts in terms of the distribution of transcript length, exon length, G+C content</td>
        <td>Confident transcripts; Reference transcripts</td>
        <td>Diagnostic plots (barplot/density)</td>
        <td>In-house scripts</td>
        <td>Our study</td>
</tr>
<tr>
        <td>Characterize Amino-acid Feature</td>
        <td>Explore the similarity between assembled and annotated transcripts in terms of the distribution of amino acid-based features used in SAT</td>
        <td>Confident transcripts; Reference transcripts</td>
        <td>Diagnostic plots (barplot/density)</td>
        <td>In-house scripts</td>
        <td>Our study</td>
</tr>
</table>


### Expression profiles of assembled transcripts
<table>
<tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
</tr>
<tr>
        <td>Analysis Condition Specificity</td>
        <td>Identify condition-specifically expressed transcripts</td>
        <td>All transcripts; High-quality RNA-seq data; Sample information</td>
        <td>Condition-specific table; Diagnostic plots(heatmap)</td>
        <td>In-house scripts</td>
        <td>Our study</td>
</tr>
<tr>
        <td>Analysis Heterogeneous</td>
        <td>Identify stablely expressed transcripts</td>
        <td>Expression matrix; Sample information</td>
        <td>Gini coefficient table; Diagnostic plot (dotplot)</td>
        <td>In-house scripts</td>
        <td>Our study</td>
</tr>
<tr>
        <td rowspan="2">Analysis Differential Expression</td>
        <td rowspan="2">Identify differentially expressed transcripts</td>
        <td rowspan="2">All transcripts; High-quality RNA-seq data; Sample information</td>
        <td rowspan="2">DE analysis table; Diagnostic plot (Volcano plot; Venn-daragram )</td>
        <td>RSEM (version 1.3.0; https://deweylab.github.io/RSEM/) ; Bowtie2 (version 2.3.4.1; http://bowtie-bio.sourceforge.net/index.shtml)</td>
        <td>[10, 11]</td>
</tr>
<tr>
        <td>RSEM (version 1.3.0; https://deweylab.github.io/RSEM/) ; Bowtie2 (version 2.3.4.1; http://bowtie-bio.sourceforge.net/index.shtml)</td>
        <td>[15]</td>
</tr>
</table>


### Function annotation of assembled transcripts
<table>
<tr>
        <td font-weight:bold>Functions</td>
        <td font-weight:bold>Applications</td>
        <td font-weight:bold>Input files</td>
        <td font-weight:bold>Main output files</td>
        <td font-weight:bold>Programs</td>
        <td font-weight:bold>References</td>
</tr>
<tr>
        <td rowspan="2">Co-expression and GO</td>
        <td rowspan="2">Invesitgate the putative function of assembled transcripts through co-expression network and GO enrichment analysis</td>
        <td rowspan="2">Invesitgate the putative function of assembled transcripts through co-expression network and GO enrichment analysis</td>
        <td rowspan="2">Co-expression result tables; Diagnostic plot (Dendrogram); GO analysis table</td>
        <td>WGCNA (version 1.63; https://cran.r-project.org/web/packages/WGCNA/index.html)</td>
        <td>[16]</td>
</tr>
<tr>
        <td>topGO (version 3.7; https://bioconductor.org/packages/release/bioc/html/topGO.html)</td>
        <td>[17]</td>
</tr>
</table>

## CAFU Docker image installation
### Docker installation and start
#### For Windows (Test on Windows 10 Enterprise version):
* Download [Docker](<https://download.docker.com/win/stable/Docker%20for%20Windows%20Installer.exe>) for windows </br>
* Double click the EXE file to open it;
* Follow the wizard instruction and complete installation;
* Search docker, select ___Docker for Windows___ in the search results and clickit.
#### For Mac OS X (Test on macOS Sierra version 10.12.6 and macOS High Sierra version 10.13.3):
* Download [Docker](<https://download.docker.com/mac/stable/Docker.dmg>) for Mac os <br>
* Double click the DMG file to open it;
* Drag the docker into Applications and complete installation;
* Start docker from Launchpad by click it.
#### For Ubuntu (Test on Ubuntu 14.04 LTS and Ubuntu 16.04 LTS):
* Go to [Docker](<https://download.docker.com/linux/ubuntu/dists/>), choose your Ubuntuversion, browse to ___pool/stable___ and choose ___amd64, armhf, ppc64el or s390x.____ Download the ___DEB___ file for the Docker version you want to install;
* Install Docker, supposing that the DEB file is download into following path:___"/home/docker-ce<version-XXX>~ubuntu_amd64.deb"___ </br>
```bash
$ sudo dpkg -i /home/docker-ce<version-XXX>~ubuntu_amd64.deb      
$ sudo apt-get install -f
```
 ### Verify if Docker is installed correctly
----------------------------------------
   Once Docker installation is completed, we can run ____hello-world____ image to verify if Docker is installed correctly. Open terminal in Mac OS X and Linux operating system and open CMD for Windows operating system, then type the following command:
```bash
$ docker run hello-world
```
   **<font color =red>Note</font>:** root permission is required for Linux operating system.
   **<font color =red>Note</font>:** considering that differences between different computers may exist, please refer to [official installation manual](https://docs.docker.com/install) if instructions above don’t work.

### CAFU installation from Docker Hub
--------------------------------
  For Mac OS X and Linux operating systems, open the terminal, for Windows operating system, open CMD. Typing the following command:
```bash
# Pull CAFU from Docker Hub
$ docker pull malab/cafu
```

### Qucikly start
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
