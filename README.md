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


## How to use CAFU

- More details about how to use CAFU are presented at https://github.com/cma2015/CAFU/blob/master/tutorial/User_manual.md
- Test datasets can de downloaded from https://github.com/cma2015/CAFU/tree/master/test_data


## How to access help
* If users encounter any bugs or issues, feel free to leave a message at Github [issues](<https://github.com/cma2015/CAFU/issues>). We will try our best to deal with all issues as soon as possible.
* In addition, if any suggestions are available, feel free to contact: __Siyuan Chen__ <chenzhuod@gmail.com> or __Jingjing Zhai__ <zhaijingjing603@gmail.com> 

## How to cite this work
Siyuan Chen#, Chengzhi Ren#, Jingjing Zhai#, Jiantao Yu#, Xuyang Zhao, Zelong Li, Ting Zhang, Wenlong Ma, Zhaoxue Han, Chuang Ma, CAFU: A Galaxy framework for exploring unmapped RNA-Seq data (submitted to Briefings in Bioinformatics).

