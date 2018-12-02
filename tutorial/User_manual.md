### Brief introduction
- CAFU is a Galaxy-based bioinformatics framework for comprehensive assembly and functional annotation of unmapped RNA-seq data from single- and mixed-species samples which integrates plenty of existing NGS analytical tools and our developed programs, and features an easy-to-use interface to manage, manipulate and most importantly, explore large-scale unmapped reads. Besides the common process of reads cleansing, reads mapping, unmapped reads generation and novel transcription assembly, CAFU optionally offers the multiple-level evidence analysis of assembled transcripts, the sequence and expression characteristics of assembled transcripts, and the functional exploration of assembled transcripts through gene co-expression analysis and genome-wide association analysis. Taking the advantages of machine learning (ML) technologies, CAFU also effectively addresses the challenge of classifying species-specific transcript assembled using unmapped reads from mixed-species samples. The CAFU project is hosted on GitHub(https://github.com/cma2015/CAFU) and can be accessed from http://bioinfo.nwafu.edu.cn:4001. In addition, in order to enable large-scale analysis, we also provided a standardized Docker image: [CAFU Docker image](https://hub.docker.com/r/malab/cafu/).


### CAFU Docker image installation
- Step 1: [Docker installation](https://github.com/cma2015/CAFU/blob/master/tutorial/Docker_installation.md)
- Step 2: CAFU installation from Docker Hub
  ```bash
  # Pull latest version of CAFU from Docker Hub
  $ docker pull malab/cafu
  ```
- Step 3: Qucikly start
  ```bash
  $ docker run -it -p 80:80 malab/cafu bash
  $ cd /home/galaxy
  $ bash run.sh
  ```
  Then you can access CAFU instance via http://localhost:80


### Upload data
#### Download CAFU test data
- Download test data from [CAFU GitHub project](https://github.com/cma2015/CAFU). Click **Clone or download** (see figure below), and download the ZIP compressed files into your local device and then uncompress it. 

  ![Download data](https://github.com/cma2015/CAFU/blob/master/CAFU_images/1.png)


- For users who installed [Git](https://gist.github.com/derhuerst/1b15ff4652a867391f03), the following command can be used to download CAFU project to local device.
  ```bash
  git clone https://github.com/cma2015/CAFU.git
  ```
#### Upload regular file
- Click **Get Data** in the homepage (see figure below) of CAFU to upload files.


  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)


  And then you will see the following interface:

  ![](https://github.com/cma2015/CAFU/blob/master/CAFU_images/3.png)


  Next, Click the button **Choose local file** and select a file you would like to upload (e.g. upload the file ```mapping_info``` in the directory ```/your directory/CAFU/test_data/SE RNA-Seq/```), you will see the following interface:

  
  ![Upload regular file](https://github.com/cma2015/CAFU/blob/master/CAFU_images/4.png)

  
  Then click **Start** to upload file.

#### Upload a collection
- Similar to **Upload regular file**, click **Get Data** first (see figure below):


  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)

  
  And then you will see the following interface:
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/5.png)

  
  Afterwards, select a list files to upload as a collection (e.g. upload all files with ZIP suffix in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/```):

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/6.png)


  Click **Start** to upload, after finishing uploading, click **Build** (see figure below):

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/7.png)

  
  Then enter a name for your collection and click **Create list** to finish.

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/8.png)


### UNMAPPED READ EXTRACTION
In this section, we provide an example for each module to show how to use to perform unmapped reads extraction.
- **Quality control**

  In this module, we implemented FastQC (Andrews *et al*., 2010) to enable users to perform quality control. In this tutorial, we will use a list of single-end RNA-Seq collection file (located in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/```) to perform quality control.
  
  To run this module correctly, upload the files in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/``` with **ZIP** suffix as a collection named as ```SE-RNA-Seq```（see section **Upload a collection** to see how to upload a list of RNA-Seq datasets as a collection).

  
  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/9.png)

  
  Finally, click **Execute** to start performing quality control.
  Once the analysis is finished, a basic text and a HTML output file will be returned. The HTML output contains:

  - Basic Statistics
  - Per base sequence quality
  - Per sequence quality scores
  - Per base sequence content
  - Per base GC content
  - Per sequence GC content
  - Per base N content
  - Sequence Length Distribution
  - Sequence Duplication Levels
  - Overrepresented sequences
  - Kmer Content

  You can access the output by click results in the right history bar (see figure below).


  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/10.png)


- **Trim raw reads**

  In this function, poly-A/T is firstly trimmmed using fqtrim (Pertea, 2015), and then high-quality reads (e.g. score > 20) are retained by Trimmomatic (Bolger et al., 2014). The reads less than 20bp (in default) are discarded. 
  
  Here, we used the same collection (```SE-RNA-Seq```) with last step (**Quality control**) to trim raw reads.


  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/11.png)


  For each FASTQ formatted file, CAFU returns two outputs including:

  - ```fqtrim_DATA_ID.fastq```: poly-A/T trimmed RNA-Seq using fqtrim.

  - ```Trimmed_DATA_ID.fastq```: High-quality RNA-Seq generated by Trimmomatic.

- **Extract Unmapped Reads**

  This function integrates several bioinformatic softwares including HISAT2 (Kim et al., 2015), SAMTools (Li et al., 2009), and BEDTools (Quinlan et al., 2010) to extract unmapped reads. Firstly, HISAT2 is used to align user-defined high-quality reads to reference genome. Secondly, SAMTools is used to extract unmapped alignments in SAM format from alignments generated by HISAT2 with parameter "-f 4". Finally, BEDTools is used to convert the unmapped alignments in SAM format to FASTQ format. In addition, CAFU also supports dual RNA-Seq analysis by aligning RNA-Seq data to multiple reference genomes.

  To run this module, three inputs (see figure below) are required. Here, we still  use single-end RNA-Seq files (upload these files as a collection named as ```SE-RNA-Seq```) in directory ```/your directory/CAFU/test_data/SE RNA-Seq/``` to show how to use this module.
 
  **Inputs description:**
  
  **Input 1:** A **collection** of reference genome. Upload files (```maize.fa.zip```) in directory ```/your directory/CAFU/test_data/genomes/``` as a collection named as ```Ref_Genome``` (user-defined name). See section **Upload collection file** to learn how to upload a collection.
  
  **Input 2:** A **collection** of trimmed single-end RNA-Seq files (```Trimmed_DATA_ID.fastq```) generated from **Trim Raw Reads**.
  
  **Input 3:** A **regular** file containing mapping-operation information. Upload files in directory ```/your directory/CAFU/test_data/SE RNA-Seq/mapping_info``` as a regular file.
  
  **Note:** **Input 3** is a semicolon seperated matrix which contains two columns. The first column contains RNA-Seq ID in each experiment. The second column is the corresponding reference genome ID for each experiment.
  ```bash
  SRR2144382;maize
  SRR2144383;maize
  SRR2144384;maize
  SRR2144385;maize
  SRR2144410;maize
  SRR2144411;maize
  SRR2144442;maize
  SRR2144443;maize
  ```
 
  ![Extrat unmapped reads](https://github.com/cma2015/CAFU/blob/master/CAFU_images/12.png)


  Then click **Execute** to start extracting unmapped reads using CAFU. After finishing this process, the final unmapped reads named as ```all_unmapped_reads.fastq``` (see figure below) will be returned.


  ![unmapped reads output](https://github.com/cma2015/CAFU/blob/master/CAFU_images/13.png)


- **Remove Contamination**

  Considering that sequences obtained from impure nucleic acid preparations may contain DNA from sources other than the sample. Those sequence contaminations are a serious concern to the quality of the data used for downstream analysis. So in this module, potential contamination sequences are removed using Deconseq (Schmieder *et al*., 2011) with user-defined coverage and identity (e.g., 0.95) by aligning input reads to a contamination database. In the current version of CAFU, 3,529 bacterial and 81 viral reference genomes (downloaded from NCBI on 2018/11/05) are provided as default, however, user-defined contamination database is also supported. 

  ![remove contamination](https://github.com/cma2015/CAFU/blob/master/CAFU_images/17.png)


  Then the clean reads with FASTQ format will be returned.


###  *DE NOVO* TRANSCRIPT ASSEMBLY OF UNMAPPED READS
- **Assemble Unmapped Reads**

  In this module, three steps including ***de novo* assembly of unmapped reads**, **remove redundancy of transcript fragments**, and **re-assembly transcript fragments** will be sequentially performed to assemble unmapped reads. In this tutorial, we will use the unmapped reads generated by module ```Extract Unmapped Reads``` as the input (see figure below).

  
  ![unmapped reads output](https://github.com/cma2015/CAFU/blob/master/CAFU_images/14.png)
  

  Then assembled transcripts named as ```Unmapped_reads_de_novo_assembled_transcripts```  will be returned (see figure below).


  ![assembled transcripts](https://github.com/cma2015/CAFU/blob/master/CAFU_images/15.png)


### EVIDENCE SUPPORT OF ASSEMBLED TRANSCRIPTS
- **Expression-level Evidence**

  This module allows users to eliminate assembled transcripts with low read coverage and/or low expression abundance, which are likely assembly artifacts. RNA-Seq reads used for Assemble Unmapped Reads are mapped with newly assembled transcripts and reference transcripts by using bowtie2 (Langmead et al., 2012). CAFU outputs the read coverage of assembled transcripts at single-base resolution using BEDTools (Quinlan et al., 2010), and estimates the expression abundance of all transcripts in terms of FPKM (Fragments Per Kilobase Million) using RSEM (Li et al., 2011). Assembled transcripts with low read coverage (e.g., less than 10) and/or low expression (e.g., FPKM less than 1) in the majority of samples (e.g., 80%) are discarded.

  The only required input in this module is the newly assembled transcripts with FASTA format from unmapped reads. Here, we use the file ```/your directory/CAFU/test_data/others/assembled_transcript.fasta``` to show CAFU's usage.


  ![expression-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/18.png)
  
  
  Then three files will be returned:

  
  **Output 1**: ```Read coverage of each transcript in each sample```, a matrix whose rows represent transcripts, and columns represent read coverage.

  **Output 2**: ```Expression abundance of each transcripts in each sample```, a matrix whose rows represent transcripts, and columns represent expression abundance.

  **Output 3**: ```Confident transcript ID by Expression-level Evidence```, the common transcript ID from Confident transcript ID filtered by read coverage results and expression abundance.


- **Genome-level Evidence**
  
  This function can be used to identify de novo-assembled transcripts missing from the existing genome annotation. All the assembled transcripts are aligned to the reference genome sequences of the species of interests using GMAP (Wu *et al*., 2005), and selects the best genomic matches with high identity (e.g., ≥ 95%) and coverage (e.g., ≥ 95%). Users can also eliminate assembled transcripts with no introns, which could represent either noise or pseudogenes.
  
  To run this module, at least three inputs are required including:

  **Input 1**: Reference genome sequences of the species of interest (input as a data collection).

  **Input 2**: Assembled transcript sequences generated from module **Assemble Unmapped Reads**.

  **Input 3**: Genome annotation (GFF/GTF format) of corresponding species with RNA-Seq samples.

  **Input 4**: A character indicating the reference genome name (e.g. maize).


  ![genome-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/19.png)


  Then four outputs will be returned:

  **Output 1**: ```Integrated GMAP results of newly assembled transcripts against all reference genome sequences```: GMAP alignment results (coverage and identity) of each assembled transcript against all reference genome sequences.

  **Output 2**: ```Confitent transcript information```: Confident transcript information filtered by high coverage and identity.

  **Output 3**: ```The same/similar-intron transcript ID```: Assembled transcript ID which possess the same/similar intron with corresponding species reference transcripts.

  **Output 4**: ```Novel transcript ID```: Novel transcripts missing in the existing genome annotation.

- **Transcript-level Evidence**

  This function can be used to select assembled transcripts with high similarity comparing to other well-annotated transcripts, such as full-length transcripts generated from single-molecule real-time sequencing and/or high-quality transcripts annotated in closely related species. After aligning assembled transcripts with other well-annotated transcripts with GAMP (Wu *et al*., 2005), CAFU outputs the best transcript alignments with high identity (e.g., ≥ 95%) and coverage (e.g., ≥ 95%).
  
  To run this module, at least two inputs are required including:

  **Input 1**: Reference sequence of well-annotated transcripts, such as full-length transcripts generated from single-molecule real-time sequencing and/or high-quality transcripts annotated in closely related species.

  **Input 2**: Assembled transcript sequences generated from module **Assemble Unmapped Reads**.


   ![genome-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/22.png) 

   
   Then three outputs will be returned:

   **Output 1**: ```Integrated GMAP results of newly assembled transcripts against all reference transcript sequences```, GMAP alignment results (coverage and identity) of each assembled transcript against all other transcript sequences.

   **Output 2**: ```Confident transcript information```, Confident transcript information filtered by high coverage and identity.

   **Output 3**: ```Confident transcript ID```, IDs of transcripts that prossess high similarity comparing to other well-annotated transcripts.

- **Protein-level Evidence**

  In this module, coding potential evidence of transcripts is fistly evaluated using CPC2 (Kang *et al*., 2017). Then for coding transcripts, Pfam (Finn *et al*., 2014) will be used to identify the protein families. 
  
  Here, we use the file ```/your directory/CAFU/test_data/others/assembled_transcript.fasta``` to execute this module.

  Then three outputs will be returned:

  **Output 1**: CPC2 output. A tab seperated matrix contains seven columns. Each column shows the sequence ID, putative peptide length, Fickett score, isoelectric point, the integrity of the orf, coding probability and the coding/noncoding classification label. More details about this output can be seen from [CPC2 official website](http://cpc2.cbi.pku.edu.cn/help.php). 

  **Output 2**: Confident transcript ID. 

  **Output 3**: A tab seperated matrix contains transcript ID, alignment start, alignment end, envelope start, envelope end, Hmm access, Hmm name, Type of domain, Hmm start, Hmm end, Hmm length, Bit score, E-value, Significance, Clan, etc.


### SPECIES ASSIGNMENT OF ASSEMBLED TRANSCRIPTS
- **Species Assignment of Transcripts**

  SAT (Species Assignment of Transcripts) is a machine learning-based toolkit used for species assignment of transcritps assembled using unmapped reads from mixed-species (eg., pathogen-host) samples. 
  
  In this tutorial, we will show how to use SAT using all the FASTA format files in directory ```/your directory/CAFU/test_data/SAT/```. To run SAT, at least three files are required including:

  **Input 1**: The positive coding sequences (CDS) with FASTA format.  

  **Input 2**: The negative coding sequences (CDS) with FASTA format.

  **Input 3**: The coding sequences (CDS) with FASTA format assembled using unmapped reads from mixed-species (e.g., pathogen-host) 

  ![assembled transcripts](https://github.com/cma2015/CAFU/blob/master/CAFU_images/16.png)

  Then click **Execute** to run this module, then four outputs will be returned:

  **Output 1**: ```Probabilistic score of each transcript```, A probabilistic score of each transcript, the higher, the more probable being positive.

  **Output 2**: ```Eight commonly used measures under specidied threshold```, A bar plot evaluating the measures (Sn, Sp, Pr, Acc, MCC, Fscore, AUC and AUPR) under specified threshold. 

  **Output 3**: ```The presicion recall curves in k-fold cross validation```, The PR curve of k-fold cross-validation. 

  **Output 4**: ```The receiver operating curves in k-fold cross validation```, The ROC curve of k-fold cross-validation. 

### SEQUENCE CHARACTERIZATION OF ASSEMBLED TRANSCRIPTS

- **Characterize Nucleic-acid Feature**

  This module allows users to explore the nucleic-acid similarity between assembled and reference transcripts in terms of the distribution of transcript length and G+C content.

  In this module, two features can be analyzed using this module (see figure below). To run this module, two inputs are required including:

  **Input 1**: The sequences of assembled transcripts (or confident/novel transcripts) derived from unmapped reads. 

  **Input 2**: The sequences of transcripts from the existing genome annotation. 


  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/21.png)

  Then three outputs will be returned:

  - **For Transcript length**

    **Assembled transcript length**: Length of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    **Reference transcript length**: Length of transcripts from the existing genome annotation.

    **Length distribution comparison**: Length distribution comparison between assembled transcripts and reference transcript.

  - **For GC content**


    **Assembled transcript GC content**: GC content of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    **Reference transcript GC content**: GC content of transcripts from the existing genome annotation.

    **Length distribution comparison**: GC content distribution comparison between assembled transcripts and reference transcript.


- **Characterize Amino-acid Feature**

  This module allows users to explore the amino-acid features similarity used in SAT between assembled and reference transcripts.

  Here, we take an example (see figure below) to show how to use this module to compare k-mer frequency of assembled and reference transcripts.

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/23.png)
  
  The outputs contain:
  
  **Note**: The detailed feature descriptions are available at http://bioinformatics.hitsz.edu.cn/BioSeq-Analysis/
  
  - **For Kmer**

    ```Assembled transcript K-mer (k = 1)```: K-mer frequency of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript K-mer (k = 1)```: K-mer frequency of transcripts from the existing genome annotation.

    ```Assembled transcript K-mer (k = 2)```: K-mer frequency of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.
    ```Reference transcript K-mer (k = 2)```: K-mer frequency of transcripts from the existing genome annotation.
  - **For DR**

    ```Assembled transcript DR```: Distance-based residues encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript DR```: Distance-based residues encoding of transcripts from the existing genome annotation.
  - **For AC**

    ```Assembled transcript AC```: Auto-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript AC```: Auto-covariance encoding of transcripts from the existing genome annotation.
  - **For CC**

    ```Assembled transcript CC```: Cross-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads
    
    ```Reference transcript CC```: Cross-covariance encoding of transcripts from the existing genome annotation.
  - **For ACC**

    ```Assembled transcript ACC```: Auto-cross-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript ACC```: Auto-cross-covariance encoding of transcripts from the existing genome annotation.
  - **PDT**

    ```Assembled transcript PDT```: Kmer of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript PDT```: Kmer of transcripts from the existing genome annotation.

  - **For PC-PseAAC**

    ```Assembled transcript PC-PseAAC```: Physicochemical distance transformation encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript PC-PseAAC```: Physicochemical distance transformation encoding of transcripts from the existing genome annotation.
  - **For SC-PseAAC**

    ```Assembled transcript SC-PseAAC```: General series correlation pseudo amino acid composition encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript SC-PseAAC```: General series correlation pseudo amino acid composition encoding of transcripts from the existing genome annotation.

- **Detect Alternative Splicing Events**

  This function is used to detect alternative splicing events in assembled transcripts using an R pacaage SGSeq (Goldstein *et al*., 2016).

  The only input of this module is GFF/GTF annotation file, here we use the test data located in ```/your directory/CAFU/test_data/other/AS_test.gtf``` to run this module.

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/28.png)

  Then two outputs will be returned:

  **Output 1**: ```Alternative splicing results```, Alternative splicing events of all transcripts.

  **Output 2**: ```Read counts of all genome features in GTF file```, Read counts supported for all genome features in GTF file provided by users. 

### EXPRESSION PROFILES OF ASSEMBLED TRANSCRIPTS
- **Analyze Condition Specificity**

  This function identifies a set of transcripts highly expressed under different conditions. The condition specificity of a transcript for condition type T is defined using the formula described in (Ma et al., 2014).

  Here, we use the test data ```assembled_transcript_expression, RNA-Seq_sample_information``` in the folder ```/your directory/CAFU/test_data/others/``` to show its usage (see figure below).

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/24.png)

   Then two outputs will be returned:

   ```Condition specific transcritps information```: Information of condition-specific transcripts, including condition-specific transcript ID and corresponding specifically expressed condition.

   ```Expression heatmap condition-specific transcritps```: Digital plot of condition-specific transcript expression.

- **Analyze Heterogeneous**

  This function examines the stability of each transcript using its expression values in all samples and Gini index (coefficient) according to (O'Hagan *et al*., 2017).

  To run this module, the only required input is the ```Transcript expression abundance matrix```, which is a tab seperated expression abundance matrix with the rows as transcripts and the columns as samples. Here, we still use the sample data ```assembled_transcript_expression, RNA-Seq_sample_information``` in the folder ```/your directory/CAFU/test_data/others/``` to show its usage (see figure below).

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/25.png)

  The two outputs will be returned:

  **Output 1**: ```Gini Coefficient of assembled transcripts```, A tab seperated matrix contains transcript ID, and corresponding Gini Coefficient.

  **Output 2**: ```Plot of Gini coefficient of assembled transcripts```: Scatter diagram of Gini coefficient of each assembled transcript.

- **Analyze Differential Expression**

  This function integrates RSEM (Li *et al*., 2011) and EBSeq (Leng *et al*., 2013) to identify differential expression transcripts.

  Then four outputs will be returned:

  ```Differential expression transcript information```: Information of differential expression (DE) transcripts filtered by user-specific criteria , including transcript ID, FDR, fold change, etc.

  ```Differential expression transcript ID```: IDs of DE transcripts filtered by differential expression criteria.

  ```Expression value of all transcripts```: Expression abundance matrix of **all transcripts** with the rows as transcripts and the columns as samples.

  ```Expression value of DE transcripts```: Expression abundance matrix of **DE transcripts**   with the rows as transcripts and the columns as samples.

### FUNCTION ANNOTATION OF ASSEMBLED TRANSCRIPTS

- Co-expression and GO

  In this module, co-expression network and GO enrichment analysis are used to annotate transcripts using "WGCNA" and "topGO", respectively.

  In this tutorial, we used the file located in ```/your directory/CAFU/test_data/differentially_expressed_transcript_expression``` to perform co-expression network and GO enrichment analysis (see figure below).

  ![Co-expression and GO](https://github.com/cma2015/CAFU/blob/master/CAFU_images/20.png)

  Then five results will be generated including:

   **Output 1**: ```Dendrograms and module colors```, Plot of dendrograms and module colors.

   **Output 2**: ```Edge File```, The edge results of co-expression network.

   **Output 3**: ```Node File```, The node results of co-expression network.

   **Output 4**: ```Hub transcript ID```, The hub transcript ID in each module of co-expression network.

   **Output 5**: ```GO results```, GO enrichment results of each module.

### OTHER TOOLS


- Remove Batch Effect

  This function can be used to remove batch effect using an R package sva (Leek *et al*., 2012).

  To run this module, there are two required inputs (see figure below) including:

  **Input 1**: ```Transcript expression abundance matrix```, Expression abundance matrix with the rows as transcripts and the columns as samples.

  **Input 2**: RNA-Seq sample batch effect information.
  
  - **Expample:** Number 1, 2, 3 represent different experiments. The batch information should be:
    ```bash
    SRR001;1
    SRR002;1
    SRR003;2
    SRR004;3
    SRR005;3
    SRR006;3
    ```

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/26.png)

  Then the corrected expression matrix named as ```Corrected transcript abundance matrix``` will be returned.


- Extract Sequences

  This function is provided for extracting transcript sequences. It includes three ways, one is extracting sub-sequences (novel transcripts) from a fasta file (assembled transcript sequences) according to target sequence ID (novel transcript ID); one is used for extracting transcript sequences from reference genome according to genome annotation file; the rest one is combine the former two ways to establish an all transcripts file.

  For **Extract sub-sequences according to sequence ID from a fasta file**, there are two required inputs including:
  
  **Input 1 **: ```All assembled transcript sequences``` (This file can be generated from the function **Assemble Unmapped Reads**).
  
  **Input 2 **: ```Novel transcript ID```: ID of novel transcripts. (This file can be generated from the module **iii) Evidence support of assembled transcripts**).
  
  Then the target transcript sequence file named as ```Sequences of novel transcripts``` will be returned.
  
  For **Extract transcripts according to GTF/GFF from reference genome sequence**, there are also two required inputs including:
  
  **Input 1 **: ```Reference genome sequence``` Sequences of reference genome including all chromosome..
  
  **Input 2 **: ```Reference genome annotation``` Corresponding reference genome annotation file with reference genome.
  
  Then the target transcript sequence file named as ```Sequences of reference transcripts``` will be returned.
  
  For **Establsh all transcript file (including novel and reference transcripts)**, the input files is the same as the former two ways.
  
  Then the target transcript sequence file named as ```Sequences of all transcripts``` will be returned.
  
