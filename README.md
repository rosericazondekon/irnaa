Integrated RNA-seq Analysis Application (IRNAA)
================
Roseric AZONDEKON

February 14, 2019

<br/>

IRNAA is a shiny application built for the purpose of facilitating Differential Gene Expression (DGE) analysis. IRNAA integrates command line tools such as `salmon`, `fastQC`, `multiQC` with DGE `R` packages such as `DESeq2`, `edgeR`, and `limma-voom`, and hides all the abstractions and coding hurdles from the end user. IRNAA makes it possible for the user to either preprocess fastq files and import read counts into the application from the fastq preprocessed files, or directly import their read counts into the application, or import RNA-seq datasets directly from the The Cancer Genome Atlas (TCGA). Finally, IRNAA also integrates other tools such as `gProfiler` and `WebGestalt` for Gene Pathway Analysis. So far, IRNAA preprocesses both **paired-end** and **single-end** fastq files. IRNAA is a work in progress and might, in the future, integrate other tools such as `STAR`, `samR`, etc.


# Installation
IRNAA is currently built to optimally run on Ubuntu (any version). It may run on other Operating Systems, especially Windows systems with the exception of the user to import read counts from fastq files.

Before installation, please, make sure the following tools are installed and can be easily "called" from the terminal:

- <a href="https://combine-lab.github.io/salmon/" target="_blank">salmon</a> is a command line tool for quantifying the expression of transcripts using RNA-seq data.
- <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">fastQC</a> is a quality control tool for high throughput sequence data built in Java.
- <a href="https://multiqc.info/" target="_blank">multiQC</a> is a command line tool to Aggregate results from bioinformatics analyses across many samples into a single report. Here, `multiQC` is used to aggregate `fastQC` reports across many samples into a single report.
- <a href="https://www.r-project.org/" target="_blank">R >= 3.5.1</a> is a free software environment for statistical computing and graphics.

You may directly download IRNAA from <a href="http://www.github.com/rosericazondekon/irnaa" target="_blank">IRNAA's Github repository</a>. If you have git installed, you may also clone IRNAA by running:

```shell
git clone https://github.com/rosericazondekon/irnaa.git
```

Next, extract IRNAA to your local drive. 

To install it on Ubuntu, from the command line interface, change your directory (using the `cd` command) to the extracted local `irnaa` folder on your local drive, and run the following script:

```shell
sudo bash linux_install.sh && source ~/.bashrc
```

On a Windows machine, double-click on the `win_install.bat` file (**as administrator**) to install IRNAA on your computer. To install it from the DOS command line interface, open your command prompt, change your directory (using the `cd` command) to the extracted local `irnaa` folder and execute the following code:

```batch
win_install
```

You should now launch the IRNAA shiny application by executing the `irnaa` command from either your shell terminal or your DOS command prompt:

```shell
irnaa
```

On windows, you may also run IRNAA by double-clicking on `irnaa.bat` located the `win` directory inside your locally extracted `irnaa` folder.


### Additional requirements for Windows Users
For Windows users, **the installation of IRNAA might require manual installation of certain R libraries**. Also, the installation of <a href="http://gnuwin32.sourceforge.net/downlinks/tar-bin.php" target="_blank">TAR</a> is required to import TCGA data.


# Read Counts from RNA-seq fastq files

### RNA-seq files format and file organization

IRNAA expects RNA-seq files to be in the `.fastq.gz` format. The following file organization is also expected from the end user (with 4 paired-end RNA-seq samples):

```
(working directory)
|
+-- data
|   +-- sample1
|   |   +-- sample1_1.fastq.gz
|   |   +-- sample1_2.fastq.gz
|   +-- sample2
|   |   +-- sample2_1.fastq.gz
|   |   +-- sample2_2.fastq.gz
|   +-- sample3
|   |   +-- sample3_1.fastq.gz
|   |   +-- sample3_2.fastq.gz
|   +-- sample4
|   |   +-- sample4_1.fastq.gz
|   |   +-- sample4_2.fastq.gz
+-- (other folders and files)
```

Similarly, with 4 single-end RNA-seq FASTQ files, the expected file organization is the following:

```
(working directory)
|
+-- data
|   +-- sample1
|   |   +-- sample1.fastq.gz
|   +-- sample2
|   |   +-- sample2.fastq.gz
|   +-- sample3
|   |   +-- sample3.fastq.gz
|   +-- sample4
|   |   +-- sample4.fastq.gz
+-- (other folders and files)
```

According to the above file organization, all RNA-seq samples are organized in the 'data' folder (which can obviously be named otherwise).

To import read counts from RNA-seq fastq files, use the IRNAA shinyapp panel to set-up your working directory. The working directory is the parent directory to the 'data' folder (as portrayed above in the file organization diagram). Once the working directory is set, IRNAA automatically detect the samples and allow the user to run Quality Check (QC) on the fastq files (QC Check menu), and visualize the QC report generated by `multiQC`.


### Transcriptome reference file and index building

IRNAA expects a reference transcriptome in the `.fa` or `.fna` format. Although IRNAA functions nominally with either RefSeq NCBI, or Ensembl genome reference transcripts (in FASTA format), we recommend the use of **RefSeq transcriptome reference FASTA file** available on the <a href="https://www.ncbi.nlm.nih.gov/genome" target="_blank">NCBI website</a>.



Now, enjoy the IRNAA and should you have any suggestion, critique, and/or recommendation, feel free to contact me at <roseric_2000@yahoo.fr>.

Thank You!

Roseric Azondekon

<br/>
Last revision: April 16, 2019
