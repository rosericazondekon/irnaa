Integrated RNA-seq Analysis Application (IRNAA)
================
Roseric AZONDEKON

February 14, 2019

<br/>

IRNAA is a shiny application built for the purpose of facilitating Differential Gene Expression (DGE) analysis. IRNAA integrates command line tools such as `salmon`, `fastQC`, `multiQC` with DGE `R` packages such as `DESeq2`, `edgeR`, and `limma-voom`, and hides all the abstractions and coding hurdles from the end user. IRNAA makes it possible for the user to either preprocess fastq files and import read counts into the application from the fastq preprocessed files, or directly import their read counts into the application, or import RNA-seq datasets directly from the The Cancer Genome Atlas (TCGA). Finally, IRNAA also integrates other tools such as `gProfiler` and `WebGestalt` for Gene Pathway Analysis. So far, IRNAA preprocesses both **paired-end** and **single-end** fastq files. IRNAA is a work in progress and might, in the future, integrate other tools such as `STAR`, `samR`, etc.

<br/>
<br/>
Installation
-----------
IRNAA is currently built to optimally run on Ubuntu (any version). It also runs on Windows systems with no to little support for fastq files pre-processing and data import from The Cancer Genome Atlas (TCGA).

You may directly download IRNAA from <a href="http://www.github.com/rosericazondekon/irnaa" target="_blank">IRNAA's Github repository</a>. In case you have `git` installed, you may also clone IRNAA by executing:

```shell
git clone https://github.com/rosericazondekon/irnaa.git
```

### GNU/Linux Ubuntu distribution
Before installation, please, make sure the following tools are installed and can be easily "called" from the terminal:

- <a href="https://combine-lab.github.io/salmon/" target="_blank">salmon</a> is a command line tool for quantifying the expression of transcripts using RNA-seq data.
- <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">fastQC</a> is a quality control tool for high throughput sequence data built in Java.
- <a href="https://multiqc.info/" target="_blank">multiQC</a> is a command line tool to Aggregate results from bioinformatics analyses across many samples into a single report. Here, `multiQC` is used to aggregate `fastQC` reports across many samples into a single report.
- <a href="https://www.r-project.org/" target="_blank">R >= 3.5.1</a> is a free software environment for statistical computing and graphics.

Next, extract IRNAA to your local drive. 

From the command line interface, change your directory (using the `cd` command) to the extracted local `irnaa` folder on your local drive, and run the following command:

```shell
sudo bash linux_install.sh && source ~/.bashrc
```

To launch `IRNAA`, execute the `irnaa` command from your shell terminal:

```shell
irnaa
```

### Windows Systems
On windows computers, before installing IRNAA, make sure that `R` is installed and can be run from the command line.

To make `R` accessible from the command line interface, you may use the following instructions: 
- Search for the `Rscript.exe` executable and copy its location path (e.g `C:\Program Files\R\R-3.5.1\bin`)
- Open the start menu and type in "View advanced system settings", click on "Environment variables"
- Under "System variables", select Path and click on edit
    - **Windows 10 and related** Click "New", and paste in the directory path to `Rscript.exe`.
    - **For earlier versions of Windows systems**, paste in the directory path at the beginning of the path line followed by a semicolon (no space before or after the semicolon)
- Open the Windows Command Prompt and run `Rscript --version` which should output the version of `R` installed on your computer.

Now, from the extracted local `irnaa` folder, right-click on the `win_install.bat` file and run it **as administrator**. When the installation is succesful, a shortcut to `IRNAA` is created on your desktop.

To launch `IRNAA`, double-click on the `IRNAA` shortcut on your desktop.

To install it from the DOS command line interface, open your command prompt (**as administrator**), change your directory (using the `cd` command) to the extracted local `irnaa` folder and execute the following code:

```batch
win_install
```

You may also run IRNAA by double-clicking on the `irnaa.bat` file located the `win` directory inside your locally extracted `irnaa` folder.

### Caution
- For Windows users, **the installation of IRNAA might require manual installation of certain R libraries**. 
- Currently, `IRNAA` is compatible with all major versions of `R 3.4` and `R 3.5`. It is yet to be compatible with the newly released `R 3.6`.

<br/>
<br/>
Read Counts from RNA-seq fastq files
-----------

<br/>

### RNA-seq files format and file organization

IRNAA expects RNA-seq files to be in the `.fastq.gz` format. `IRNAA` expects the **working directory** to be the **parent directory of the folder** containing the samples FASTQ files. The following file organization is expected from the end user (e.g. with 4 paired-end RNA-seq samples):

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

To import read counts from RNA-seq fastq files, use the IRNAA shinyapp panel to set-up your working directory. According to the file organization portrayed above, **the working directory is the parent directory to the 'data' folder**.

Once the working directory is set, IRNAA automatically detects the samples and allows the user to run Quality Check (QC) on the fastq files, and to visualize the QC report generated by `multiQC`.

<br/>

### Transcriptome reference file and index building

IRNAA expects a reference transcriptome in the `.fa` or `.fna` format. Although IRNAA functions nominally with either RefSeq NCBI, or Ensembl transcriptome references (in FASTA format), we recommend the use of **RefSeq transcriptome reference FASTA file** available on the <a href="https://www.ncbi.nlm.nih.gov/genome" target="_blank">NCBI website</a>.

### Read Quantification
For read quantification from FASTQ files, it is important to choose the correct reference source when importing Gene level estimates. By default, `IRNAA` provides gene annotation for Human, Mouse, and Zebrafish. When working with a different species, the species annotation database (e.g. `org.Bt.eg.db` for Bovine) can be specified and `IRNAA` will automatically install it and use it to import gene level estimates. A full list of the Genome wide annotation database packages is available on the <a href="https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData" target="_blank">Bioconductor AnnotationData page</a>.

<br/>
Slider Input
-----
You may use your mouse and the `Left/Right` directional keys on your keyboard to fine tune the numbers of the sliders in the `IRNAA` shinyApp.

<br/>
Stopping IRNAA
-----
To stop `IRNAA`, close your Shell Terminal or Windows Command Prompt.

<br/>
Now, enjoy `IRNAA` and should you have any suggestion, critique, and/or recommendation, feel free to reach out to <roseric_2000@yahoo.fr>.

Thank You!

Roseric Azondekon

<br/>
<hr>
Last revision: May 21, 2019
