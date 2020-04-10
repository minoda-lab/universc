![Docker Automated build](https://img.shields.io/docker/automated/tomkellygenetics/universc)
![Docker Build Status](https://img.shields.io/docker/build/tomkellygenetics/universc)
![Docker Image Size (tag)](https://img.shields.io/docker/image-size/tomkellygenetics/universc/latest)
[![GitHub Views](http://hits.dwyl.com/tomkellygenetics/universc.svg)](http://hits.dwyl.com/tomkellygenetics/universc)

![Docker Stars](https://img.shields.io/docker/stars/tomkellygenetics/universc)
![Docker Pulls](https://img.shields.io/docker/pulls/tomkellygenetics/universc)
![Docker Image Version (latest by date)](https://img.shields.io/docker/v/tomkellygenetics/universc)

### UniverSC version 0.3.0

Single-cell processing across technologies.

Provides a conversion script to run multiple technologies and custom libraries with cellranger (10x Genomics analysis tool).


#### Tom Kelly (RIKEN IMS) and Kai Battenberg (RIKEN CSRS/IMS)

## Purpose

We've developed a bash script that will run cellranger on FASTQ files for these technologies. See below for details on how to use it.

If you use this tools, please [cite](#Citation) to acknowledge the efforts of the authors. You can report problems and request
new features to the maintainers with and [issue](#Issues) on GitHub. Details on how to [install](#Install) and [run](#Usage) are provided
below. Please see the [help](#Help) and [examples](#Examples) to try solve your problem before submitting an issue.

### Supported Technologies

In principle, any technology with a cell barcode and unique molecular identifier (UMI) can be supported.

The following technologies have been tested to ensure that they give the expected results: 10x Genomics, Nadia (DropSeq), iCELL8 version 3 

We provide the following preset configurations for convenience based on published data and configurations used by other pipelines 
(e.g, DropSeqPipe and Kallisto/Bustools). To add further support for other technologies or troubleshoot problems, please submit an Issue
to the GitHub repository: https://github.com/TomKellyGenetics/universc/issues] as described in [Bug Reports](#Issues) below.

Some changes to the cellranger install are required to run other technologies. Therefore we provide settings for 10x Genomics
which restores settings for the Chromium instrument. We therefore recommend using 'convert' for processing all data from different
technologies as the tool manages these changes. Please note that multiple technologies cannot be run on the same install of cellranger
at the same time (the tool will also check for this to avoid causing problems with existing runs). Multiple samples of the same technology
can be run simultaneously.

#### Pre-set configurations

-  10x Genomics (version automatically detected): 10x, chromium
--  10x Genomics version 2 (16bp barcode, 10bp UMI): 10x-v2, chromium-v2
--  10x Genomics version 3 (16bp barcode, 12bp UMI): 10x-v3, chromium-v3
-  CEL-Seq (8bp barcode, 4bp UMI): celseq
-  CEL-Seq2 (6bp UMI, 6bp barcode): celseq2
-  Drop-Seq (12bp barcode, 8bp UMI): nadia, dropseq
-  iCell8 version 3 (11bp barcode, 14bp UMI): icell8 or custom
-  inDrops version 1 (19bp barcode, 8bp UMI): indrops-v1, 1cellbio-v1
-  inDrops version 2 (19bp barcode, 8bp UMI): indrops-v2, 1cellbio-v2
-  inDrops version 3 (8bp barcode, 6bp UMI): indrops-v3, 1cellbio-v3
-  Quartz-Seq2 (14bp barcode, 8bp UMI): quartzseq2-384
-  Quartz-Seq2 (15bp barcode, 8bp UMI): quartzseq2-1536
-  Sci-Seq (8bp UMI, 10bp barcode): sciseq
-  SCRB-Seq (6bp barcode, 10bp UMI): scrbseq, mcscrbseq
-  SeqWell (12bp barcode, 8bp UMI): seqwell
-  Smart-seq2-UMI, Smart-seq3 (11bp barcode, 8bp UMI): smartseq
-  SureCell (18bp barcode, 8bp UMI): surecell, ddseq, biorad

#### Custom inputs

Custom inputs are also supported by giving the name "custom" and length of barcode and UMI separated by "_"

 e.g. Custom (16bp barcode, 10bp UMI): custom_16_10

Custom barcode files are also supported for preset technologies. These are particularly useful for well-based
technologies to demutliplex based on the wells.

## Release

At the moment, we have not released the script publicly but we do intend to. We welcome any feedback on it. 
Hopefully it will save people time as make it easier to compare technologies.

We plan to make this open-source with the agreement of everyone in the project.

### Citation <span id="Citation"><span>

### Bug Reports <span id="Issues"><span>

#### Reporting issues

To add further support for other technologies or troubleshoot problems, please submit an Issue 
to the GitHub repository: https://github.com/TomKellyGenetics/universc/issues

### Requesting new technologies

Where possible, please provide an minimal example of the first few lines of each FASTQ file for testing purposes.

It is also helpful to describe the technology, such as:

- length of barcode
- length of UMI
- which reads they're on
- whether there is a known barcode whitelist available
- whether adapters or linkers are required
- whether a preprint, publication, or company specifications are available

Technologies that may be difficult to support are those with:

- barcodes longer than 16bp or varying length 
- combinatorial indexing
- dual indexing 

Please bear this in mind when submitting requests. We will consider to add further technologies but
it could take significant resources to add support for these.

## Installation <span id="Install"><span>

This script requires cellranger to be installed and exported to the PATH (version 3.0.0 of higher recommended).
The script itself is exectuable and does not require installation to run but you can put it in your PATH or 
bin of your cellranger install if you wish to do so.

The conversion process requires write-access to to the cellranger install directory so
an install on your user directory is recommended.

This script will run in bash on any OS (but it has only been tested on Linux Debian). Running cellranger 
with this configuration requires a lot of memory (40Gb) so running on server is recommended.
SGE job modes are supported to run cellranger with multiple threads.

```
make
```


```
sudo make install manual
```

```
make install prefix=$HOME/local
```

-> migrate to Makefile?

### System Requirements

...Unix systems with cellranger installed.

### Docker image

...

## Usage <span id="Usage"><span>

The script will:

- give you a help guide

`bash launch_universc.sh -h`

- convert R1 files so that barcodes and UMIs are where they're expected to be for 10x (this can take some time for larger files)

- runs cellranger with the same parameters as for 10x and treats samples exactly the same

- the barcode whitelists are changed and some checks on barcodes disabled (requires a writeable install of cellranger in your user directory)

- it can run cellranger in parallel in SGE mode on the server if you use `--jobmode "sge"` and set up an `sge.template` file

- it can also restore the original cellranger settings for running 10x samples

`bash launch_universc.sh --setup --technology "10x"`

### Valid barcodes

Please note that this script alters the barcode whitelist. Known iCELL8 barcodes are supported but this is not possible with Nadia or DropSeq chemistry so 100% valid barcodes will be returned.

### Manual <span id="Help"><span>

#### Locally instaleld manual

You can display a manual from the locally installed universc directory with:

```
 man man/launch_universc.sh 
```

Note that the working directory must be `universc` or the full path to the man directory must be given.

#### Installing the manual with root priviliges:

The manual can be installed as follows on Mac and Linux.

```
# add manual directory to PATH if not already found
## check config for Linux
if [[ -f /etc/manpath.config ]]
    then CONFIG="/etc/manpath.config"
fi
## check config for Mac
if [[ -f /etc/manpaths ]]
    then CONFIG="/etc/manpaths"
fi
if [[ ! -z $CONFIG ]]
    then MANDIR=`tail -n 1 ${CONFIG}`
else if [[ ! -z $MANPATH ]]
    then
    SHELL_RC=`echo ~/.${0}rc`
    echo "export MANPATH=/usr/local/man" >> $SHELL_RC
    MANDIR=`echo ${MANPATH} | cut -d: -f1`
    fi
fi
sudo mkdir -p ${MANDIR}/man1
cp man/launch_universc.sh man/launch_universc.sh.1
sudo mv man/launch_universc.sh.1 ${MANDIR}/man1/launch_universc.sh.1
gzip ${MANDIR}/man1/launch_universc.sh.1
```

Alternatively the man can be installed with:

```
cp man/launch_universc.sh man/launch_universc.sh.1
sudo install -g 0 -o 0 -m 0644 man/launch_universc.sh.1 ${MANDIR}/man1
```

The manual can then be called from any directory as follows:

```
man launch_universc.sh
```
`
-> change to make?

#### Help menu

You can access the following help menu with `launch_universc.sh --help` in the terminal.

```
Usage:
  bash launch_universc.sh --testrun -t THECHNOLOGY
  bash launch_universc.sh -t TECHNOLOGY --setup
  bash launch_universc.sh -R1 FILE1 -R2 FILE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -R1 READ1_LANE1 READ1_LANE2 -R2 READ2_LANE1 READ2_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -f SAMPLE_LANE -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -f SAMPLE_LANE1 SAMPLE_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -v
  bash launch_universc.sh -h

Convert sequencing data (FASTQ) from Nadia or iCELL8 platforms for compatibility with 10x Genomics and run cellranger count

Mandatory arguments to long options are mandatory for short options too.
       --testrun                Initiates a test trun with the test dataset
  -R1, --read1 FILE             Read 1 FASTQ file to pass to cellranger (cell barcodes and umi)
  -R2, --read2 FILE             Read 2 FASTQ file to pass to cellranger
  -f,  --file NAME              Path and the name of FASTQ files to pass to cellranger (prefix before R1 or R2)
                                  e.g. /path/to/files/Example_S1_L001

  -i,  --id ID                  A unique run id, used to name output folder
  -d,  --description TEXT       Sample description to embed in output files.
  -r,  --reference DIR          Path of directory containing 10x-compatible reference.
  -t,  --technology PLATFORM    Name of technology used to generate data.
                                Supported technologies:
                                  10x Genomics (version automatically detected): 10x, chromium
                                  10x Genomics version 2 (16bp barcode, 10bp UMI): 10x-v2, chromium-v2
                                  10x Genomics version 3 (16bp barcode, 12bp UMI): 10x-v3, chromium-v3
                                  CEL-Seq (8bp barcode, 4bp UMI): celseq
                                  CEL-Seq2 (6bp UMI, 6bp barcode): celseq2
                                  Drop-Seq (12bp barcode, 8bp UMI): nadia, dropseq
                                  iCell8 version 3 (11bp barcode, 14bp UMI): icell8 or custom
                                  inDrops version 1 (19bp barcode, 8bp UMI): indrops-v1, 1cellbio-v1
                                  inDrops version 2 (19bp barcode, 8bp UMI): indrops-v2, 1cellbio-v2
                                  inDrops version 3 (8bp barcode, 6bp UMI): indrops-v3, 1cellbio-v3
                                  Quartz-Seq2 (14bp barcode, 8bp UMI): quartzseq2-384
                                  Quartz-Seq2 (15bp barcode, 8bp UMI): quartzseq2-1536
                                  Sci-Seq (8bp UMI, 10bp barcode): sciseq
                                  SCRB-Seq (6bp barcode, 10bp UMI): scrbseq, mcscrbseq
                                  SeqWell (12bp barcode, 8bp UMI): seqwell
                                  Smart-seq2-UMI, Smart-seq3 (11bp barcode, 8bp UMI): smartseq
                                  SureCell (18bp barcode, 8bp UMI): surecell, ddseq, biorad
                                Custom inputs are also supported by giving the name "custom" and length of barcode and UMI separated by "_"
                                  e.g. Custom (16bp barcode, 10bp UMI): custom_16_10
  -b,  --barcodefile FILE       Custom barcode list in plain text (with each line containing a barcode)

  -c,  --chemistry CHEM         Assay configuration, autodetection is not possible for converted files: SC3Pv2 (default), SC5P-PE, or SC5P-R2
  -n,  --force-cells NUM        Force pipeline to use this number of cells, bypassing the cell detection algorithm.
  -j,  --jobmode MODE           Job manager to use. Valid options: local (default), sge, lsf, or a .template file
       --localcores NUM         Set max cores the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --localmem NUM           Set max GB the pipeline may request at one time.
                                    Only applies when --jobmode=local.
       --mempercore NUM         Set max GB each job may use at one time.
                                    Only applies in cluster jobmodes.

  -p,  --per-cell-data          Generates a file with basic run statistics along with per-cell data

       --setup                  Set up whitelists for compatibility with new technology and exit
       --as-is                  Skips the FASTQ file conversion if the file already exists

  -h,  --help                   Display this help and exit
  -v,  --version                Output version information and exit
       --verbose                Print additional outputs for debugging

For each fastq file, follow the naming convention below:
  <SampleName>_<SampleNumber>_<LaneNumber>_<ReadNumber>_001.fastq
  e.g. EXAMPLE_S1_L001_R1_001.fastq
       Example_S4_L002_R2_001.fastq.gz

For custom barcode and umi length, follow the format below:
  custom_<barcode>_<UMI>
  e.g. custom_16_10 (which is the same as 10x)

Files will be renamed if they do not follow this format. File extension will be detected automatically.

```

### Examples <span id="Examples"><span>

### Licensing

This package is provided open-source on a GPL-3 license. This means that you are free to use and 
modify this code provided that they also contain this license.

Please note that we are third-party developers releasing it for use by users like ourselves.
We are not affiliated with 10x Genomics, Dolomite Bio, Takara Bio, or any other vendor of
single-cell technologies. This software is not supported by 10x Genomics and only changes
data formats so that other technologies can be used with the cellranger pipeline.

Cellranger (version 2.0.2, 2.1.0, 2.1.0, and 3.0.2) has been released open source on and MIT
license on GitHub. We use this version of cellranger for testing and running our tools.
Note that the code that generates the 'cloupe' files is not included in this release.
The Cloupe browser uses files generated by proprietary closed-source software and is
subject to the 10x Genomics End-User License Agreement which does not allow use with
data generated from other platforms.

Therefore 'launch_universc.sh' does not support Cloupe files and you should not use them with
technologies other than 10x Genomics.  

