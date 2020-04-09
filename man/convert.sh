.\" Manpage for UniverSC
.\" Contact tom.kelly@riken.jp to correct errors or typos.
.TH man 1 "08 April 2020" "0.3" "launch_universc.sh man page"
.SH NAME
launch_universc.sh \- single-cell processing across technologies
.SH SYNOPSIS
 bash launch_universc.sh [--version] [--help] [--setup] [-t <technology>] [-i <id>]
           [-r <reference>] [--option <OPT>]

  bash launch_universc.sh --testrun -t TECHNOLOGY
  bash launch_universc.sh -t TECHNOLOGY --setup
  bash launch_universc.sh -R1 FILE1 -R2 FILE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -R1 READ1_LANE1 READ1_LANE2 -R2 READ2_LANE1 READ2_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -f SAMPLE_LANE -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -f SAMPLE_LANE1 SAMPLE_LANE2 -t TECHNOLOGY -i ID -r REFERENCE [--option OPT]
  bash launch_universc.sh -v
  bash launch_universc.sh -h
.SH DESCRIPTION
Provides a conversion script to run multiple technologies and custom libraries with cellranger (10x Genomics analysis tool).
.SH OPTIONS
       --testrun
            Initiates a test trun with the test dataset. The technology and id must be specified.

                e.g., bash launch_universc.sh -t "10x" -i "test-10x" --testrun

  -R1, --read1 FILE
            Read 1 FASTQ file to pass to cellranger (contains the cell barcodes and umi).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                .fastq .fq .fastq.gz .fq.gz

            Compressed files will be opened automatically. Files will be renamed for
            compatibility with cellranger:

                e.g.,  SRR1873277_R1.fastq will be renamed to SRR1873277_S1_L001_R1_001.fastq

            Names for multiple files can be given, for example multiple lanes:

                --read1 Sample_S1_L001_R1_001.fastq Sample_S1_L002_R1_001.fastq

            Apart from inDrops-v2 or inDrops-v3, all technologies expect barcodes in Read 1.

  -R2, --read2 FILE
            Read 2 FASTQ file to pass to cellranger (contains the transcript reads).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                 .fastq .fq .fastq.gz .fq.gz

            Compressed files will be opened automatically. Files will be renamed for
            compatibility with cellranger:

                 e.g.,  SRR1873277_R2.fastq will be renamed to SRR1873277_S1_L001_R2_001.fastq

            Names for multiple files can be given, for example multiple lanes:

                 --read2 Sample_S1_L001_R2_001.fastq Sample_S1_L002_R2_001.fastq

  -f,  --file NAME
            Path and the name of FASTQ files to pass to cellranger (prefix before R1 or R2)

                e.g. /path/to/files/Example_S1_L001

            This enables giving a prefix instead of "read1" and "read2". This requires
            that there are fastq files ending with the following suffixes:

               ${NAME}_R1_001.fastq and ${NAME}_R2_001.fastq

            Automatic renaming of files and detection of file type may not work in this mode.
            Multiple inputs are still supported for multiple lanes:

                  e.g,. --file Example_S1_L001 Example_S2_L002

                 for files: Example_S1_L001_R1_001.fastq Example_S2_L002_R1_001.fastq
                            Example_S1_L001_R2_001.fastq Example_S2_L002_R2_001.fastq

  -i,  --id ID
            A unique run id, used to name output folder. Must be a string that doesn't
            contain special characters or an existing filename.

  -d,  --description TEXT
            Sample description to embed in output files, passes to cellranger HTML output.

  -r,  --reference DIR
            Path of directory containing 10x-compatible reference.
            See cellranger documentation on how to generate custom "transcriptomes" or
            download human and mouse references from the 10x Genomics website.

  -t,  --technology PLATFORM
            Name of technology used to generate data.

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

           A barcode whitelist is provided for all beads or wells for the following technologies:

                 10x Genomics, iCell8, inDrops-v2, and QuartzSeq2

            Where no known barcodes are available all possible barcodes of the expected length are
            generated and converted if the permutations have not been computed already.

  -b,  --barcodefile FILE
            Custom barcode list in plain text (with each line containing a barcode). Please provide
            the name of a text file in the working directory or the path to it.

  -c,  --chemistry CHEM
            Assay configuration, autodetection is not possible for converted files:

                SC3Pv2 (default), SC3Pv3, SC5P-PE, or SC5P-R2

            Chemistry can only be automatically detected for 10x Genomics Chromium as it relies
            on matches to a barcode whitelist. For other technologies we do not recommend changing
            the chemistry input. All samples are converted to contain the barcode and UMI in Read1
            as used for SC3Pv2. SC3Pv3 is only used for technologies with longer UMI.

  -n,  --force-cells NUM
            Force pipeline to use this number of cells, bypassing the cell detection algorithm.

  -j,  --jobmode MODE
           Job manager to use. Valid options: local (default), sge, lsf, or a .template file

           We recommend to use a cluster configuration when submitting jobs in to a job scheduler.
           DO NOT submit jobs in "local" mode to a slurm, SGE, or LSF cluster as cellranger runs
           multiple threads in parallel by default. Performance is significantly improved using a
           cluster mode with a job scheduler.

          See the cellranger documentation on how to set up a cluster mode with a template file.

       --localcores NUM
           Set max cores the pipeline may request at one time.
           Only applies when --jobmode=local.

       --localmem NUM
           Set max GB the pipeline may request at one time.
           Only applies when --jobmode=local.

       --mempercore NUM
           Set max GB each job may use at one time.
           Only applies in cluster jobmodes.

  -p,  --per-cell-data
           Generates a file with basic run statistics along with per-cell data (additional output to cellranger).
           Recommended but disabled by default due to additional runtime required to parse BAM files.
           This provides more accurate summary statistics than cellranger (which uses an average across cells
           that are filtered out).

       --setup
           Set up whitelists for compatibility with new technology. Called automatically when a new
           technology is run and no other technology is running. Recommended to run before submitting
           multiple samples of the same technology. Example:

              bash launch_universc.sh -t "dropseq" --setup

       --as-is
           Skips the FASTQ file conversion if the file already exists and run cellranger on pre-converted file.

  -h,  --help
           Prints the usage and a list of the most commonly used commands.

  -v,  --version                Output version information and exit
           Prints the version of 'convert' and the 'cellranger' version that will be called from the PATH.

       --verbose
           Print additional information to standard out for debugging purposes.

       --version
           Prints the version of 'convert' and the 'cellranger' version that will be called from the PATH.

.SH SEE ALSO
cellranger
.SH BUGS
No known bugs.
.SH AUTHOR
S. Thomas Kelly (tom.kelly [at] riken.jp)
Kai Battenberg (kai.batenberg [at] riken.jp)
.SH LICENSE
Released under the GNU General Public License GPLv3
