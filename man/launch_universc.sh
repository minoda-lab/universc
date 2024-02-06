na.\" Manpage for UniverSC
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
Provides a conversion script to run multiple technologies and custom libraries with Cell Ranger (10x Genomics analysis tool).
.SH OPTIONS
       --testrun
            Initiates a test trun with the test dataset. The technology and id must be specified.

                e.g., bash launch_universc.sh -t "10x" -i "test-10x" --testrun

  -R1, --read1 FILE
            Read 1 FASTQ file to pass to Cell Ranger (contains the cell barcodes and umi).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                .fastq .fq .fastq.gz .fq.gz

            Compressed files will be opened automatically. Files will be renamed for
            compatibility with Cell Ranger:

                e.g.,  SRR1873277_R1.fastq will be renamed to SRR1873277_S1_L001_R1_001.fastq

            Names for multiple files can be given, for example multiple lanes:

                --read1 Sample_S1_L001_R1_001.fastq Sample_S1_L002_R1_001.fastq

            Apart from inDrops-v2 or inDrops-v3, all technologies expect barcodes in Read 1.

  -R2, --read2 FILE
            Read 2 FASTQ file to pass to Cell Ranger (contains the transcript reads).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                 .fastq .fq .fastq.gz .fq.gz

            Compressed files will be opened automatically. Files will be renamed for
            compatibility with Cell Ranger:

                 e.g.,  SRR1873277_R2.fastq will be renamed to SRR1873277_S1_L001_R2_001.fastq

            Names for multiple files can be given, for example multiple lanes:

                 --read2 Sample_S1_L001_R2_001.fastq Sample_S1_L002_R2_001.fastq

            UniverSC will automatically switch R1 and R2 to pass barcodes in R1 to Cell Ranger for:

                   inDrops-v2, SLiPT-Seq, SLiPT-Seq2, STRT-Seq-2018

  -I1, --index1 FILE
            Index (I1) FASTQ file to pass to Cell Ranger (OPTIONAL). Contains the indexes 
            for each sample. (In the case of Illumina paired-ends these are the i7 indexes).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                  .fastq .fq .fastq.gz .fq.gz

             Compressed files will be opened automatically. Files will be renamed for
             compatibility with Cell Ranger:

                  e.g.,  SRR1873277_I1.fastq will be renamed to SRR1873277_S1_L001_I1_001.fastq

             Names for multiple files can be given, for example multiple lanes:

                 --index1 Sample_S1_L001_I1_001.fastq Sample_S1_L002_I1_001.fastq

             If index files are not given but are contained in the same directory
             as the files for the reads, they will be inferred from them.

            For example is a file Sample_S1_L001_I1_001.fastq is in the same directory,
            it will be passed to Cell Ranger when launch_universc.sh is called:

                bash launch_universc.sh -t "dropseq" -R1 Sample_S1_L001_R1_001.fastq -R2 Sample_S1_L001_R2_001.fastq

            It is still advisable to demultiplex samples with Illumina bcl2fastq or 'cellranger mkfastq'
            before passing them to UniverSC. Index files are passed to Cell Ranger for QC. For example:

                'cellranger mkfastq' --run=/path/to/illumina/bcls --id=sample-name  --sample-sheet=/path/to/SampleSheet.csv --lanes=1,2

            Or

                /usr/local/bin/bcl2fastq -v --runfolder-dir "/path/to/illumina/bcls"  --output-dir "./Data/Intensities/BaseCalls"\
                                            --sample-sheet "/path/to/SampleSheet.csv" --create-fastq-for-index-reads

            Index 1 file is required for the following technologies, in addition to those requiring Index 2.
            UniverSC will attempt to extract them from Read 1 headers if not found:

                   10x-v1, inDrops-v3, STRT-Seq-C1

           Note that 10x v1 chemistry and inDrops-v3 have unusual naming conventions:
               For 10x v1 the I1 file is the sample index and the R3 file containing
               the UMI sequences is required.
               For inDrops-v3 the R3 file is the sample index and the R4 file containing
               part of the barcode sequence is required.
           Do not rename these files. UniverSC will detect the correct file based on filenames.

  -I2, --index2 FILE
            Index (I2) FASTQ file to pass to Cell Ranger (OPTIONAL). Contains the indexes 
            for each sample. (In the case of Illumina paired-ends these are the i5 indexes).
            Please provide the name of FASTQ file in the working directory or the path to it.
            String must match the name of an exiting file. Files can have any of the
            following extensions:

                  .fastq .fq .fastq.gz .fq.gz

             Compressed files will be opened automatically. Files will be renamed for
             compatibility with Cell Ranger:

                  e.g.,  SRR1873277_I2.fastq will be renamed to SRR1873277_S1_L001_I2_001.fastq

             Names for multiple files can be given, for example multiple lanes:

                 --index1 Sample_S1_L001_I2_001.fastq Sample_S1_L002_I2_001.fastq

             If index files are not given but are contained in the same directory
             as the files for the reads, they will be inferred from them.

            For example is a file Sample_S1_L001_I2_001.fastq is in the same directory,
            it will be passed to Cell Ranger when launch_universc.sh is called:

                bash launch_universc.sh -t "dropseq" -R1 Sample_S1_L001_R1_001.fastq -R2 Sample_S1_L001_R2_001.fastq

            This is sufficent to pass files Sample_S1_L001_I1_001.fastq and Sample_S1_L001_I2_001.fastq
            to Cell Ranger if they are in the same directory.

            It is still advisable to demultiplex samples with Illumina bcl2fastq or 'cellranger mkfastq'
            before passing them to UniverSC. Index files are passed to Cell Ranger for QC. For example:

                'cellranger mkfastq' --run=/path/to/illumina/bcls --id=sample-name  --sample-sheet=/path/to/SampleSheet.csv\
                                   --lanes=1,2 --use-bases-mask y26n,I8n,I8n,Y50n

            Or

                /usr/local/bin/bcl2fastq  -v --runfolder-dir "/path/to/illumina/bcls"  --output-dir "./Data/Intensities/BaseCalls"\
                                             --sample-sheet "/path/to/SampleSheet.csv" --create-fastq-for-index-reads\
                                             --use-bases-mask Y26n,I8n,I8n,Y50n  --mask-short-adapter-reads 0\
                                             --minimum-trimmed-read-length 0

            For dual-indexed technologies such as inDrops-v3, Sci-Seq, SmartSeq3 it is advised to use "bcl2fastq"

            Note that dual indexes are not supported by Cell Ranger. Manually demultiplexing as above into separate
            FASTQ files before processing should work as multiple samples are supported. For example, files names as:

                Sample[ABCD]_S[1234]_L00[12]_R[12]_001.fastq

            These can be processed separately and aggregated together to include all cell barcodes.

            Note: processing dual-indexed files is not stable. If behaviour is not as you expect,
            we welcome you to contact us on GitHub to help you out.

            Index 1 and Index 2 files are required for the following technologies
            UniverSC will attempt to extract them from Read 1 headers if not found:

                   SCI-RNA-Seq, SCI-RNA-Seq3, scifi-seq, Smart-Seq2, Smart-Seq3, STRT-Seq-2i

  -f,  --file NAME
            Path and the name of FASTQ files to pass to Cell Ranger (prefix before R1 or R2)

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
            Sample description to embed in output files, passes to Cell Ranger HTML output.

  -r,  --reference DIR
            Path of directory containing 10x-compatible reference.
            See Cell Ranger documentation on how to generate custom "transcriptomes" or
            download human and mouse references from the 10x Genomics website.

            For convenience we provide pre-generated references for the human genome and
            various model species available for download:

            https://genomec.gsc.riken.jp/gerg/UniverSC/Premade_references/


  -t,  --technology PLATFORM
            Name of technology used to generate data.

                Supported technologies:

                                  10x Genomics (version 2 or 3 automatically detected): 10x, chromium
                                  10x Genomics version 1 (14 bp barcode, 10 bp UMI): 10x-v1, chromium-v1
                                  10x Genomics version 2 (16 bp barcode, 10 bp UMI): 10x-v2, chromium-v2
                                  10x Genomics version 3 (16 bp barcode, 12 bp UMI): 10x-v3, chromium-v3
                                  BD Rhapsody (27 bp barcode, 8 bp UMI): bd-rhapsody
                                  BD Rhapsody v2 enhanced beads (27 bp barcode, 8 bp UMI): bd-rhapsody-v2
                                  CEL-Seq (8 bp barcode, 4 bp UMI): celseq
                                  CEL-Seq2 (6 bp UMI, 6 bp barcode): celseq2
                                  C1 Fluidigm (16 bp barcode, No UMI): c1, fluidgm-c1
                                  C1 CAGE (16 bp, No UMI): c1-cage
                                  C1 RamDA-Seq (16 bp, No UMI): c1-ramda-seq
                                  Drop-Seq (12 bp barcode, 8 bp UMI): dropseq
                                  ICELL8 version 2 (11 bp barcode, No UMI): icell8-non-umi, icell8-v2
                                  ICELL8 version 3 (11 bp barcode, 14 bp UMI): icell8 or custom
                                  ICELL8 5′ scRNA with TCR OR kit (10bp barcode, NO bp UMI): icell8-5-prime
                                  ICELL8 full-length scRNA with Smart-Seq (16 bp barcode, No UMI): icell8-full-length
                                  inDrops version 1 (19 bp barcode, 6 bp UMI): indrops-v1, 1cellbio-v1
                                  inDrops version 2 (19 bp barcode, 6 bp UMI): indrops-v2, 1cellbio-v2
                                  inDrops version 3 (16 bp barcode, 6 bp UMI): indrops-v3, 1cellbio-v3
                                  Nadia (12 bp barcode, 8 bp UMI): nadia, dropseq
                                  MARS-Seq (6 bp barcode, 10 bp UMI): marsseq, marsseq-v1
                                  MARS-Seq2 (7 bp barcode, 8 bp UMI): marsseq2, marsseq-v2
                                  Microwell-Seq (18 bp barcode, 6 bp UMI): microwell
                                  PIP-Seq version 0 (24 bp barcode, 8 bp UMI): pip-seq-v0
                                  PIP-Seq version 1 (16 bp barcode, 6 bp UMI): pip-seq-v1
                                  PIP-Seq version 2 (24 bp barcode, 12 bp UMI): pip-seq-v2
                                  PIP-Seq version 3 (28 bp barcode, 12 bp UMI): pip-seq-v3, fluent-bio-v3
                                  PIP-Seq version 4 (28 bp barcode, 12 bp UMI): pip-seq-v4, fluent-bio-v4
                                  QuartzSeq (6 bp index, no UMI): quartz-seq
                                  Quartz-Seq2 (14 bp barcode, 8bp UMI): quartzseq2-384
                                  Quartz-Seq2 (15 bp barcode, 8bp UMI): quartzseq2-1536
                                  RamDA-Seq (6 bp index, no UMI): ramda-seq
                                  SCI-Seq 2-level indexing (30 bp barcode, 8 bp UMI): sciseq2
                                  SCI-Seq 3-level indexing (40 bp barcode, 8 bp UMI): sciseq3
                                  SCIFI-Seq (27 bp barcode, 8 bp UMI
                                  SCRB-Seq (6 bp barcode, 10 bp UMI): scrbseq, mcscrbseq
                                  SeqWell (12 bp barcode, 8 bp UMI): plexwell, seqwell, seqwells3
                                  Smart-seq (16 bp barcode, No UMI): smartseq
                                  Smart-seq2 (16 bp barcode, No UMI): smartseq2
                                  Smart-seq2-UMI, Smart-seq3 (16 bp barcode, 8 bp UMI): smartseq3
                                  SPLiT-Seq (8 bp UMI, 24 bp barcode): splitseq
                                  SPLiT-Seq v2.1 (10 bp UMI, 24 bp barcode): splitseq2
                                  STRT-Seq (6 bp barcode, no UMI): strt-seq
                                  STRT-Seq-C1 (8 bp barode, 5 bp UMI): strt-seq-c1
                                  STRT-Seq-2i (13 bp barcode, 6 bp UMI): strt-seq-2i
                                  STRT-Seq-2018 (8 bp barcode, 8 bp UMI): strt-seq-2019
                                  SureCell (18 bp barcode, 8 bp UMI): surecell, ddseq, biorad
                                  VASA-plate (6 bp UMI, 6 bp barcode): vasa-plate
                                  VASA-drop (6 bp UMI, 16 bp barcode): vasa-drop
                                Custom inputs are also supported by giving the name "custom" and length of barcode and UMI separated by "_"
                                  e.g. Custom (16 bp barcode, 10 bp UMI): custom_16_10

            A barcode whitelist is provided for all beads or wells for the following technologies:

                  10x Genomics, ICELL8, inDrops-v2, inDrops-v3, SCI-Seq (2-level), SCI-Seq3, SmartSeq3, and QuartzSeq2

            Where no known barcodes are available all possible barcodes of the expected length are
            generated and converted if the permutations have not been computed already.

            Linkers are automatically removed from the following technologies:

                  BD Rhapsody, inDrops-v1, Microwell-Seq, SCI-Seq3 Split-Seq, Smart-Seq2, Smart-Seq3, SureCell

            The following technologies default to non-UMI parameters (others can be forced):

                  ICELL8-v2, RamDA-Seq, Quartz-Seq, Smart-Seq, Smart-Seq2

            The following technologies require Index 1 or Index 2 sequences (see above):

                  Fluidigm C1, ICELL8 full-length, inDrops-v3, Quartz-Seq, RamDA-Seq,
                  SCI-RNA-Seq, SCI-RNA-Seq3, scifi-seq, Smart-Seq2, Smart-Seq3, STRT-Seq-2i, STRT-Seq-C1

            For the following full-length technologies the --chemistry argument can be used to refine parameters:

                  SmartSeq, SmartSeq2


  -b,  --barcodefile FILE
            Custom barcode list in plain text (with each line containing a barcode). Please provide
            the name of a text file in the working directory or the path to it.

  -c,  --chemistry CHEM
            Assay configuration, autodetection is not possible for converted files:

                SC3Pv2 (default), SC3Pv3, SC5P-PE, SC5P-R1, or SC5P-R2

            Chemistry can only be automatically detected for 10x Genomics Chromium v2 or v3 as it
            relies on matches to a barcode whitelist. Setting 'SC3Pv1' for 10x version 1 chemistry
            is recommended. For other technologies we do not recommend changing the chemistry input.
            All samples are converted to contain the barcode and UMI in Read1 as used for 'SC3Pv2'.
            'SC3Pv3' is only used for technologies with longer UMI.

            5′ scRNA-Seq ('SC5P-PE') is available only for 10x Genomics, ICELL8, SmartSeq3, and
            STRT-Seq technologies. All other technologies default to 3′ scRNA-Seq parameters.
            Only 10x Genomics and ICELL8 allow choosing which chemistry parameter to use.

           For full-length single-cell technologies (SmartSeq and SmartSeq2) the chemistry input
           configures the following settings:

              auto or SC3Pv2 (default): 5' tag sequences are removed and all reads are mapped in 3' chemistry
              SC5P-PE: 5' tag sequences are replaced by the 10x TSO sequence and all reads are mapped in 5' chemistry (paired-ends)
              SC5P-R1: reads not containing a 5' tag are filtered out and 5' ends are mapped using read 1 only

  -n,  --force-cells NUM
            Force pipeline to use this number of cells, bypassing the cell detection algorithm.

  -j,  --jobmode MODE
           Job manager to use. Valid options: local (default), sge, lsf, or a .template file

           We recommend to use a cluster configuration when submitting jobs in to a job scheduler.
           DO NOT submit jobs in "local" mode to a slurm, SGE, or LSF cluster as Cell Ranger runs
           multiple threads in parallel by default. Performance is significantly improved using a
           cluster mode with a job scheduler.

          See the Cell Ranger documentation on how to set up a cluster mode with a template file.

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
           Generates a file with basic run statistics along with per-cell data (additional output to Cell Ranger).
           Recommended but disabled by default due to additional runtime required to parse BAM files.
           This provides more accurate summary statistics than Cell Ranger (which uses an average across cells
           that are filtered out).

       --non-umi or --read-only
           Force counting reads by adding or replacing UMI with a mock sequence so that each read
           has a unique UMI.

       --setup
           Set up whitelists for compatibility with new technology. Called automatically when a new
           technology is run and no other technology is running. Recommended to run before submitting
           multiple samples of the same technology. Example:

              bash launch_universc.sh -t "dropseq" --setup

       --as-is
           Skips the FASTQ file conversion if the file already exists and run Cell Ranger on pre-converted file.

  -h,  --help
           Prints the usage and a list of the most commonly used commands.

  -v,  --version                Output version information and exit
           Prints the version of 'UniverSC' and the 'Cell Ranger' version that will be called from the PATH.

       --verbose
           Print additional information to standard out for debugging purposes.

       --version
           Prints the version of 'UniverSC' and the 'Cell Ranger' version that will be called from the PATH.

.SH SEE ALSO
cellranger
.SH BUGS
No known bugs.
.SH AUTHOR
S. Thomas Kelly (tom.kelly [at] riken.jp)
Kai Battenberg (kai.batenberg [at] riken.jp)
.SH LICENSE
Released under the GNU General Public License GPLv3
