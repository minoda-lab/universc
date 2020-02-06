### UniverSC version 0.2.2

Conversion script to run Nadia and iCELL8 libraries with cellranger (10x Genomics analysis tool)

#### Tom Kelly (RIKEN IMS) and Kai Battenberg (RIKEN CSRS/IMS)

## Purpose

We've developed a bash script that will run cellranger on FASTQ files for these technologies. See below for details on how to use it.

## Release

At the moment, we have not released the script publicly but we do intend to. We welcome any feedback on it. 
Hopefully it will save people time as make it easier to compare technologies.

We plan to make this open-source with the agreement of everyone in the project.

## Installation

This script requires cellranger to be installed and exported to the PATH (version 3.0.0 of higher recommended).
The script itself is exectuable and does not require installation to run but you can put it in your PATH or 
bin of your cellranger install if you wish to do so.

The conversion process requires write-access to to the cellranger install directory so
an install on your user directory is recommended.

This script will run in bash on any OS (but it has only been tested on Linux Debian). Running cellranger 
with this configuration requires a lot of memory (40Gb) so running on server is recommended.
SGE job modes are supported to run cellranger with multiple threads.

## Usage

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

This is a work-in-progress and documentation with examples will be added in the future. The script is stable and functional.
Please send feedback, comments, or issues to Kai Battenberg <[kai.battenberg@riken.jp](mailto:kai.battenberg@riken.jp)>
 or Tom Kelly <[tom.kelly@riken.jp](mailto:tom.kelly@riken.jp)>
