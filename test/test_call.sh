#!/bin/bash

# run tests in universc directory (parent of test directory)
cd $(dirname ${BASH_SOURCE[0]})/..
pwd

## check cellranger installed
cellranger count ---version

git pull --ff-only origin $(git branch --show-current)

## check convert is installed
bash launch_universc.sh -v

