#FROM ubuntu:bionic
FROM tomkellygenetics/cellranger_clean:latest

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y \
 git \
 git-lfs \
 make \
 gzip \
 rename

RUN apt-get install -y \
  curl \
  zsh \
  nano \
  && sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh) --unattended"

RUN git clone "https://github.com/TomKellyGenetics/universc.git"

RUN cd universc/test/cellranger_reference/cellranger-tiny-ref/ \
 && git lfs pull

RUN cd universc \
 && make reference \
 && cd .. 

RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-fastq \
 && ln -s /universc/test/shared/cellranger-tiny-fastq/1.2.0 /cellranger-3.0.2.9001/cellranger-tiny-fastq \
 && ln -s /universc/test/shared/cellranger-tiny-fastq/3.0.0 /cellranger-3.0.2.9001/cellranger-tiny-fastq

RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-ref \
 && ln -s /universc/test/cellranger_reference/cellranger-tiny-ref/1.2.0 /cellranger-3.0.2.9001/cellranger-tiny-ref \ 
 && ln -s /universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0 /cellranger-3.0.2.9001/cellranger-tiny-ref

ENV PATH universc:$PATH

RUN ln -s /universc/launch_universc.sh /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/bin/universc

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get install -y \
 python3-pip \
 && python3 -m pip install --upgrade pip setuptools wheel \
 && pip3 install multiqc

RUN git clone https://github.com/linsalrob/fastq-pair.git \
 && cd fastq-pair \
 && mkdir build \
 && cd build \
 && gcc -std=gnu99   ../main.c ../robstr.c ../fastq_pair.c ../is_gzipped.c  -o fastq_pair \
 && cp fastq_pair /bin/fastq_pair
