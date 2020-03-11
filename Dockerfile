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

RUN git clone "https://github.com/TomKellyGenetics/cellranger_convert.git"

RUN cd cellranger_convert/test/cellranger_reference/cellranger-tiny-ref/ \
 && git lfs pull

RUN cd cellranger_convert \
 && make reference \
 && cd .. 

RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-fastq \
 && ln -s /cellranger_convert/test/shared/cellranger-tiny-fastq/1.2.0 /cellranger-3.0.2.9001/cellranger-tiny-fastq \
 && ln -s /cellranger_convert/test/shared/cellranger-tiny-fastq/3.0.0 /cellranger-3.0.2.9001/cellranger-tiny-fastq

RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-ref \
 && ln -s /cellranger_convert/test/cellranger_reference/cellranger-tiny-ref/1.2.0 /cellranger-3.0.2.9001/cellranger-tiny-ref \ 
 && ln -s /cellranger_convert/test/cellranger_reference/cellranger-tiny-ref/3.0.0 /cellranger-3.0.2.9001/cellranger-tiny-ref

ENV PATH cellranger_convert:$PATH

RUN ln -s /cellranger_convert/convert.sh /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/bin/conversion
