FROM ubuntu:18.04
ENV TIMEZONE America/Los_Angeles

RUN apt-get update && \
    apt-get install -y \
    git build-essential curl libz-dev libbz2-dev libcurl4-gnutls-dev libgsl0-dev liblzma-dev libncurses5-dev libperl-dev libssl-dev zlib1g-dev autoconf automake wget gawk parallel python-dev python-pip

# Star aligner 
ENV star_version 2.5.4b
ADD https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz /usr/bin/
RUN tar -xzf /usr/bin/${star_version}.tar.gz -C /usr/bin/ && \
    cp /usr/bin/STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin

# Other tools
RUN cd /tmp && git clone -b v1.3 https://github.com/lh3/seqtk.git && \
    cd seqtk && make && make install && \
    cd /tmp && git clone -b 1.9 https://github.com/samtools/htslib && \
    cd htslib && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && git clone -b 1.9 https://github.com/samtools/samtools && \
    cd samtools && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && git clone -b v0.7.17 https://github.com/lh3/bwa.git && \
    cd bwa && make && cp bwa /usr/local/bin && \
    cd /tmp && wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && tar -zxvf bedtools-2.29.1.tar.gz && cd bedtools2 && make && cp ./bin/* /usr/local/bin && \
    cd /tmp && rm -rf * 

# Python
RUN apt-get install -y python3-dev libopenblas-dev python3-pip && \
    pip3 install numpy cython scipy pandas click

RUN apt-get --purge remove -y build-essential git autoconf automake wget

RUN apt-get -y install locales && \
    sed -i -e 's/# \(en_US\.UTF-8 .*\)/\1/' /etc/locale.gen && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV PYTHONUNBUFFERED=0

# Install miniconda to /miniconda
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda create --name hifi_conda 
RUN conda install -n hifi_conda -c bioconda fastp
RUN conda install -n hifi_conda -c bioconda pear

#
COPY ./hifi_* /usr/local/bin/
COPY ./bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/* && mkdir /hifi
WORKDIR /hifi