FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    wget \
    unzip \
    gzip \
    make \
    zstd \
    gnupg \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libreadline-dev \
    libffi-dev \
    libsqlite3-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libxslt-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libpng-dev \
    libxt-dev \
    libxaw7-dev \
    libglu1-mesa-dev \
    libglfw3-dev \
    libarchive-dev \
    libzip-dev \
    xvfb \
    muscle \
    python3-pip \
    git \
    libssl-dev \
    libcurl4-openssl-dev \
    openjdk-16-jre-headless \
    libnss-sss \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install nextflow and make sure it is configured
RUN wget -qO- https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mv nextflow /usr/local/bin
ENV NXF_HOME=/scratch/.nextflow

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Install mamba
RUN conda install -y -n base -c conda-forge mamba

# Install python libraries
RUN mamba install -y -n base -c anaconda -c conda-forge polars numpy geopy scipy

# Install R channel packages
# (also including UShER because it wants access to BioconductoR)
RUN mamba install -y -n base -c conda-forge -c r -c bioconda r-base \
    r-tidyverse \
    r-geosphere \
    r-viridis \
    r-forcats \
    r-sp \
    r-vcfr \
    r-mass

RUN Rscript -e "install.packages('hzar', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library')"

# install biology-specific CLI packages
RUN mamba install -y -n base -c defaults -c bioconda -c conda-forge \
    csvkit \
    samtools \
    bcftools \
    vcftools \
    bwa \
    nextflow

# make sure shells are bash
CMD ["/bin/bash"]
