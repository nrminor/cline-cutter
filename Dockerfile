FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set home directory
ENV HOME /opt

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
    perl \
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
    r-mass \
    r-devtools \
    bioconductor-rhdf5

RUN Rscript -e "devtools::install_github('GrahamDB/hzar', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library', clean = TRUE)"
RUN Rscript -e "install.packages('SimDesign', repos = 'https://cloud.r-project.org/', lib='/opt/conda/lib/R/library', clean = TRUE)"

# install biology-specific CLI packages
RUN mamba install -y -n base -c defaults -c bioconda -c conda-forge \
    xlsx2csv \
    csvkit \
    samtools \
    bcftools \
    vcftools \
    bwa \
    popgen-entropy

# pull in entropy source code to get access to the package's perl, python,
# and r scripts, and add them to $PATH
RUN mkdir /opt/entropy && \
    wget -q https://bitbucket.org/buerklelab/mixedploidy-entropy/get/246ccf1003c4.zip && \
    mv 246ccf1003c4.zip /opt/entropy/ && \
    cd /opt/entropy && \
    unzip 246ccf1003c4.zip && \
    cd buerklelab-mixedploidy-entropy-246ccf1003c4 && \
    find . -type f -name "*.pl" -print0 | xargs -0 chmod +x && \
    find . -type f -name "*.sh" -print0 | xargs -0 chmod +x && \
    find . -type f -name "*.py" -print0 | xargs -0 chmod +x && \
    find . -type f -name "*.R" -print0 | xargs -0 chmod +x
ENV PATH $PATH:/opt/entropy/buerklelab-mixedploidy-entropy-246ccf1003c4:/opt/entropy/auxfiles/buerklelab-mixedploidy-entropy-246ccf1003c4:/opt/entropy/buerklelab-mixedploidy-entropy-246ccf1003c4/simfiles/diploid

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-linux-x86_64.tar.gz && \
    tar -xzf julia-1.10.0-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/julia-1.10.0/bin/julia /usr/local/bin/julia && \
    rm julia-1.10.0-linux-x86_64.tar.gz

# Copy Julia dependencies to precompile as a module
ENV JULIA_DEPOT_PATH /opt/.julia
ENV JULIA_SCRATCH_TRACK_ACCESS 0
ENV JULIA_HISTORY /scratch/.julia_history
COPY Project.toml /opt/.julia/environments/v1.10/Project.toml
COPY Manifest.toml /opt/.julia/environments/v1.10/Manifest.toml
RUN julia -e 'using Pkg; \
            Pkg.activate(joinpath(DEPOT_PATH[1], "environments", "v1.10")); \
            Pkg.instantiate(); \
            Pkg.precompile()'
RUN julia -e 'using CSV, DataFrames, Pipe, VariantCallFormat, VCFTools'
ENV PATH=$PATH:/opt/julia-1.10.0/bin:/scratch/.julia/compiled/v1.10

# make sure bin files are executable
RUN chmod +x /usr/local/bin/* && \
    chmod +x /opt/julia-1.10.0/bin/* && \
    rm -f /opt/conda/bin/cpp && \
    chmod +x /opt/conda/bin/* && \
    chmod +rw /opt/ && \
    chmod -R +rwx /opt/ && \
    chmod -R +rwx /opt/.julia/ && \
    chmod -R +rw /opt/.julia/logs/ && \
    chmod -R +rwx /opt/.julia/logs/* && \
    chmod -R +rwx /opt/.julia/compiled/ && \
    chmod -R +rwx /opt/.julia/compiled/v1.10/* && \
    chmod -R +rwx /opt/.julia/logs/* && \
    chmod -R +rwx /opt/.julia/packages/

# make sure shells are bash
CMD ["/bin/bash"]
