ARG VERSION=v202305
# start with miniconda3 as build environment
FROM continuumio/miniconda3 AS build

# Update, install mamba and conda-pack:
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes mamba conda-pack

# Install unfate deps from bioconda

RUN mamba create -c conda-forge -c bioconda -c defaults -c anaconda \
    -n unfate --yes python=3.11 blast=2.14.1 spades=3.15.5 exonerate=2.4.0 hmmer=3.3.2 trimal=1.4 mafft=7.520 \
    && conda clean -a -y 
RUN mamba install -n unfate --yes click=8.1.7
RUN mamba install -n unfate --yes -c conda-forge parallel=20240722
RUN mamba install -n unfate --yes -c conda-forge openjdk=22.0.1
RUN mamba install -n unfate --yes -c conda-forge biopython=1.80
RUN mamba install -n unfate --yes -c conda-forge pebble=5.0.3 seaborn=0.13.2 matplotlib=3.9.1 progressbar2=4.4.2 scipy=1.14.0 pandas=2.2.2 psutil=6.0.0 pbzip2=1.1.13
RUN mamba install -n unfate --yes -c bioconda  bbmap=39.06 bwa=0.7.18 diamond=2.1.9 samtools=1.19.2
RUN mamba install -n unfate --yes -c conda-forge ete3=3.1.3
RUN mamba install -n unfate --yes -c bioconda hybpiper=2.1.8

# Since we want the most recent, install from repo
SHELL ["conda", "run", "-n", "unfate", "/bin/bash", "-c"]
RUN git clone https://github.com/claudioametrano/UnFATE.git
RUN cd UnFATE

# add it to the PATH and add env variables
ENV PATH /opt/conda/envs/unfate/bin:$PATH
ENV PATH /UnFATE:$PATH

SHELL ["/bin/bash", "-c"]
