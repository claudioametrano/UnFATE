ARG VERSION=v202305
# start with miniconda3 as build environment
FROM continuumio/miniconda3 AS build

# Update, install mamba and conda-pack:
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes mamba conda-pack

# Install unfate deps from bioconda

RUN mamba create -c conda-forge -c bioconda -c defaults \
    -n unfate --yes "python>=3.6,<3.7" blast spades exonerate hmmer trimal mafft \
    && conda clean -a -y
RUN mamba install -n unfate --yes biopython==1.76 pandas seaborn click
RUN mamba install -n unfate --yes -c conda-forge parallel
RUN mamba install -n unfate --yes -c etetoolkit ete3 ete_toolchain && conda clean -a -y

# Since we want the most recent, install from repo
SHELL ["conda", "run", "-n", "unfate", "/bin/bash", "-c"]
RUN git clone https://github.com/claudioametrano/UnFATE.git
RUN cd UnFATE

# Install java via apt-get
RUN apt-get update && apt-get install -y default-jre

# add it to the PATH and add env variables
ENV PATH /opt/conda/envs/unfate/bin:$PATH
ENV PATH /UnFATE:$PATH

# When image is run, run the code with the environment
SHELL ["/bin/bash", "-c"]
CMD python3 main_wrap.py