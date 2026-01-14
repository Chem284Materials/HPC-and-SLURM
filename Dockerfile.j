FROM rancavil/slurm-jupyter:19.05.5-1

USER root

RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    cmake \
    mpich \
    libmpich-dev \
    libblas-dev \
    liblapack-dev \
    python3 \
    python3-pip \
    && apt-get clean

USER admin