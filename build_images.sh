#!/bin/sh

docker build -t slurm-node-hpc -f Dockerfile.n .
docker build -t slurm-master-hpc -f Dockerfile.m .
docker build -t slurm-jupyter-hpc -f Dockerfile.j .
