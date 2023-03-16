#!/usr/bin/env bash

ARG1=${1:-none}

if [ $ARG1 == "UMAP" ]; then
    if [ ! -d "data/" ]; then
        mkdir data/
    fi

    if [ ! -d "data/raw/" ]; then
        mkdir data/raw/
    fi

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520356&format=file&file=GSM5520356%5Fs9%5Frosa26WT%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*/*.gz

    rm data/raw/dl.tar.gz
fi

python src/vis/UMAP.py