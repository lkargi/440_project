#!/usr/bin/env bash

ARG1=${1:-none}

if [ ! -d "data/" ]; then
    mkdir data/
fi

if [ ! -d "data/raw/" ]; then
    mkdir data/raw/
fi

if [ $ARG1 == "UMAP" ]; then
    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520356&format=file&file=GSM5520356%5Fs9%5Frosa26WT%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*/*.gz

    rm data/raw/dl.tar.gz

elif [ $ARG1 == "Data" ]; then
    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520356&format=file&file=GSM5520356%5Fs9%5Frosa26WT%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*.gz

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520357&format=file&file=GSM5520357%5Fs9%5Frosa26A%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*/*.gz

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520358&format=file&file=GSM5520358%5Fs9%5Frosa26GA%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*/*.gz

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520359&format=file&file=GSM5520359%5Fs9%5Frosa26GAP%5Fp8%2Etar%2Egz"
    tar -xf data/raw/dl.tar.gz -C data/raw/
    gunzip data/raw/*/*.gz


    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520360&format=file&file=GSM5520360%5FS9WT%2DP15%2Etar%2Egz"
    mkdir data/raw/s9_rosa26WT_p15
    tar -xf data/raw/dl.tar.gz -C data/raw/s9_rosa26WT_p15/
    gunzip data/raw/*/*.gz

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520361&format=file&file=GSM5520361%5FS9Rosa26A%2DP15%2Etar%2Egz"
    mkdir data/raw/s9_rosa26A_p15
    tar -xf data/raw/dl.tar.gz -C data/raw/s9_rosa26A_p15/
    gunzip data/raw/*/*.gz

    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520362&format=file&file=GSM5520362%5FS9Rosa26GA%2DP15%2Etar%2Egz"
    mkdir data/raw/s9_rosa26GA_p15
    tar -xf data/raw/dl.tar.gz -C data/raw/s9_rosa26GA_p15/
    gunzip data/raw/*/*.gz
    
    curl -o data/raw/dl.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5520363&format=file&file=GSM5520363%5FS9Rosa26GAP%2DP15%2Etar%2Egz"
    mkdir data/raw/s9_rosa26GAP_p15
    tar -xf data/raw/dl.tar.gz -C data/raw/s9_rosa26GAP_p15/
    gunzip data/raw/*/*.gz

    rm data/raw/dl.tar.gz

else
    echo "test"
fi



python src/vis/UMAP.py