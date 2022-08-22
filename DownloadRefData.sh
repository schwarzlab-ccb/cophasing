#!/usr/bin/env bash

mkdir -p refs
cd refs
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
echo "unzip"
gunzip hg38.fa.gz
echo "rezip"
bgzip hg38.fa