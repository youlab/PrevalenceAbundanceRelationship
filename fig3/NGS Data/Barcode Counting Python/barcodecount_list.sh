#!/bin/bash

"""
By: Andrea Weiss
Last edited: 15 April 2020

This code is designed to run the python script (specifically the barcode matching algorithm) 
through all files in the directory specified in IN1 while matching to target barcode sequences in IN2


INPUTS 
IN1 = path to sequencing reads to serach for matching barcodes
In2 = path to location of reference barcodes to match

Inside of loop specify matching algorithm to use

Usage:
./barcodecount_list.sh
NGS/HGT_Trial3/Data/DeMu/Strains/Read2
"""


IN1=$(ls ~/Documents/NGS/KK_16Nov23/1/*.fastqsanger)
IN2=~/Documents/NGS/KK_16Nov23/barcodes.txt

for f in $IN1
do
    ./Naive_withMismatch_Barcode_V2.py $f $IN2 2
done