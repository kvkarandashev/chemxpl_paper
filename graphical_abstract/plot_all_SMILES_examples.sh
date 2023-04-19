#!/bin/bash

cd example_SMILES

version=2

for s in *.txt
do
    python ../draw_from_example_SMILES_$version.py $s $(echo $s | cut -d'.' -f1).tex
done