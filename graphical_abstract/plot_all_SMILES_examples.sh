#!/bin/bash

cd example_SMILES

for s in *.txt
do
    python ../draw_from_example_SMILES_1.py $s $(echo $s | cut -d'.' -f1).tex
done