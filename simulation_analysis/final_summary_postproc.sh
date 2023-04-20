#!/bin/bash

sed -i 's/ \{1,\}/ /g' $1
sed -i 's/& \\phantom{\\_}/ /g' $1
awk 'BEGIN {FS="{"; precol="\\\midrule\\\multicolumn"; precol1="\\\multicolumn"}
    {
        if ($1 == " "precol || $1 == precol || $1 == " "precol1 || $1 == precol1)
            {print $0, "\\\midrule"}
        else
            {print $0}}' $1 > final_$1
