#!/bin/bash

sed 's/ \{1,\}/ /g' $1 > final_$1
sed -i 's/& \\phantom{\\_}/ /g' final_$1
