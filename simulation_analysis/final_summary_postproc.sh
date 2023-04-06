#!/bin/bash

awk 'BEGIN {FS="}"; S="\\midrule\\multicolumn{6"} {if ($1=="  "S || $1==S){print $1"}"$2"}"$3"}"$4"}"$5"}"$6"}"$7"}\\\\ \\midrule"} else {print $0}}' $1 > final_$1
