#!/bin/bash

awk -F '  ' 'NR>138 {print $7}' SrFeO-test.ucf | sort | uniq -c | sort -nr

awk -F '\t' 'NR>138 {print $7}' BaFeO-original.ucf | sort | uniq -c | sort -nr

awk -v OFS='\t' '$7 ~ /-5.0068125e-21/ {print $2, $3}' BaFeO-original.ucf
