#!/bin/bash

awk -F '\t' 'NR>138 {print $7}' * | sort | uniq -c | sort -nr 
