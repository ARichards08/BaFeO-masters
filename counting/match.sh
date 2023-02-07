#!/bin/bash

awk -F '  ' 'NR>138 {print $7}' SrFeO-test.ucf | sort | uniq -c | sort -nr
