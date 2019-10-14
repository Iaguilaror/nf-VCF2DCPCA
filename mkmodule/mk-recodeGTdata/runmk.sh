#!/bin/bash

find -L . \
  -type f \
  -name "*.vcf" \
| sed "s#.vcf#.GTrecoded.tsv#" \
| xargs mk
