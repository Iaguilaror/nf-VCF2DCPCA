#!/bin/bash

find -L . \
  -type f \
  -name "*.vcf" \
| sed "s#.vcf#.GTrecoded.txt#" \
| xargs mk
