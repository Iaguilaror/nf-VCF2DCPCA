#!/bin/bash

find -L . \
  -type f \
  -name "*.txt" \
| sed "s#.txt#.DCPCA#" \
| xargs mk
