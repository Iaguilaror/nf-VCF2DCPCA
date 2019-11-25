#!/bin/bash

find -L . \
  -type f \
  -name "*.GTrecoded.txt" \
| sed "s#.GTrecoded.txt#.DCPCA#" \
| xargs mk
