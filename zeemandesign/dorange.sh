#!/bin/bash

SERIES=$1

for i in {5..15}; do
  filename="${SERIES}_AWG12Coat_${i}.npz"
  python getlayers.py $filename > ${filename}.log
done
