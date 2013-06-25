#!/bin/bash
# $1: ratio of original image file (in percentage)
cd image/original;
for f in *jpg;
do convert -resize $1% $f ../$f;
done
