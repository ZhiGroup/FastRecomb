#!/bin/bash
python3.7 res.py "$2" "$1"
python3.7 interpolate.py result."$1".map
python3.7 compRecomb.py result."$1".map.comp
python3.7 regions.py result."$1".map.comp.recomb
python3.7 pearson.py result."$1".map.comp.recomb.color
#Rscript pearson.R result."$1".map.comp.recomb.color
rm result."$1".map result."$1".map.comp result."$1".map.comp.recomb result."$1".map.comp.recomb.color
#mv result."$1".map.comp.recomb.color.png FastRecomb."$1".png
