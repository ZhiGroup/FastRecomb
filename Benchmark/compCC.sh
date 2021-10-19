#!/bin/bash
python3 res.py "$2" "$1"
python3 interpolate.py result."$1".map
python3 compRecomb.py result."$1".map.comp
python3 regions.py result."$1".map.comp.recomb
python3 pearson.py result."$1".map.comp.recomb.color
#Rscript pearson.R result."$1".map.comp.recomb.color
rm result."$1".map result."$1".map.comp result."$1".map.comp.recomb result."$1".map.comp.recomb.color
#mv result."$1".map.comp.recomb.color.png FastRecomb."$1".png
