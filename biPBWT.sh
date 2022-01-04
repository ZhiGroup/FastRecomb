#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: Bidirectional PBWT (Positional Burrows-Wheeler Transform)"
	echo ""
	echo "Contact: Shaojie Zhang [shzhang@cs.ucf.edu] or Degui Zhi [degui.zhi@uth.tmc.edu]"
	echo ""
	echo "Usage: ./biPBWT.sh [options] parameters"
	echo ""
	echo "Required Parameters:"
	echo -e "\t--readVcf <file>\tVCF file"
	echo ""
	echo "Optional Parameters:"
	echo -e "\t--writeTo <filename>\tOutput filename and location (parameter can be full file path or just filename) [Default = VCF filename]"
	echo -e "\t--length_p <integer>\tError Correcting Block length (in units of sites) [Default = 20]"
	echo -e "\t--width_p <integer>\tError Correcting Block width [Default = 5]"
	echo -e "\t--gap <integer>\t\tGap Size (site) [Default = 1]"
	echo -e "\t--MAF <float>\t\tMinimum Allele Frequency for Error Correction [Default = 0.05]"
	echo -e "\t--checkpoint <integer>\tConsole output every n sites [Default = 100000]"
	exit 1
fi

OPTIONS=c:r:o:g:f:
LONGOPTS=checkpoint:,readVcf:,writeTo:,length_p:,width_p:,gap:,MAF:

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
eval set -- "$PARSED" 

checkpoint=100000 readVcf="" writeTo="" length_p=20 width_p=5 gap=1 MAF=0.05
while true; do
	case "$1" in
		-c|--checkpoint)
			checkpoint="$2"
			shift 2
			;;
		-r|--readVcf)
			readVcf="$2"
			shift 2
			;;
		-o|--writeTo)
			writeTo="$2"
			shift 2
			;;
		--length_p)
			length_p="$2"
			shift 2
			;;
		--width_p)
			width_p="$2"
			shift 2
			;;
		-g|--gap)
			gap="$2"
			shift 2
			;;
		-f|--MAF)
			MAF="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
	esac
done

if [ "$readVcf" == "" ]; then
	echo "The VCF input file must be specified with the required option --readVcf <file>"
	exit 1
fi

basename=$(basename $readVcf)
filename="${basename%.*}"
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running rPBWT..."
./rPBWT "$readVcf" "$writeTo" "$checkpoint"
echo "Running PBWT..."
./PBWT "$readVcf" "$writeTo" "$checkpoint" "$length_p" "$width_p" "$gap" "$MAF"
rm "${writeTo}.rpbwt" "${writeTo}.sites" #"${writeTo}.meta"
echo "biPBWT Finished."

