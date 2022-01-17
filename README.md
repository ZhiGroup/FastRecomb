# FastRecomb

## Introduction
FastRecomb is an efficient method to infer recombiantion rates using a biobank-scale haplotype panel. We recommend using P-smoother before running FastRecomb, epseically if the genotyping error rate is not low (e.g. > 0.1%).  The input file format for FastRecomb is phased VCF.

## Dependencies
- C++ (at least GCC 5)  
- GNU Make  

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/FastRecomb.git`

Enter the Debug folder and compile the program:

`cd FastRecomb/Debug`  
`make`


## Usage:

<pre>
Usage: ./FastRecomb -i [vcf_file] -o [output_folder] -L [genetic_length_cM_total] -d [min_length] -r [num_iterations] -w [window_size]

Required parameters:
-i [vcf_file]
     Path to the VCF input file.
-o [output_folder]
        Output file will be written in output_folder including any intermediate file.
-L [genetic_length_cM_total]
        The expected total length of the chromosome in cM.

Optional parameters:

-d [min_length]
        Minimum genetic length of the haplotype matches. The minium length is 0.5 cM by default.

-r [num_iterations]
        Number of iterations. The values is set to 5 by default.
        
-w [window_size]
        Window size in bps. It is set to 5000 by default. A recombination rate is calculated for each window across the chromosome.

 </pre>

