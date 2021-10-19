

#ifndef VCFPARSER_H_
#define VCFPARSER_H_
#include <string>
#include <vector>
#include <iostream>
using namespace std;
static const int OFFSET_SAMPLES_INDEX = 9; // stores the starting index of the sample names in vcf line.
static const int GENOMIC_LOCATION_INDEX = 1; // stores index of genomic location in vcf line.
static const int REF_ALLELE_INDEX = 3; //
static const int ALT_ALLELE_INDEX = 4; //
static const int CHR_INDEX = 0; //

static const int FORMAT_INDEX = 8; // stores the format type index in vcf line.

class VCFParser {
public:
	VCFParser();
	virtual ~VCFParser();
	int parse_line(string& line, int &maf, vector<char> &output,string &chr_name); // returns the genomic location if the line contains a valid genotype data, If so the minor allele frequency
	// is computed and the values are stored in the output vector. Otherwise, -1

	void set_num_samples(string& line);
	vector<string> sample_ids;
	int sample_size;
	vector< vector< char> > window;

private:
	string delimiters = "\t";
};




#endif /* VCFPARSER_H_ */
