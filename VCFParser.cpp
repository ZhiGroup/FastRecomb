

#include "VCFParser.h"

VCFParser::VCFParser() {
	sample_size = 0;
}

VCFParser::~VCFParser() {
}



void VCFParser::set_num_samples(string &line){


	if (line.size() < 1 || (this->sample_size > 0 && line[0] == '#')){
		return;
	}

	string::size_type lastPos = line.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = line.find_first_of(delimiters, lastPos);
	int index_counter = 0;
	bool sample_line = false;

	if (string::npos != pos || string::npos != lastPos)
	{
		if (this->sample_size == 0){
			string header =line.substr(lastPos, pos - lastPos);
			if ("#CHROM" == header){
				sample_line = true;
			}

			if (sample_line == false and '#' == line[0]){
				return;
			}
			if (sample_line){
				while (string::npos != pos || string::npos != lastPos){

					if (index_counter >= OFFSET_SAMPLES_INDEX){
						this->sample_ids.push_back(line.substr(lastPos, pos - lastPos));
						this->sample_size++;
					}
					lastPos = line.find_first_not_of(delimiters,pos);
					pos = line.find_first_of(delimiters, lastPos);
					index_counter++;
				}
			}

		//	cout << this->sample_ids.size() << "\n";
		}
	}
}



int VCFParser::parse_line(string &line,int &maf, vector<unsigned char>& output,string &chr_name){

	string::size_type tmp_pos;

	int genomic_location = -1;
	int num_zeros = 0;
	int num_ones = 0;
	int hap_num = sample_size*2;
	if (line.size() < 1 || (this->sample_size > 0 && line[0] == '#')){
		return genomic_location;
	}

	string::size_type lastPos = line.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = line.find_first_of(delimiters, lastPos);
	int index_counter = 0;

	int idx = 0;
	if (string::npos != pos || string::npos != lastPos)
	{

		while (string::npos != pos || string::npos != lastPos){

			if (index_counter >= OFFSET_SAMPLES_INDEX){
				string _values = line.substr(lastPos, pos - lastPos);
				tmp_pos = _values.find_first_of("|");
				unsigned char first_hap =_values.substr(0,tmp_pos)[0];
				unsigned char second_hap =_values.substr(tmp_pos+1)[0];

				if (first_hap == '1') num_ones++;

				else if (first_hap == '0') num_zeros++;
				else
					return -1;

				if (second_hap == '1') num_ones++;

				else if (second_hap== '0') num_zeros++;
				else
					return -1;

				if (idx >= hap_num)
					return -1;
				output[idx++] =  first_hap ;
				output[idx++] =  second_hap;
			}
			else if (index_counter == GENOMIC_LOCATION_INDEX){
				genomic_location = atoi(line.substr(lastPos, pos - lastPos).c_str());
			}

			else if (index_counter == FORMAT_INDEX){
				string format = line.substr(lastPos, pos - lastPos);
				if (format != "GT"){
					return -1;
				}

			}

/*			else if (index_counter == REF_ALLELE_INDEX){
				string ref_allele = line.substr(lastPos, pos - lastPos);
				if (ref_allele.size() > 1){
					return -1;
				}

			}

			else if (index_counter == ALT_ALLELE_INDEX){
				string alt_allele = line.substr(lastPos, pos - lastPos);
				if (alt_allele.size() > 1){
					return -1;
				}

			}*/
			else if (index_counter == CHR_INDEX){
				chr_name = line.substr(lastPos, pos - lastPos);


			}

			lastPos = line.find_first_not_of(delimiters,pos);
			pos = line.find_first_of(delimiters, lastPos);
			index_counter++;
		}
	}

	if (genomic_location >= 0){
		if (num_ones > num_zeros){
			maf = num_zeros;
		} else{
			maf = num_ones;
		}

		return genomic_location;
	}

	return -1;
}



