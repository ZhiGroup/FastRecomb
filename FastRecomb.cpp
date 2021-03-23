
#include <iostream>
#include <fstream>
#include <vector>
#include "VCFParser.h"
#include <map>
#include <bitset>
#include <algorithm>
using namespace std;




struct IncGenerator {
	int current_;
	IncGenerator (int start) : current_(start) {}
	int operator() () { return current_++; }
};


/***
 * Create Prefix, Divergence arrays
 */
void build_prefix_divergence_array(vector<unsigned int>& ak,vector<unsigned int>& dk,vector<unsigned char>&alleles, int k, vector<unsigned int> &ak1,vector<unsigned int>&dk1){

	int M = alleles.size();
	//cout << "M is " << M << "\n";
	unsigned int* u = new unsigned int[M];
	unsigned int* v = new unsigned int[M];
	unsigned int* d = new unsigned int[M];
	unsigned int* e = new unsigned int[M];
	unsigned int p = k +1,q = k + 1;
	int a = 0;
	int b = 0;
	dk1[0] = k+1;
	char allele = '0';

	for (int i = 0; i < M ;i++){
		allele = alleles[ak[i]];
		if (dk[i] > p ){
			p = dk[i];
		}if (dk[i] > q){
			q = dk[i];
		}
		if (allele == '0'){
			u[a] = ak[i];
			d[a] = p;
			a++;
			p = 0;
		}else if (allele == '1'){
			v[b] = ak[i];
			e[b] = q;
			b++;
			q = 0;
		}
		else{
			cout << "allele [" <<allele << "] not valid!\n";
		}
	}
	int counter = 0;
	for (int t = 0; t < a; t++){
		//	cout << u[t] << ",";
		ak1[counter] = u[t];
		counter++;
	}
	for (int t = 0; t < b; t++){
		//	cout << v[t] << ",";
		ak1[counter] = v[t];
		counter++;

	}

	counter = 1;
	for (int t = 0; t < a; t++){
		//	cout << u[t] << ",";
		dk1[counter] = d[t];
		counter++;
	}
	for (int t = 0; t < b; t++){
		//	cout << v[t] << ",";
		dk1[counter] = e[t];
		counter++;

	}
	delete[] u;
	delete[] v;
	delete[] d;
	delete[] e;
	//	cout << "\n";

}

double interpolate(vector<int> & x_values, vector<double>& y_values, double x, bool extrapolate )
{
	int size = x_values.size();

	int i = 0;
	if ( x >= x_values[size - 2] )
	{
		i = size - 2;
	}
	else
	{
		while ( x > x_values[i+1] ) i++;
	}
	double x_left = x_values[i], y_left = y_values[i], x_right = x_values[i+1], y_right = y_values[i+1];
	if ( !extrapolate )
	{
		if ( x < x_left ) y_right = y_left;
		if ( x > x_right ) y_left = y_right;
	}

	double dydx = ( y_right - y_left ) / ( x_right - x_left );

	return y_left + dydx * ( x - x_left );
}


int main(int argc, char* argv[]){

	string input_vcf;
	string output_file;
	double centiMoragen_length = 70;
	int num_iterations = 1;
	double minLength = 0.1;
	int window_size = 1000;
	bool smooth_itr = false;
	bool only_end_point = false;
	if (argc < 2){
		cout << "Usage: ./FastRecomb -i <vcf_file> -o output_pbwt -L <genetic_length_cM_total> -d <min_length> -r <num_iteations>  -e [only endpoint] -s [smooth iteration]\n";
		return 1;
	}

	for (int i = 1 ; i < argc; i++){
		string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")) {
			cout << "Usage: ./FastRecomb <vcf_file>\n";
			return 0;
		} else if ((arg == "-i") || (arg == "--input")) {
			if (i + 1 < argc) {
				input_vcf = argv[++i];
			}
		}else if ((arg == "-o") || (arg == "--ouput")) {
			if (i + 1 < argc) {
				output_file = argv[++i];
			}
		}

		else if ((arg == "-L") || (arg == "--cm")) {
			if (i + 1 < argc) {
				centiMoragen_length = atof(argv[++i]);
			}
		}

		else if ((arg == "-r") || (arg == "--iterations")) {
			if (i + 1 < argc) {
				num_iterations = atoi(argv[++i]);
			}
		}

		else if ((arg == "-d") || (arg == "--min_length")) {
			if (i + 1 < argc) {
				minLength = atof(argv[++i]);
			}
		}


		else if ((arg == "-w") || (arg == "--window")) {
			if (i + 1 < argc) {
				window_size = atoi(argv[++i]);
			}
		}
		else if ((arg == "-s") || (arg == "--smooth_itr")) {
			smooth_itr = true;

		}
		else if ((arg == "-e") || (arg == "--only_endpoint")) {
			only_end_point = true;

		}

	}

	VCFParser vcf_parser;
	ifstream input_vcf_file(input_vcf.c_str());
	string line;
	while (getline(input_vcf_file,line)){

		vcf_parser.set_num_samples(line);

		if (vcf_parser.sample_size > 0){
			break;
		}
	}

	int k = 0;
	int  M = vcf_parser.sample_size*2;
	vector<unsigned int> ak(M,0);
	vector<unsigned int> dk(M+1,1);
	vector<unsigned int> ak1(M,0);
	vector<unsigned int> dk1(M+1,0);
	vector<int> genomic_positions;
	int N;
	int maf = 0;
	string chr_name = "1";
	vector<unsigned char> alleles(M,0);
	IncGenerator g(0);

	std::generate( ak.begin(), ak.end(), g);
	//cout << "ak ," << ak[2] << "\n";
	ofstream out(output_file + ".pbwt", ios::binary);
	ofstream out_raw(output_file + ".raw", ios::binary);

	for (int j = 0; j < M; j++) {
		out.write((char*)&ak[j], sizeof ak[j]);
		//if (alleles[j] == '1')
		//	bitVals[j] = 1;
		//else
		//	bitVals[j] = 0;
	}


	for (int j = 0; j < M+1; j++) {
		out.write((char*)&dk[j], sizeof dk[j]);
	}


	vector<char> bitVals(M);
	while (getline(input_vcf_file,line)){


		int pos = vcf_parser.parse_line(line,maf,alleles,chr_name);
		if (pos < 0){
			continue;
		}
		genomic_positions.push_back(pos);

		for (int j = 0; j < M; j++) {
			if (alleles[j] == '1')
				bitVals[j] = '1';
			else
				bitVals[j] = '0';
		}

		//cout << "k:" << k << "\n";
		//cout << "MAF:" << maf << "\n";
		k++;
		build_prefix_divergence_array(ak,dk,alleles,k,ak1,dk1);

		ak = ak1;
		dk = dk1;
		// output
		for (int i = 0; i < M; ++i) {
			//cout << ak[i] << ",";
			out.write((char*)&ak[i], sizeof ak[i]);
		}

		for (int i = 0; i < M+1; ++i) {
			out.write((char*)&dk[i], sizeof dk[i]);
		}

		for (int i = 0; i < M; ++i) {

			// unsigned long n = bitVals.to_ulong() ;
			out_raw.write((char*)&bitVals[i], sizeof(bitVals[i]) ) ;
		}

		//	cout << "\n";

	}

	out.close();
	out_raw.close();

	N = genomic_positions.size();
	double* genetic_map = new double[N+1];
	double* tmp_gen_map = new double[N+1];
	vector<int> median_x_vals;


	for (int w = 0; w < genomic_positions[genomic_positions.size()-1] ; w=w+window_size){
		//recomb_rate_window_tmp[w] = 0;
		int median = (w + w + window_size) / 2;
		median_x_vals.push_back(median);
		//median_y_vals.push_back(0);
	}

	vector<double> median_y_vals(median_x_vals.size(),0);

	for (int k = 0 ; k  < N ; k++){
		genetic_map[k] = (double)genomic_positions[k]/double(1000000);

	}

	genetic_map[N] = genetic_map[N-1];



	unsigned int* ak_array = new unsigned int[M];
	unsigned int* dk_array = new unsigned int[M+1];

	double min_length = minLength;
	int i0 = 0;
	int na = 0;
	int nb = 0;
	int total_global_matches = 0;
	char* bitVals_array = new char[M];
	int iteration_counter = 1;


	std::map<int, double> recomb_rate_window_tmp;
	std::map<int, double> recomb_rate_window;
	std::map<int,double> recomb_rate_window_old;

	for (int w = 0; w < genomic_positions[genomic_positions.size()-1] ; w=w+window_size){
		int median = (w + w + window_size) / 2;
		recomb_rate_window_old[median] = 1;
	}



	while (iteration_counter <= num_iterations){

		for (int w = 0; w < genomic_positions[genomic_positions.size()-1] ; w=w+window_size){
			recomb_rate_window_tmp[w/window_size] = 0;
			//cout << w << "\n";
		}


		//cout << "read the pbwt file..\n";
		ifstream input(output_file + ".pbwt", ios::binary);
		std::ifstream input_raw(output_file + ".raw",std::ios::binary);




		total_global_matches =0;

		// Compute total matches!
		for (int i = 0 ; i < N; i++){
			na = 0;
			nb = 0;
			input_raw.read((char*)bitVals_array, sizeof(bitVals[0])*M  ) ;
			input.seekg(i* (sizeof(ak[i])*M + sizeof(dk[i])*(M+1)),ios_base::beg);
			input.read((char*)ak_array, M * sizeof(ak[0]));
			input.read((char*)dk_array, (M+1) * sizeof(dk[0]));
			for (int j = 0; j < M ; j++){
				if (genetic_map[dk_array[j]] > genetic_map[i] - min_length){
					if (na && nb){
						int total_match_here = min(na,nb);
						//	cout << "i:" << i << "," << "m, "<<total_match_here << "\n";
						total_global_matches += total_match_here;
						/*					for ( int ia = i0; ia < i;ia++){
						for ( int ib = ia+1; ib < i ;ib++){

						}
					}*/
					}

					na = 0 ; nb = 0 ; i0 = j ;
				}
				if (bitVals_array[ak_array[j]] == '0')
					na++;
				else
					nb++;

			}
			/*			if (na && nb){ // report matches
				int t_m = 0;

				for (unsigned int ia = i0; ia < M;ia++){
					for (unsigned int ib = ia+1; ib < M ;ib++){

						if (bitVals_array[ak_array[ia] != bitVals_array[ak_array[ib]]]){
							t_m++;
						}


					}
				}
				total_global_matches += t_m;
			}*/

		}

		cout << "total matches: " << total_global_matches << "\n";
		input_raw.clear();
		input_raw.seekg(0, ios::beg);
		input.clear();
		input.seekg(0, ios::beg);
		na = 0;
		nb = 0;
		//input_raw.close();
		//input_raw(output_file + ".raw",std::ios::binary);
		//int total_matches_2 = 0;


		for (int i = 0 ; i < N; i++){

			na = 0;
			nb = 0;
			input_raw.read((char*)bitVals_array, sizeof(bitVals[0])*M  ) ;
			input.seekg(i* (sizeof(ak[i])*M + sizeof(dk[i])*(M+1)),ios_base::beg);
			input.read((char*)ak_array, M * sizeof(ak[0]));
			input.read((char*)dk_array, (M+1) * sizeof(dk[0]));
			vector<int> block_matches_per_site;
			for (int j = 0; j < M ; j++){
				if (genetic_map[dk_array[j]] > genetic_map[i] - min_length){
					if (na && nb){
						int total_match_here = min(na,nb);
						bool na_true = true;
						if (na != total_match_here)
							na_true = false;
						//	cout << "i:" << i << "," << "m2, "<<total_match_here << "\n";

						block_matches_per_site.push_back(total_match_here);

						for ( int ia = i0; ia < j;ia++){
							if (!only_end_point  && na_true && bitVals_array[ak_array[ia]] == '0'){
								double d_p = dk_array[ia];

								int w_pre_pos = genomic_positions[d_p]/window_size;
								recomb_rate_window_tmp[w_pre_pos] += 1;
							}
							else if (!only_end_point  && !na_true && bitVals_array[ak_array[ia]] == '1') {
								double d_p = dk_array[ia];
								int w_pre_pos = genomic_positions[d_p]/window_size;
								recomb_rate_window_tmp[w_pre_pos] += 1;

							}

						}

						//cout << total_match_here << "\n";
						/*					for ( int ia = i0; ia < i;ia++){
						for ( int ib = ia+1; ib < i ;ib++){

						}
					}*/
					}

					na = 0 ; nb = 0 ; i0 = j ;
				}
				if (bitVals_array[ak_array[j]] == '0')
					na++;
				else
					nb++;
			}


			/*			if (na && nb){ // report matches
				int t_m = 0;

				for (unsigned int ia = i0; ia < M;ia++){
					for (unsigned int ib = ia+1; ib < M ;ib++){

						if (bitVals_array[ak_array[ia] != bitVals_array[ak_array[ib]]]){
							t_m++;
						}


					}
				}
				block_matches_per_site.push_back(t_m);
			}*/
			int sum_t = 0;
			for (size_t l = 0; l < block_matches_per_site.size(); l++){
				sum_t += block_matches_per_site[l];
			}

			//cout << sum_t << "\n";
			int w_pos = genomic_positions[i]/window_size;
			recomb_rate_window_tmp[w_pos] = recomb_rate_window_tmp[w_pos] + sum_t;
		}
		//cout << "last site:" << genomic_positions[genomic_positions.size()-1] << "\n";
		int counter = 0;
		if (!only_end_point){
			total_global_matches = total_global_matches * 2;
		}
		for (int w = 0; w < genomic_positions[genomic_positions.size()-1] ; w=w+window_size){
			int median = (w + w + window_size) / 2;

			recomb_rate_window[median] = (recomb_rate_window_tmp[w/window_size]/double(total_global_matches)) * (centiMoragen_length/(double(window_size)/1000000.0));
			//cout << "rate in window: " << median <<"\t" << recomb_rate_window[median] <<"\t" << recomb_rate_window_tmp[w] << "\n";
			double prev_gen_location = 0 ;

			if (counter > 0){
				prev_gen_location = median_y_vals[counter-1];
			}

			if (smooth_itr){

				recomb_rate_window[median] = 0.5*(recomb_rate_window[median] + recomb_rate_window_old[median]);
				recomb_rate_window_old[median] = recomb_rate_window[median];
			}
			double genetic_map_value = (recomb_rate_window[median] * (double)window_size/1000000.0) + prev_gen_location;
			//cout << recomb_rate_window_tmp[w/window_size] << ",";
			median_y_vals[counter] = genetic_map_value;
			counter++;

		}
		//cout << "\n";

		for (size_t s = 0; s < genomic_positions.size(); s++){
			double y = interpolate(median_x_vals,median_y_vals,genomic_positions[s],false);
			genetic_map[s] = y;
		}


		iteration_counter++;
		input.close();
		input_raw.close();
	}

	ofstream out_f(output_file + ".map", ios::binary);

	out_f << "Genomic Position\tRecombination Rate(cM/Mpbs)\tGenetic Location (cM)\n";
	for (unsigned int j = 0; j < median_x_vals.size() ; j ++){
		out_f << median_x_vals[j] << "\t" <<recomb_rate_window[median_x_vals[j]]  <<"\t"<< median_y_vals[j] << "\n";
	}

	out.close();
	//cout << "total matches: " << total_global_matches << "\n";
	//input.close();
	out_f.close();

	delete[] genetic_map;
	delete[] tmp_gen_map;
	delete[] bitVals_array;
}
