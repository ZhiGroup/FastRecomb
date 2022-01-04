#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <cassert>

using namespace std;

int M = 0; // # of sequences

// skip meta-info lines and header line
void skipMeta(ifstream& in) {
	string s;
	while (getline(in, s)) {
		if ((int)s.size() < 2 || s[0] != '#' || s[1] != '#') break;
	}
}

// find M using vcf's info field - resets input stream pointer when finished
void getM(ifstream& in) {
	int startPos = in.tellg(); 
	string s; getline(in, s);
	int tabs = 0;
	for (int i = 0; i < (int)s.size(); ++i) {
		if (s[i] == '\t') ++tabs;
		if (tabs >= 9 && s[i] == '|') ++M;
	}
	in.seekg(startPos);
	M *= 2;
}

int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0); cin.tie(0);

	ifstream in(argv[1]);
	ofstream out(string(argv[2]) + ".rpbwt", ios::binary), sites(string(argv[2]) + ".sites"), meta(string(argv[2]) + ".meta");

	int checkpoint = atoi(argv[3]);
	skipMeta(in);
	getM(in);
	
	int site = 0;
	vector<int> positions;
	vector<int> pre(M), div(M); // prefix and divergence array
	iota(pre.begin(), pre.end(), 0);
	vector<int> a(M), b(M), d(M), e(M);
	char s[2 * M + 5000]; // assumes fixed fields take up less than 5000 characters

	in.seekg(-1, ios::end);
	if (in.peek() == '\n') in.seekg(-1, ios::cur); // the process won't work if it starts on a newline character
	in.seekg(-2 * M, ios::cur);
	while (in.tellg() >= 0) {
		// process for reading file backwards
		while (in.peek() != '\n') in.seekg(-1, ios::cur);
		streampos pos = in.tellg();
		in.seekg(1, ios::cur);
		in.getline(s, 2 * M + 5000);
		in.seekg(pos - (streampos)1);
		in.seekg(-2 * M, ios::cur);
		if (s[0] == '#' && s[1] == 'C') break; // header line that starts with "#CHROM POS ID etc."

		// skip fixed fields and get chromosome position
		string position;
		int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
		while (tabs < 9) {
			if (tabs == 1 && s[offset] != '\t') position += s[offset];
			if (s[offset] == '\t') ++tabs;
			++offset;
		}
		positions.push_back(stoi(position));

		// pbwt algorithm
		int u = 0, v = 0, p = site + 1, q = site + 1;
		for (int i = 0; i < M; ++i) {
			int id = pre[i];
			if (div[i] > p) p = div[i];
			if (div[i] > q) q = div[i];
			assert(s[offset + (id / 2) * 4 + 1] == '|'); // sanity check
			if (s[offset + (id / 2) * 4 + (id % 2) * 2] == '0') {
				a[u] = id;
				d[u] = p;
				++u;
				p = 0;
			}
			else {
				b[v] = id;
				e[v] = q;
				++v;
				q = 0;
			}
		}
		for (int i = 0; i < u; ++i) {
			pre[i] = a[i];
			div[i] = d[i];
		}
		for (int i = 0; i < v; ++i) {
			pre[u + i] = b[i];
			div[u + i] = e[i];
		}

		// output
		for (int i = 0; i < M; ++i) {
			out.write((char*)&pre[i], sizeof pre[i]);
			out.write((char*)&div[i], sizeof div[i]);
		}

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
		++site;
	}

	// create file with site positions
	for (int i = (int)positions.size() - 1; i >= 0; --i) {
		sites << positions[i] << '\n';
	}

	// create meta file
	meta << M << ' ' << site << '\n'; // M and N

	in.close();
	out.close();
	sites.close();
	meta.close();

	return 0;
}
