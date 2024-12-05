#include <string>
#include <functional>
#include <queue>
#include <utility>
#include <unordered_set>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>

using std::string;
using std::hash;
using std::priority_queue;
using std::pair;
using std::make_pair;
using std::unordered_set;
using std::max;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::vector;
using std::unordered_map;
using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using std::stoi;

template <typename T>
unordered_set<size_t> bottom_k_sketch(const unordered_set<T>& S, hash<T> my_hash,  size_t k);

unordered_set<size_t> bottom_k_sketch(const unordered_set<size_t>& S, size_t k);

unordered_set<size_t> bottom_k_sketch_fromfile(const string& filename,  size_t kmer_size, hash<string> my_hash,  size_t k);

double estimate_jaccard(const unordered_set<size_t>& A, const unordered_set<size_t>& B, size_t k=0);

unordered_set<string> make_kmer_set(const string& s, size_t k);

void make_similarity_matrix(const string& filename, size_t kmer_size, hash<string> my_hash, size_t k, ostream& out=cout);

bool is_nucleotide(char c);