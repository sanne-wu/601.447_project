#include <string>
#include <functional>
#include <queue>
#include <utility>
#include <unordered_set>
#include <fstream>

using std::string;
using std::hash;
using std::priority_queue;
using std::pair;
using std::make_pair;
using std::unordered_set;
using std::max;
using std::ifstream;

unordered_set<string> bottom_k_sketch(const unordered_set<std::string>&, hash<std::string>, size_t);

unordered_set<string> bottom_k_sketch_fromfile(const string&, size_t, hash<string>,  size_t);

double estimate_jaccard(const unordered_set<std::string>&, const unordered_set<std::string>&, hash<std::string>, size_t);

unordered_set<string> make_kmer_set(const string&, size_t);