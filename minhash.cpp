#include "minhash.h"

unordered_set<string> bottom_k_sketch(const unordered_set<string>& kmers, hash<string> my_hash,  size_t k) {
    if (k >= kmers.size()) {
        return kmers;
    }
    priority_queue<pair<size_t, string>> pq;
    for (const auto& kmer : kmers) {
        pair<size_t, string> temp = make_pair(my_hash(kmer), kmer);
        if (pq.size() < k) {
            pq.push(temp);
        } else if (temp < pq.top()) {
            pq.pop();
            pq.push(temp);
        }
    }
    unordered_set<string> sketch;
    while (!pq.empty()) {
        sketch.insert(pq.top().second);
        pq.pop();
    }
    return sketch;
}

unordered_set<string> bottom_k_sketch_fromfile(const string& filename,  size_t kmer_size, hash<string> my_hash,  size_t k) {
    ifstream genome_input;
    genome_input.open(filename);

    string kmer;
    char c;
    
    priority_queue<pair<size_t, string>> pq;

    while (genome_input.get(c)) {
        if (!is_nucleotide(c)) continue;
        if (kmer.length() < kmer_size) {
            kmer.push_back(c);
        } else {
            kmer.erase(kmer.begin());
            kmer.push_back(c);
            pair<size_t, string> temp = make_pair(my_hash(kmer), kmer);
            if (pq.size() < k) {
                pq.push(temp);
            } else if (temp < pq.top()) {
                pq.pop();
                pq.push(temp);
            }
        }
    }
    genome_input.close();

    unordered_set<string> sketch;
    while (!pq.empty()) {
        sketch.insert(pq.top().second);
        pq.pop();
    }
    return sketch;
}

double estimate_jaccard(const unordered_set<string>& A, const unordered_set<string>& B, hash<string> my_hash, size_t k) {
    if (k == 0) {
        k = max(A.size(), B.size());
    }
    unordered_set<string> union_set(A);
    for (const auto& kmer : B) {
        union_set.insert(kmer);
    }
    
    unordered_set<string> union_set_sketch = bottom_k_sketch(union_set, my_hash, k);
    double intersect_size = 0;
    for (const auto& kmer : union_set_sketch) {
        if (A.find(kmer) != A.end() && B.find(kmer) != B.end()) {
            intersect_size++;
        }
    }
    return intersect_size/union_set_sketch.size();
}

unordered_set<string> make_kmer_set(const string& s, size_t k) {
    unordered_set<string> kmer_set;
    for (size_t i = 0; i < s.length() - k + 1; i++) {
        kmer_set.insert(s.substr(i, k));
    }
    return kmer_set;
}

void make_similarity_matrix(const string& filename, size_t kmer_size, hash<string> my_hash, size_t k, ostream& out) {
    // out << "Source: " << filename << endl;
    ifstream file_list;
    file_list.open(filename);
    string file;
    string name;
    unordered_map<string, string> name_to_file;
    vector<string> names;
    while (file_list >> file >> name) {
        name_to_file[name] = file;
        names.push_back(name);
    }
    size_t n = names.size();
    vector<vector<double>> similarity_matrix(n, vector<double>(n, 0));
    unordered_map<string, unordered_set<string>> sketches;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            if (i == j) {
                similarity_matrix[i][j] = 1;
            } else {
                string A_name = names[i];
                string B_name = names[j];
                if (sketches.find(A_name) == sketches.end()) {
                    sketches[A_name] = bottom_k_sketch_fromfile(name_to_file[A_name], kmer_size, my_hash, k);
                }
                if (sketches.find(B_name) == sketches.end()) {
                    sketches[B_name] = bottom_k_sketch_fromfile(name_to_file[B_name], kmer_size, my_hash, k);
                }
                double jaccard_estimate = estimate_jaccard(sketches[A_name], sketches[B_name], my_hash, k);
                similarity_matrix[i][j] = jaccard_estimate;
                similarity_matrix[j][i] = jaccard_estimate;
            }
        }
    }
    out << "##NAMES";
    for (size_t i = 0; i < n; i++) {
        out << "\t" << names[i];
    }
    out << endl;
    for (size_t i = 0; i < n; i++) {
        out << names[i];
        for (size_t j = 0; j < n; j++) {
            out << "\t" << similarity_matrix[i][j];
        }
        out << endl;
    }
}

bool is_nucleotide(char c) {
    return (c == 'A' || c == 'T' || c == 'C' || c == 'G');
}