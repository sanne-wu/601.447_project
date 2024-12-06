#include "minhash.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "usage: ./minhash <list_name> <kmer size> <bottom-k parameter> <output_name (optional)> <-d (optional)>" << endl;
        return 1;
    }
       
    // hash function
    hash<string> myhash;
    string list_name;
    int t;
    int k;
    try {
        list_name = argv[1];
        t = stoi(argv[2]);
        k = stoi(argv[3]);
    } catch(...) {
        cerr << "one or more type errors in input" << endl;
        return 1;
    }
    if (argc == 4) {
        make_similarity_matrix(list_name, t, myhash, k, cout);
        return 0;
    }
    if (argc == 5) {
        bool flag = (string(argv[argc - 1]) == "-d");
        if (flag) {
            make_similarity_matrix(list_name, t, myhash, k, cout, true);
            return 0;
        } else {
            ofstream temp;
            temp.open(argv[4]);
            make_similarity_matrix(list_name, t, myhash, k, temp);
            temp.close();
            return 0;
        }
    }
    if (argc == 6) {
        if (string(argv[argc - 1]) != "-d") {
            cerr << "unknown last argument" << endl;
            return 1;
        }
        ofstream temp;
        temp.open(argv[4]);
        make_similarity_matrix(list_name, t, myhash, k, temp, true);
        temp.close();
        return 0;
    }
    cerr << "number of arguments incorrect" << endl;
    return 1;


    // testing for large length strings
    /*
    vector<string> strings = {"AAACTGGTCAACTCGGGTAAAGCTGCTCAGGGTCGAGGTTCATTACTCCAAAGTAATTGCTGGCTCAGCGGGTTCTACCAAAAGAACACGCATTCGGGACCGTCACTAACTGGAACGCGTTTATAAAATACTATACCTTCGCCGGATTGTTTCTTTGATATGGCATTCCGTAACGTATCACGTACAGCGTATAGGACTGGACAGACCATCTAGCAACGCGCTAAAGTCAAGCGCGGAACCGCACTCAGCGTTCAAAACCGTCCTTACAACGGTGGGCGAACACATTCAGGTCCGGTTGAGTTTGCTTTGATCCAGCTCTGAATTAGTAGCCTTGTAGGGTTGTAATGTGTAGTATAACTTTAATACTGGGCGCCCTGTTACGAGAGGACCAGTCCACATGGTGGTCCCCCCAGTCCATTCGCCACACTAATGGGTGTGCAGCCGTCATTGGTTTGGGAAACTCCAGCGTAACTCTGGCGCTTGTGCGGTACAAGAGGTCCCTCAGTGGCACGTTTTTAGTCTTTCACTTGGGTCTGTCAACGGACTGAGTATCCTCCGGAATCGAATATTTCACCGTTGTCAAGTTTTGAACCTCGAGCGGGAAAATGATCAATGTGTAACCTAACTTGGCTTCAGGGATTGCGAGGCACATTCTTATTTGTCCGGTAAGAGACACAGCGAGAAGGGGTCAGGACGGATCGCCGTGACTATCAAAAGTATCAAGGTAACCTGCCTCTTAGTCATCAACCCGGACCCTGAGTTAGGTTCCACTCCCATGGCAAGGGATGAGGATCTCCAAACTGAGTCTCCATACTTGTTCTACCAACCACTGTGACGGTAGGGACACACGTTGCGTCTATCGTCCGAGGCCACGCGGGCTGAGATGCTACATCGCTCGTTAGCCTCATCTCGTTCTTGCAGGGCCAGGCACGTTACCAGTAGTGCTCTAGGATGCATTCCCGCTCAAAAGTACAACGACTACCGGCATATTATGGAGACC"
    , "AAACTGGTCAACTCGGGTAAAGCTGCTCAGGGTCGAGGTTCATTACTCCAAAGTAATTGCTGGCTCAGCGGGTTCTACCAAAAGAACACGCATTCGGGACCGTCACTAACTGGAACGCGTTTATAAAATACTATACCTTCGCCGGATTGTTTCTTTGATATGGCATTCCGTAACGTATCACGTACAGCGTATAGGACTGGACAGACCATCTAGCAACGCGCTAAAGTCAAGCGCGGAACCGCACTCAGCGTTCAAAACCGTCCTTACAACGGTGGGCGAACACATTCAGGTCCGGTTGAGTTTGCTTTGATCCAGCTCTGAATTAGTAGCCTTGTAGGGTTGTAATGTGTAGTATAACTTTAATACTGGGCGCCCTGTTACGAGAGGACCAGTCCACATGGTGGTCCCCCCAGTCCATTCGCCACACTAATGGGTGTGCAGCCGTCATTGGTTTGGGAAACTCCAGCGTAACTCTGGCGCTTGTGCGGTACAAGAGGTCCCTCAGTGGCACGTTTTTAGTCTTTCACTTGGGTCTGTCAACGGACTGAGTATCCTCCGGAATCGAATATTTCACCGTTGTCAAGTTTTGAACCTCGAGCGGGAAAATGATCAATGTGTAACCTAACTTGGCTTCAGGGATTGCGAGGCACATTCTTATTTGTCCGGTAAGAGACACAGCGAGAAGGGGTCAGGACGGATCGCCGTGACTATCAAAAGTATCAAGGTAACCTGCCTCTTAGTCATCAACCCGGACCCTGAGTTAGGTTCCACTCCCATGGCAAGGGATGAGGATCTCCAAAGACTTATAATTGCTAGACCGATGTCCGGGAAAAAGGTTGGTCTCGTTCGAGGATTCACCCACCGTATCTGGTTTTAACGAACCTATGGTTGAGAAGACAGTTAGGGGCTATTCGGGCCAGCTCTACTTTTTCGGTGTCTGTTACTGACTCCTACTGGGTAGCGGTCTTTCGTATTGCCAAGGGCCATCCACCGATGACGG"
    , "AAACTGGTCAACTCGGGTAAAGCTGCTCAGGGTCGAGGTTCATTACTCCAAAGTAATTGCTGGCTCAGCGGGTTCTACCAAAAGAACCGGATATGGGTGTCGTTTTGTCCGTTCAGGGCGCGGTAATTGAGGATGTGGCCACAGTATGAGATTGTATGCGATTGTAGTGAAGTCGTGTGGAACACTGTAATGCTACTGACAGGGCGTGCTGCGATCGACTTTACGGACTATGGCGAACGGTTGATGCAATGCCAGATGGCGTTAACCATTCGATGACCTGCGCCGTGAGAGCGCTTTTGTTTATTCAGTATTAGCCGTACGGAGGCATTAGCATTGGTGTTGCCATGTCGTCCGTGTTCGTTATGATGGGCGCCCTGTTACGAGAGGACCAGTCCACATGGTGGTCCCCCCAGTCCATTCGCCACACTAATGGGTGTGCAGCCGTCATTGGTTTGGGAAACACATACGTAGTCCGCGTGTTGGAATGATGATTCCGGGATAGGCTGTGGAGCATTCTGCAAGGGTTAGGTAGGTCGTATGTTGTGGTAGCGCATGGCTTCTCATTGTTGTAGACCTGAGGTGACGAAGAGTGTGGAACCGCAGGAGTCGGAAGGGATTACATTTATGCCTGCGTGTGACGGAATCGAATTCTGTAGTACCGTATGCCAGGGTGGTGAGCCTCTCGTAGGACGTTTTAGGTTTGATAGTATGTACATCAGTAGTTAGAGGTCTGGCAGGGAATTCGCGCTGTCAAGTGGGACGGTATCAACGCAGTTCTCATTGTATTCATTTGCATGCAAGTCCGGTAGGACTGCCGGCTGGACCAACGATGTGTCGGGAAGGAGTTTGCAATAGGTAGTCCGTAGTCCTTCACTGGAAGCCGCGGGTGACAGTGCTCTGTTGATCCGATATTCATGGACTATTTTACGTGGAGACGGGGCCTTAGTATCCTACTGGGTAGCGGTCTTTCGTATTGCCAAGGGCCATCCACCGATGACGG"
    }; 

    vector<unordered_set<string>> kmer_sets;
    int n = strings.size();
    vector<vector<double>> similarity(n, vector<double>(n, 0));
    for (int i = 0 ; i < n; i++) {
        kmer_sets.push_back(make_kmer_set(strings[i], k));
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                similarity[i][j] = 1;
            } else {
                double intersect_size = 0;
                for (const auto& kmer : kmer_sets[i]) {
                    if (kmer_sets[j].find(kmer) != kmer_sets[j].end()) {
                        intersect_size++;
                    }
                }
                double temp =  intersect_size/(kmer_sets[i].size() + kmer_sets[j].size() - intersect_size);
                similarity[i][j] = temp;
                similarity[j][i] = temp;
            }
        }
    }
    cout << "actual" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << similarity[i][j] << "\t";
        }
        cout << endl;
    }
    */
}