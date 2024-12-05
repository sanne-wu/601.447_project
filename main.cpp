#include "minhash.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "missing list filename" << endl;
        return 1;
    }
    string list_name = argv[1];
    
    int k = 31;
    hash<string> myhash;
    int s = 200;

    if (argc == 2) {
        make_similarity_matrix(list_name, k, myhash, s);
        return 0;
    } else if (argc == 3) {
        ofstream temp;
        temp.open(argv[2]);
        make_similarity_matrix(list_name, k, myhash, s, temp);
        temp.close();
        return 0;
    }
    cout << "number of arguments incorrect" << endl;
    return 1;
    
}