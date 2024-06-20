#include "differentialMaps.hpp"
#pragma GCC optimize("O2")

using namespace std;
using ld = long double;

int main(){
// takes in the number of crossings followed by all the crossings in
// space separated planar diagram notation in input.txt, and outputs
// all of the differential matrices in output.txt

    bool outputMatrixAsIntegers = 0;
    // if set to true, it will output columns as integers
    bool reducedHomology = 0;
    // if set to true, it will use reduced homology, with
    // the marked point on strand 1
    bool timeOutput = 0;
    // if set true, then it will simply time how long it takes
    // to compute the differential maps

    freopen("input.txt", "r", stdin);
    if (!timeOutput) freopen("output.txt", "w", stdout);

    int n;
    std::cin >> n;

    ld tic = clock();
    PD D = readPlanarDiagram(n);
    
    std::vector<std::vector<std::vector<bool>>> maps;
    if (reducedHomology) maps = reducedDifferentialMaps(D);
    else maps = differentialMaps(D);

    if (timeOutput){
        ld tac = clock();
        cerr << (tac - tic) / CLOCKS_PER_SEC * 1000 << " ms" << endl;
        return 0;
    }
    
    if (!outputMatrixAsIntegers){
        for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
            for (ll columnIndex = 0; (ull) columnIndex < maps[i][0].size(); columnIndex++){
                for (ll rowIndex = 0; (ull) rowIndex < maps[i].size(); rowIndex++){
                    std::cout << maps[i][rowIndex][columnIndex] << ' ';
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    else{
        for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
            for (ll columnVectorIndex = 0; (ull) columnVectorIndex < maps[i].size(); columnVectorIndex++){
                ll val = 1;
                ll ret = 0;
                for (ll row = 0; (ull) row < maps[i][columnVectorIndex].size(); row++){
                    ret += val * maps[i][columnVectorIndex][row];
                    val <<= 1;
                }
                std::cout << ret << ' ';
            }
            std::cout << std::endl << std::endl;
        }
    }
    
    

    // std::cerr << "done";
    
    
    // for (ll i = 0; i < n; i++){
    //     std::cout << "Distance from " << i << " 1-resolutions to " << i + 1 << " 1-resolutions is "
    //     << distance(maps[i]) << std::endl;
    // }

    return 0;
}

// 3 1 4 2 3 2 4 1
// 5 1 8 4 2 6 1 5 3 6 2 7 8 3 7 4
// 1 3 2 4 2 3 1 4