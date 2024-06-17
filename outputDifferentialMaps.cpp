#include "differentialMaps.h"
using namespace std;

int main(){
// takes in the number of crossings followed by all the crossings in
// space separated planar diagram notation in input.txt, and outputs
// all of the differential matrices in output.txt

    bool outputMatrixAsIntegers = 0;
    // if set to true, it will output columns as integers
    bool reducedHomology = 0;
    // if set to true, it will use reduced homology, with
    // the marked point on strand 1

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int n;
    cin >> n;
    
    PD D = readPlanarDiagram(n);
    
    vector<vector<vector<bool>>> maps;
    if (reducedHomology) maps = reducedDifferentialMaps(D);
    else maps = differentialMaps(D);

    if (!outputMatrixAsIntegers){
        for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
            for (ll columnIndex = 0; columnIndex < maps[i][0].size(); columnIndex++){
                for (ll rowIndex = 0; rowIndex < maps[i].size(); rowIndex++){
                    cout << maps[i][rowIndex][columnIndex] << ' ';
                }
                cout << endl;
            }
            cout << endl;
        }
    }
    else{
        for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
            for (ll columnVectorIndex = 0; columnVectorIndex < maps[i].size(); columnVectorIndex++){
                ll val = 1;
                ll ret = 0;
                for (ll row = 0; row < maps[i][columnVectorIndex].size(); row++){
                    ret += val * maps[i][columnVectorIndex][row];
                    val <<= 1;
                }
                cout << ret << ' ';
            }
            cout << endl << endl;
        }
    }
    
    
    // for (ll i = 0; i < n; i++){
    //     cout << "Distance from " << i << " 1-resolutions to " << i + 1 << " 1-resolutions is "
    //     << distance(maps[i]) << endl;
    // }

    return 0;
}

// 3 1 4 2 3 2 4 1
// 5 1 8 4 2 6 1 5 3 6 2 7 8 3 7 4
// 1 3 2 4 2 3 1 4