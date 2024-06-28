#include "differentialMaps.hpp"
#pragma GCC optimize("O2")

using namespace std;
using ld = long double;
using ll = long long;

void outputAsIntegers(vector<vector<vector<bool>>> &maps){
    ll n = maps.size();
    for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
        cout << maps[i].size() << endl;
        for (ll columnVectorIndex = 0; (ull) columnVectorIndex < maps[i].size(); columnVectorIndex++){
            ll val = 1;
            ll ret = 0;
            for (ll row = 0; (ull) row < maps[i][columnVectorIndex].size(); row++){
                ret += val * maps[i][columnVectorIndex][row];
                val <<= 1;
            }
            cout << ret << ' ';
        }
        cout << endl << endl;
    }
}

void outputMatrix(vector<vector<vector<bool>>> &maps){
    ll n = maps.size();
    for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
        for (ll columnIndex = 0; (ull)columnIndex < maps[i][0].size(); columnIndex++){
            for (ll rowIndex = 0; (ull)rowIndex < maps[i].size(); rowIndex++){
                cout << maps[i][rowIndex][columnIndex] << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }
}

void outputMatrixGapNotation(vector<vector<vector<bool>>> &maps){
    ll n = maps.size();
    // compare the output on line i+1 with line n+i to get a distance
    
    for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
        cout << '[';
        for (ll columnIndex = 0; (ull)columnIndex < maps[i][0].size(); columnIndex++){
            cout << '[';
            for (ll rowIndex = 0; (ull)rowIndex < maps[i].size(); rowIndex++){
                cout << maps[i][rowIndex][columnIndex];
                if (rowIndex < maps[i].size()-1) cout << ',';
            }
            cout << "]";
            if (columnIndex < maps[i][0].size() - 1) cout << ',';
        }
        cout << ']';
        cout << endl;
    }
}

signed main(){
// takes in the number of crossings followed by all the crossings in
// space separated planar diagram notation in input.txt, and outputs
// all of the differential matrices in output.txt

    int matrixFormat = 0;
    // 0 for 0 and 1 matrices, 1 for gap notation, 2 for integers

    bool reducedHomology = 1;
    // if set to true, it will use reduced homology, with
    // the marked point on strand 1
    bool timeOutput = 0;
    // if set true, then it will simply time how long it takes
    // to compute the differential maps

    bool getMaxSize = 1;
    // if set true, gets the max matrix size

    freopen("input.txt", "r", stdin);
    if (!timeOutput) freopen("output.txt", "w", stdout);

    vector<vector<int>> input;
    int x;
    while (cin >> x){
        if (!input.size() || input[input.size()-1].size() == 4){
            input.push_back({x});
        }
        else{
            input[input.size()-1].push_back(x);
        }
    }

    ld tic = clock();
    PD D = createPlanarDiagram(input);
    ll n = D.size();
    
    vector<vector<vector<bool>>> maps;
    if (reducedHomology) maps = reducedDifferentialMaps(D);
    else maps = regularDifferentialMaps(D);


    if (getMaxSize){
        ll maxSize = 0;
        for (ll i = 0; i < maps.size(); i++){
            maxSize = max(maxSize, (ll)maps[i].size());
        }
        cerr << maxSize << endl;
        return 0;
    }

    if (timeOutput){
        ld tac = clock();
        cerr << (tac - tic) / CLOCKS_PER_SEC * 1000 << " ms" << endl;
        return 0;
    }
    
    if (matrixFormat == 0){
        outputMatrix(maps);
    }
    else if (matrixFormat == 2){
        cout << n << endl;
        outputAsIntegers(maps);
    }
    else{
        outputMatrixGapNotation(maps);
    }

    for (ll i = 0; i < n; i++){
        maps[i] = takeTranspose(maps[i]);
    }

    if (matrixFormat == 0){
        cout << endl << "Transposes Below" << endl;
        outputMatrix(maps);
    }
    else if (matrixFormat == 2){
        outputAsIntegers(maps);
    }
    else{
        outputMatrixGapNotation(maps);
    }
    
    

    // cerr << "done";
    
    
    // for (ll i = 0; i < n; i++){
    //     cout << "Distance from " << i << " 1-resolutions to " << i + 1 << " 1-resolutions is "
    //     << distance(maps[i]) << endl;
    // }

    return 0;
}

// 3 1 4 2 3 2 4 1
// 5 1 8 4 2 6 1 5 3 6 2 7 8 3 7 4
// 1 3 2 4 2 3 1 4