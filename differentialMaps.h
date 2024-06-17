#include <iostream>
#include <set>
#include <vector>
#include <map>

using namespace std;

using ll = long long;

class PD{
    public:
        vector<vector<int>> crossings; // crossing[i].size() should always be 4
        int size(){
            return crossings.size();
        }
        PD(int n = 0){
            crossings.resize(n);
            for (int i = 0; i < n; i++){
                crossings[i].resize(4);
            }
        }
};

void setMerge(set<int>& a, set<int>& b){ // puts b into a
    for (int x : b){
        a.insert(x);
    }
}

template<typename T> set<T> setComplement(set<T> a, set<T> b){ // a setminus b
    for (T x : b){
        a.erase(x);
    }
    return a;
}

set<set<int>> resolutionCircles(PD diagram, int resolution){
    vector<set<int>> circles;
    for (int i = 0; i < diagram.size(); i++){
        vector<pair<int, int>> strandPairings(2);
        if (!(resolution & 1 << i)){
            // pair crossing[0] with crossing[1], crossing[2] with crossing[3]
            strandPairings[0] = {diagram.crossings[i][0], diagram.crossings[i][1]};
            strandPairings[1] = {diagram.crossings[i][2], diagram.crossings[i][3]};
        }
        else{
            // pair crossing[0] with crossing[3], crossing[1] with crossing[2]
            strandPairings[0] = {diagram.crossings[i][0], diagram.crossings[i][3]};
            strandPairings[1] = {diagram.crossings[i][1], diagram.crossings[i][2]};
        }
        for (int j = 0; j < 2; j++){
            bool found = 0;
            vector<set<int>> toMerge; // merge circles that have overlapping strands
            for (int i = 0; i < circles.size(); i++){
                set<int> circle = circles[i];
                if (circle.count(strandPairings[j].first) || circle.count(strandPairings[j].second)){
                    found = 1;
                    toMerge.push_back(circles[i]);
                    circles.erase(circles.begin() + i);
                    i--;
                }
            }
            while (toMerge.size() > 1){
                setMerge(toMerge[0], toMerge.back());
                toMerge.pop_back();
            }
            if (toMerge.size()){
                toMerge[0].insert(strandPairings[j].first);
                toMerge[0].insert(strandPairings[j].second);
                circles.push_back(toMerge[0]);
            }
            if (!found){
                set<int> newCircle{strandPairings[j].first, strandPairings[j].second};
                circles.push_back(newCircle);
            }
        }
    }
    set<set<int>> ret;
    for (auto s : circles) ret.insert(s);
    return ret;
}

PD readPlanarDiagram(int n){ // reads planar diagram, given n crossings
// space-separated
    PD D(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 4; j++){
            cin >> D.crossings[i][j];
        }
    }
    return D;
}

vector<vector<vector<bool>>> differentialMaps(PD D){
// returns n matrices, each mapping from k 1-resolutions to k+1 1-resolutions for 0 <= k < n
// works with the unreduced Khovanov homology to obtain differentials

    // int n;
    // cout << "Input the number of crossings (no more than 60):" << endl;
    // cin >> n;
    // cout << "Input " << 4 * n << " space-separated numbers in planar diagram notation." << endl;

    // PD D(n);
    // for (int i = 0; i < n; i++){
    //     for (int j = 0; j < 4; j++){
    //         cin >> D.crossings[i][j];
    //     }
    // }

    // testing if resolutionCircles produces correct circles
    // set<set<int>> res = resolutionCircles(D, 2);
    // cout << res.size() << endl;

    // construct resolution cube
    int n = D.size();
    vector<set<set<int>>> resolutionCube(1ll << n);
    for (ll i = 0; i < (1ll << n); i++){
        resolutionCube[i] = resolutionCircles(D, i);
        // cerr << i << ": " << resolutionCube[i].size() << endl;
    }

    vector<ll> ordering((1 << n));
    vector<ll> circleStartingIndex((1 << n));
    // ordering[k] takes in a resolution k (binary string) and gives the
    // zero-indexed order for all elements
    // that have the same number of bits as k

    // circleStartingIndex[k] takes in a resolution k (binary string)
    // and gives the starting index for basis elements of a
    // particular circle compared to all circles with the same number of bits
    vector<ll> bitCount(n+1, 0);
    vector<ll> basisStartCount(n+1, 0);
    for (ll i = 0; i < (1ll << n); i++){
        ordering[i] = bitCount[__builtin_popcountll(i)]++;
        circleStartingIndex[i] = basisStartCount[__builtin_popcountll(i)];
        basisStartCount[__builtin_popcountll(i)] += (1ll << resolutionCube[i].size());
    }
    


    vector<vector<vector<bool>>> differentialMap(n);
    // differentialMap[i] gives the differential map from i 1-resolutions to
    // i+1 1-resolutions
    // differentialMap[i] stores the matrix as a (horizontal) vector of column vectors
    
    // resize differentialMap matrices
    for (int i = 0; i < n; i++){
        differentialMap[i] = vector<vector<bool>>(basisStartCount[i], vector<bool>(basisStartCount[i+1]));
    }

    vector<map<set<int>, ll>> resolutionCircleIndices(1ll << n);
    // get the circle indices of everything in a resolution

    vector<vector<set<int>>> resolutionCirclesVector(1ll << n);
    // store the circles of a resolution in order in a vector

    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        ll count = 0;
        for (auto x : resolutionCube[resolution]){
            resolutionCircleIndices[resolution][x] = count++;
            resolutionCirclesVector[resolution].push_back(x);
        }
    }

    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        for (int j = 0; j < n; j++){
            if ((resolution & (1ll << j)) == 0){ // jth bit not yet set
                ll newResolution = resolution | (1ll << j);
                set<set<int>> oldCircles = resolutionCube[resolution];
                set<set<int>> newCircles = resolutionCube[newResolution];

                set<set<int>> oldDiff = setComplement<set<int>>(oldCircles, newCircles);
                set<set<int>> newDiff = setComplement<set<int>>(newCircles, oldCircles);

                // // get the circle indices of everything in the set of oldCircles
                // map<set<int>, ll> oldCircleIndices;
                // ll count = 0;
                // for (auto x : oldCircles)
                //     oldCircleIndices[x] = count++;

                // // get the circle indices of everything in the set of newCircles
                // map<set<int>, ll> newCircleIndices;
                // count = 0;
                // for (auto x : newCircles)
                //     newCircleIndices[x] = count++;

                // vector<set<int>> oldCirclesVector;
                // for (auto x : oldCircles){
                //     oldCirclesVector.push_back(x);
                // }

                for (ll oldCirclesSubset = 0; oldCirclesSubset < (1ll << oldCircles.size()); oldCirclesSubset++){
                    ll oldIndex = circleStartingIndex[resolution] + oldCirclesSubset;
                    // (-) <-> 0, (+) <-> 1
                    if (oldDiff.size() == 2){ // exactly 2 circles in the old resolution not in the new resolution
                        // must be a merge
                        // (-) x (-) -> (-); (-) x (+) = (+) x (-) -> (+), (+) x (+) -> (-)
                        bool circleOneStatus = ((oldCirclesSubset 
                        & (1ll << resolutionCircleIndices[resolution][*(oldDiff.begin())])) != 0);
                        bool circleTwoStatus = ((oldCirclesSubset 
                        & (1ll << resolutionCircleIndices[resolution][*(++oldDiff.begin())])) != 0);

                        bool newCircleStatus = circleOneStatus ^ circleTwoStatus;
                        // rule based on Audoux's notation

                        ll newCircleIndex = circleStartingIndex[newResolution];
                        // update index for the new merged circle
                        if (newCircleStatus) newCircleIndex += (1ll << resolutionCircleIndices[newResolution][*newDiff.begin()]);
                        for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size(); oldCirclesIndex++){
                            if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                if (newCircles.count(resolutionCirclesVector[resolution][oldCirclesIndex])){
                                    newCircleIndex += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex]]);
                                }
                            }
                        }

                        differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex] = 1;
                    }
                    else if (oldDiff.size() == 1){ // must be a split
                        // (+) -> (+)(+) + (-)(-); (-) -> (-)(+) + (+)(-)
                        // based on Audoux's notation

                        ll newCircleIndex1 = circleStartingIndex[newResolution];
                        ll newCircleIndex2 = circleStartingIndex[newResolution];
                        for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size(); oldCirclesIndex++){
                            if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                if (newCircles.count(resolutionCirclesVector[resolution][oldCirclesIndex])){
                                    newCircleIndex1 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex]]);
                                    newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex]]);
                                }
                            }
                        }

                        // implementing (+) ->
                        if (oldCirclesSubset & (1ll << resolutionCircleIndices[resolution][*oldDiff.begin()])){ // (+) ->
                            // newCircleIndex1: (-)(-), newCircleIndex2: (+)(+)
                            newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][*newDiff.begin()]);
                            newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][*(++newDiff.begin())]);

                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex1] = 1;
                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex2] = 1;
                        }
                        else{ // (-) ->
                            // newCircleIndex1: (+)(-), newCircleIndex2: (-)(+)
                            newCircleIndex1 += (1ll << resolutionCircleIndices[newResolution][*newDiff.begin()]);
                            newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][*(++newDiff.begin())]);

                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex1] = 1;
                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex2] = 1;
                        }
                    }
                }
            }
        }
    }

    // for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
    //     for (ll columnIndex = 0; columnIndex < differentialMap[i][0].size(); columnIndex++){
    //         for (ll rowIndex = 0; rowIndex < differentialMap[i].size(); rowIndex++){
    //             cout << differentialMap[i][rowIndex][columnIndex] << ' ';
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    return differentialMap;
}

// 3 1 4 2 3 2 4 1