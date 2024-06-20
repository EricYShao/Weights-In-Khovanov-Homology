#include <iostream>
#include <set>
#include <vector>
#include <map>
#ifndef DIFFERENTIAL_MAPS

using ll = long long;
using ull = unsigned long long;

class PD{
    public:
        std::vector<std::vector<int>> crossings; // crossing[i].size() should always be 4
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

void setMerge(std::set<int>& a, std::set<int>& b){ // puts b into a
    for (int x : b){
        a.insert(x);
    }
}

template<typename T> std::set<T> setComplement(std::set<T> a, std::set<T> b){ // a setminus b
    for (T x : b){
        a.erase(x);
    }
    return a;
}

std::set<std::set<int>> resolutionCircles(PD diagram, int resolution){
    std::vector<std::set<int>> circles;
    for (int i = 0; i < diagram.size(); i++){
        std::vector<std::pair<int, int>> strandPairings(2);
        if (!(resolution & (1ll << i))){
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
            std::vector<std::set<int>> toMerge; // merge circles that have overlapping strands
            for (int i = 0; (ull) i < circles.size(); i++){
                std::set<int> circle = circles[i];
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
                std::set<int> newCircle{strandPairings[j].first, strandPairings[j].second};
                circles.push_back(newCircle);
            }
        }
    }
    std::set<std::set<int>> ret;
    for (auto s : circles) ret.insert(s);
    return ret;
}

PD readPlanarDiagram(int n){ // reads planar diagram, given n crossings
// space-separated
    PD D(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 4; j++){
            std::cin >> D.crossings[i][j];
        }
    }
    return D;
}

PD createPlanarDiagram(std::vector<std::vector<int>> &crossings){ // crossings[i].size() should be 4
    PD D(crossings.size());
    for (int i = 0; i < crossings.size(); i++){
        for (int j = 0; j < 4; j++){
            D.crossings[i][j] = crossings[i][j];
        }
    }
    return D;
}

std::vector<std::vector<std::vector<bool>>> differentialMaps(PD D){
// returns n matrices, each mapping from k 1-resolutions to k+1 1-resolutions for 0 <= k < n
// works with the unreduced Khovanov homology to obtain differentials
// complexity: O(n * 4^n)

    // int n;
    // cout << "Input the number of crossings (no more than 60):" << endl;
    // std::cin >> n;
    // cout << "Input " << 4 * n << " space-separated numbers in planar diagram notation." << endl;

    // PD D(n);
    // for (int i = 0; i < n; i++){
    //     for (int j = 0; j < 4; j++){
    //         std::cin >> D.crossings[i][j];
    //     }
    // }

    // testing if resolutionCircles produces correct circles
    // std::set<std::set<int>> res = resolutionCircles(D, 2);
    // cout << res.size() << endl;

    // construct resolution cube
    int n = D.size();
    std::vector<std::set<std::set<int>>> resolutionCube(1ll << n);
    for (ll i = 0; i < (1ll << n); i++){
        resolutionCube[i] = resolutionCircles(D, i);
        // cerr << i << ": " << resolutionCube[i].size() << endl;
    }

    std::vector<ll> ordering((1ll << n));
    std::vector<ll> circleStartingIndex((1ll << n));
    // ordering[k] takes in a resolution k (binary string) and gives the
    // zero-indexed order for all elements
    // that have the same number of bits as k

    // circleStartingIndex[k] takes in a resolution k (binary string)
    // and gives the starting index for basis elements of a
    // particular circle compared to all circles with the same number of bits
    std::vector<ll> bitCount(n+1, 0);
    std::vector<ll> basisStartCount(n+1, 0);
    for (ll i = 0; i < (1ll << n); i++){
        ordering[i] = bitCount[__builtin_popcountll(i)]++;
        circleStartingIndex[i] = basisStartCount[__builtin_popcountll(i)];
        basisStartCount[__builtin_popcountll(i)] += (1ll << resolutionCube[i].size());
    }
    


    std::vector<std::vector<std::vector<bool>>> differentialMap(n);
    // differentialMap[i] gives the differential map from i 1-resolutions to
    // i+1 1-resolutions
    // differentialMap[i] stores the matrix as a (horizontal) std::vector of column std::vectors
    
    // resize differentialMap matrices
    for (int i = 0; i < n; i++){
        differentialMap[i] = std::vector<std::vector<bool>>(basisStartCount[i], std::vector<bool>(basisStartCount[i+1]));
    }

    std::vector<std::map<std::set<int>, ll>> resolutionCircleIndices(1ll << n);
    // get the circle indices of everything in a resolution

    std::vector<std::vector<std::set<int>>> resolutionCirclesVector(1ll << n);
    // store the circles of a resolution in order in a std::vector

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
                std::set<std::set<int>> oldCircles = resolutionCube[resolution];
                std::set<std::set<int>> newCircles = resolutionCube[newResolution];

                std::set<std::set<int>> oldDiff = setComplement<std::set<int>>(oldCircles, newCircles);
                std::set<std::set<int>> newDiff = setComplement<std::set<int>>(newCircles, oldCircles);

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
                        for (ll oldCirclesIndex = 0; (ull) oldCirclesIndex < oldCircles.size(); oldCirclesIndex++){
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
                        for (ll oldCirclesIndex = 0; (ull) oldCirclesIndex < oldCircles.size(); oldCirclesIndex++){
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

    return differentialMap;
}

std::vector<std::vector<std::vector<bool>>> reducedDifferentialMaps(PD D){
    // returns the reduced differential map of a planar diagram D
    // assumes that strand 1 is the marked strand
    
    // int n; std::cin >> n;
    // PD D = readPlanarDiagram(n);

    int n = D.size();

    // construct resolution cube
    // int n = D.size();
    std::vector<std::set<std::set<int>>> resolutionCube(1ll << n);
    for (ll i = 0; i < (1ll << n); i++){
        resolutionCube[i] = resolutionCircles(D, i);
        // cerr << i << ": " << resolutionCube[i].size() << endl;
    }

    std::vector<ll> ordering((1ll << n));
    std::vector<ll> circleStartingIndex((1ll << n));
    // ordering[k] takes in a resolution k (binary string) and gives the
    // zero-indexed order for all elements
    // that have the same number of bits as k

    // circleStartingIndex[k] takes in a resolution k (binary string)
    // and gives the starting index for basis elements of a
    // particular circle compared to all circles with the same number of bits
    std::vector<ll> bitCount(n+1, 0);
    std::vector<ll> basisStartCount(n+1, 0);
    for (ll i = 0; i < (1ll << n); i++){
        ordering[i] = bitCount[__builtin_popcountll(i)]++;
        circleStartingIndex[i] = basisStartCount[__builtin_popcountll(i)];
        basisStartCount[__builtin_popcountll(i)] += (1ll << (resolutionCube[i].size() - 1));
        // the -1 comes from forstd::cing strand 1 to be labelled as X, no choice -> halves dimension
    }

    std::vector<std::vector<std::vector<bool>>> differentialMap(n);
    // differentialMap[i] gives the differential map from i 1-resolutions to
    // i+1 1-resolutions
    // differentialMap[i] stores the matrix as a (horizontal) std::vector of column std::vectors
    
    // resize differentialMap matrices
    for (int i = 0; i < n; i++){
        differentialMap[i] = std::vector<std::vector<bool>>(basisStartCount[i], std::vector<bool>(basisStartCount[i+1]));
    }

    std::vector<std::map<std::set<int>, ll>> resolutionCircleIndices(1ll << n);
    // get the circle indices of everything in a resolution

    std::vector<std::vector<std::set<int>>> resolutionCirclesVector(1ll << n);
    // store the circles of a resolution in order in a std::vector

    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        ll count = -1;
        for (auto x : resolutionCube[resolution]){
            resolutionCircleIndices[resolution][x] = count++;
            resolutionCirclesVector[resolution].push_back(x);
        }
    }

    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        for (int j = 0; j < n; j++){
            if ((resolution & (1ll << j)) == 0){ // jth bit not yet set
                ll newResolution = resolution | (1ll << j);
                std::set<std::set<int>> oldCircles = resolutionCube[resolution];
                std::set<std::set<int>> newCircles = resolutionCube[newResolution];

                std::set<std::set<int>> oldDiff = setComplement<std::set<int>>(oldCircles, newCircles);
                std::set<std::set<int>> newDiff = setComplement<std::set<int>>(newCircles, oldCircles);

                bool containsX = 0;
                for (auto x : oldDiff){
                    if (x.count(1)) containsX = 1;
                }

                if (!containsX){
                    for (ll oldCirclesSubset = 0; oldCirclesSubset < (1ll << (oldCircles.size()-1)); oldCirclesSubset++){
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
                            for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size() - 1ll; oldCirclesIndex++){
                                if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                    if (newCircles.count(resolutionCirclesVector[resolution][oldCirclesIndex + 1])){
                                        newCircleIndex += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
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
                            for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size()-1ll; oldCirclesIndex++){
                                if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                    if (newCircles.count(resolutionCirclesVector[resolution][oldCirclesIndex + 1])){
                                        newCircleIndex1 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
                                        newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
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
                else { // contains X in the merge/split
                    for (ll oldCirclesSubset = 0; oldCirclesSubset < (1ll << (oldCircles.size()-1)); oldCirclesSubset++){
                        ll oldIndex = circleStartingIndex[resolution] + oldCirclesSubset;
                        if (oldDiff.size() == 2){ // must be a merge
                            ll newCircleIndex = circleStartingIndex[newResolution];
                            for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size() - 1ll; oldCirclesIndex++){
                                if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                    if (!oldDiff.count(resolutionCirclesVector[resolution][oldCirclesIndex + 1])){
                                        newCircleIndex += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
                                    }
                                }
                            }
                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex] = 1;
                        }
                        else if (oldDiff.size() == 1){ // must be a split
                            ll newCircleIndex1 = circleStartingIndex[newResolution];
                            ll newCircleIndex2 = circleStartingIndex[newResolution];
                            for (ll oldCirclesIndex = 0; oldCirclesIndex < oldCircles.size() - 1ll; oldCirclesIndex++){
                                if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                                    newCircleIndex1 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
                                    newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex + 1]]);
                                }
                            }
                            newCircleIndex2 += (1ll << resolutionCircleIndices[newResolution][*(++newDiff.begin())]);
                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex1] = 1;
                            differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex2] = 1;
                        }
                    }
                }
            }
        }
    }

    return differentialMap;

    // output
    // for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
    //     for (ll columnIndex = 0; columnIndex < differentialMap[i][0].size(); columnIndex++){
    //         for (ll rowIndex = 0; rowIndex < differentialMap[i].size(); rowIndex++){
    //             cout << differentialMap[i][rowIndex][columnIndex] << ' ';
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
}

#endif
#define DIFFERENTIAL_MAPS