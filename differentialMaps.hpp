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

std::set<std::set<int>> resolutionCircles(PD diagram, ll resolution){
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

std::vector<std::vector<bool>> takeTranspose(std::vector<std::vector<bool>> initial){
    std::vector<std::vector<bool>> ret(initial[0].size(), std::vector<bool>(initial.size()));
    for (ll i = 0; i < initial.size(); i++)
        for (ll j = 0; j < initial[0].size(); j++)
            ret[j][i] = initial[i][j];
    return ret;
}

std::vector<std::vector<std::vector<bool>>> planarDiagramToMaps(bool reducedHomology){
    // reads planar diagram notation from input.txt and returns the differential maps
    freopen("input.txt", "r", stdin);
    std::vector<std::vector<int>> input;
    int x;
    while (std::cin >> x){
        if (!input.size() || input[input.size()-1].size() == 4){
            input.push_back({x});
        }
        else{
            input[input.size()-1].push_back(x);
        }
    }

    PD D = createPlanarDiagram(input);
    
    std::vector<std::vector<std::vector<bool>>> maps;
    if (reducedHomology) maps = reducedDifferentialMaps(D);
    else maps = differentialMaps(D);

    return maps;
}


// annular stuff below

void insertVector(std::vector<ll> &basis, ll mask) {
	for (ll i = 0; i < basis.size(); i++) {
		if (!(mask & (1ll << i))) continue;
		if (!basis[i]) {
			basis[i] = mask;
		}
		mask ^= basis[i];
	}
}

bool isLinearlyIndependent(const std::vector<ll> &basis, ll mask){
    for (ll i = 0; i < basis.size(); i++){
        if (!(mask & (1ll << i))) continue;
        if (!basis[i]) return 1;
        mask ^= basis[i];
    }
    return 0;
}

bool containsPuncture(const std::set<int> &circle, const std::vector<ll> &basis1, const std::vector<ll> &basis2){
    ll mask = 0;
    for (auto x : circle){
        mask += (1ll << (x-1));
    }
    assert(!isLinearlyIndependent(basis2, mask));
    return isLinearlyIndependent(basis1, mask);
}

std::vector<bool> mapVVtoA(bool first, bool second){
    if (first ^ second) return {0, 1}; // means assign V-
    return {};
}

bool mapVAtoV(bool first, bool second){ 
    return first;
}

std::vector<bool> annularMerge(const std::set<int> &circle1, bool circle1Status, const std::set<int> &circle2, bool circle2Status, const std::vector<ll> &basis1, const std::vector<ll> &basis2){
    bool circle1HasPuncture = containsPuncture(circle1, basis1, basis2);
    bool circle2HasPuncture = containsPuncture(circle2, basis1, basis2);
    if (circle1HasPuncture && circle2HasPuncture){
        return mapVVtoA(circle1Status, circle2Status);
    }
    if (!circle1HasPuncture && !circle2HasPuncture){ // regular Audoux map
        return {circle1Status != circle2Status};
    }
    if (circle2HasPuncture) std::swap(circle1Status, circle2Status);
    return {mapVAtoV(circle1Status, circle2Status)};
}

std::vector<std::pair<bool, bool>> mapAtoVV(bool first){ // returns a sum
    if (first) return {{1, 0}, {0, 1}};
    return {};
}

std::vector<std::pair<bool, bool>> mapVtoVA(bool first){
    return {{first, 0}, {first, 1}};
}

std::vector<std::pair<bool, bool>> annularSplit(const std::set<int> &circle, bool circleStatus, const std::set<int> &res1, const std::set<int> &res2, const std::vector<ll> &basis1, const std::vector<ll> &basis2){
    bool origHasPuncture = containsPuncture(circle, basis1, basis2);
    bool res1HasPuncture = containsPuncture(res1, basis1, basis2);
    bool res2HasPuncture = containsPuncture(res2, basis1, basis2);
    if (!origHasPuncture && !res1HasPuncture && !res2HasPuncture){ // regular split
        if (circleStatus){ // (+) -> (-)(-) + (+)(+)
            return {{0, 0}, {1, 1}};
        }
        else{ // (-) -> (-)(+) + (+)(-)
            return {{0, 1}, {1, 0}};
        }
    } 
    if (origHasPuncture){
        if (res1HasPuncture){
            return mapVtoVA(circleStatus);
        }
        else if (res2HasPuncture){
            auto ret = mapVtoVA(circleStatus);
            std::swap(ret[0].first, ret[0].second);
            std::swap(ret[1].first, ret[1].second);
            return ret;
        }
        assert(0); // split configuration invalid
    }
    else{
        if (res1HasPuncture && res2HasPuncture){
            return {{1, 0}, {0, 1}};
        }
        assert(0); // split configuration invalid
    }
    assert(0); // split configuration invalid
    return {};
}

std::vector<std::vector<std::vector<bool>>> annularDifferentialMap(PD D, std::vector<std::vector<int>> faces){
    // freopen("input.txt", "r", stdin);
    // int n; std::cin >> n;
    // assume edges are always 1-indexed
    // PD D = readPlanarDiagram(n);
    int n = D.size();
    // int f; std::cin >> f;

    std::vector<ll> basis1(63), basis2(63);
    ll specialMask;
    for (int i = 0; i < faces.size(); i++){
        ll mask = 0;
        for (int j = 0; j < faces[i].size(); j++){
            ll x = faces[i][j];
            mask += (1ll << (x-1));
        }
        if (!i) specialMask = mask;
        else insertVector(basis1, mask);
    }
    basis2 = basis1;
    insertVector(basis2, specialMask);

    std::vector<std::set<std::set<int>>> resolutionCube(1ll << n);
    for (ll i = 0; i < (1ll << n); i++){
        resolutionCube[i] = resolutionCircles(D, i);
        // test puncture detection
        // for (auto circle : resolutionCube[i]){
        //     if (containsPuncture(circle, basis1, basis2)){
        //         for (auto x : circle) cerr << x << ' ';
        //         cerr << endl;
        //     }
        // }
    }
    
    std::vector<ll> ordering((1ll << n));
    std::vector<ll> circleStartingIndex((1ll << n));

    std::vector<ll> bitCount(n+1, 0);
    std::vector<ll> basisStartCount(n+1, 0);
    for (ll i = 0; i < (1ll << n); i++){
        ordering[i] = bitCount[__builtin_popcountll(i)]++;
        circleStartingIndex[i] = basisStartCount[__builtin_popcountll(i)];
        basisStartCount[__builtin_popcountll(i)] += (1ll << resolutionCube[i].size());
    }


    std::vector<std::vector<std::vector<bool>>> differentialMap(n);
    for (int i = 0; i < n; i++){
        differentialMap[i] = std::vector<std::vector<bool>>(basisStartCount[i], std::vector<bool>(basisStartCount[i+1]));
    }

    std::vector<std::map<std::set<int>, ll>> resolutionCircleIndices(1ll << n);
    std::vector<std::vector<std::set<int>>> resolutionCirclesVector(1ll << n);

    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        ll count = 0;
        for (auto x : resolutionCube[resolution]){
            resolutionCircleIndices[resolution][x] = count++;
            resolutionCirclesVector[resolution].push_back(x);
        }
    }


    for (ll resolution = 0; resolution < (1ll << n); resolution++){
        for (int j = 0; j < n; j++){
            if ((resolution & (1ll << j)) != 0) continue; // jth bit already set
            // jth bit not yet set
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

                    bool circleOneStatus = ((oldCirclesSubset 
                    & (1ll << resolutionCircleIndices[resolution][*(oldDiff.begin())])) != 0);
                    bool circleTwoStatus = ((oldCirclesSubset 
                    & (1ll << resolutionCircleIndices[resolution][*(++oldDiff.begin())])) != 0);

                    std::vector<bool> newCircleStatuses = annularMerge(*oldDiff.begin(), circleOneStatus,
                    *(++oldDiff.begin()), circleTwoStatus, basis1, basis2);

                    ll newCircleStartIndex = circleStartingIndex[newResolution];
                    // update index for the new merged circle
                    for (bool newCircleStatus : newCircleStatuses){
                        ll newCircleIndex = newCircleStartIndex;
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
                }
                else if (oldDiff.size() == 1){ // must be a split
                    ll newCircleStartIndex = circleStartingIndex[newResolution];
                    for (ll oldCirclesIndex = 0; (ull) oldCirclesIndex < oldCircles.size(); oldCirclesIndex++){
                        if (oldCirclesSubset & (1ll << oldCirclesIndex)){
                            if (newCircles.count(resolutionCirclesVector[resolution][oldCirclesIndex])){
                                newCircleStartIndex += (1ll << resolutionCircleIndices[newResolution][resolutionCirclesVector[resolution][oldCirclesIndex]]);
                            }
                        }
                    }

                    bool oldCircleStatus = oldCirclesSubset & (1ll << resolutionCircleIndices[resolution][*oldDiff.begin()]);
                    std::vector<std::pair<bool, bool>> newCircleStatuses = annularSplit(*oldDiff.begin(), oldCircleStatus,
                    *newDiff.begin(), *(++newDiff.begin()), basis1, basis2);

                    int newCircle1Index = resolutionCircleIndices[newResolution][*newDiff.begin()];
                    int newCircle2Index = resolutionCircleIndices[newResolution][*(++newDiff.begin())];
                    for (std::pair<bool, bool> newCircleStatus : newCircleStatuses){
                        ll newCircleIndex = newCircleStartIndex;
                        if (newCircleStatus.first) newCircleIndex += (1ll << newCircle1Index);
                        if (newCircleStatus.second) newCircleIndex += (1ll << newCircle2Index);
                        differentialMap[__builtin_popcountll(resolution)][oldIndex][newCircleIndex] = 1;
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
    //             std::cout << differentialMap[i][rowIndex][columnIndex] << ' ';
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    // return 0;
}

#endif
#define DIFFERENTIAL_MAPS