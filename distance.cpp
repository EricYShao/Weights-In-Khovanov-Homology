// program to find distances using planar diagram notation directly from input.txt
// should modify N to be the maximum matrix size + 10 to be safe

// for regular reduced homology, input format should contain 4*n space separated integers
// each referring to a strand index in a crossing

// for annular, input format should be as follows in input.txt
// line 1: contains two space-separated numbers n, f. n is the the number of crossings, f the number of faces
// lines 2 to n+1: crossing information in planar diagram notation
// line n+2 to n+f+1: the first number contains the number of mini-strands s that bound a face. On the same line, there are s more numbers, each denoting a mini-strand index
// the program assumes that line n+2 describes the face with the puncture

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <bitset>

#include "differentialMaps.hpp"
#include "matrices.hpp"
#include "sl3Calculation.hpp"

using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back

const int N = 401;
const int timeLimit = 30; // seconds to finish calculation at one degree

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;
using pll = pair<ll, ll>;
using ld = long double;
using num = bitset<N>;
using vn = vector<num>;

#define forn(i, n) for (int i = 0; i < int(n); i++)

const ll INF = 1e18;

num nextPerm(num v){
    num ret = v;
    ll count = 0;
    forn(i, N-1){
        if (v[i]){
            if (!v[i+1]){
                ret[i] = 0;
                ret[i+1] = 1;
                return ret;
            }
            else{
                ret[i] = 0;
                ret[count++] = 1;
            }
        }
    }
    ret = num(0);
    return ret;
}

// inline ll next_bit_perm(ll v) { // doesn't work for big v
//     ll t = v | (v - 1);
//     return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1));
// }

void insertVector(ll &rank, vn &basis, num mask) {
	for (ll i = 0; i < basis.size(); i++) {
		if (!(mask[i])) continue;
		if (basis[i] == num(0)) {
			basis[i] = mask;		
            rank++;
		}
		mask ^= basis[i];
	}
}

bool isLinearlyIndependent(const vn &basis, num mask){
    for (ll i = 0; i < basis.size(); i++){
        if (!(mask[i])) continue;
        if (basis[i] == num(0)) return 1;
        mask ^= basis[i];
    }
    return 0;
}

int minDist(vn oldMap, vn newMap, bool outputLengths, bool outputHomologyDimensions, bool outputCount){
    vn basis1(N), basis2(N);
    ll rank1, rank2;
    rank1 = rank2 = 0;
    for (auto x : oldMap) insertVector(rank1, basis1, x);
    for (auto x : newMap) insertVector(rank2, basis2, x);
    // cerr << "Dimension: " << newMap.size() << ", Rank: " << rank2 << endl;
    if (outputLengths){
        return newMap.size();
    }
    if (outputHomologyDimensions){
        return newMap.size() - rank2 - rank1;
    }
    if (rank1 == newMap.size() - rank2){
        return 0;
    }
    vn vectors = newMap;
    ll n = vectors.size();
    ld tic = clock();

    ll count = 0;

    for (int k = 1; k <= n; ++k) {
        num w(0);
        forn(i, k) w[i] = 1;
        for (; !w[n]; w = nextPerm(w)) {
            ld tac = clock();
            if ((tac - tic) / CLOCKS_PER_SEC > timeLimit){
                return -k;
            }
            num mask = 0;
            bool good = 0;
            forn(j, n-14){
                if (w[j]){
                    mask ^= vectors[j+14];
                    good = 1;
                }
            }
            if (!good) break;
            if (mask == num(0)){
                num newW;
                forn(i, N){
                    if (w[i]) newW[i+14] = 1;
                }
                if (isLinearlyIndependent(basis1, newW)){
                    if (outputCount) count++;
                    else return k;
                }
            }
        }
        if (count && outputCount) return count;
    }
    
    assert(0); // shouldn't reach here
    return 0;
}

void getAllDistances(vector<Matrix> &maps, bool outputLengths, bool outputHomologyDimension, bool outputDistance, bool outputCounts){
    ll maxMatrixSize = 0;
    forn(i, maps.size()){
        maxMatrixSize = max(maxMatrixSize, (ll)(maps[i].r));
    }
    if (maxMatrixSize + 5 > N){
        cerr << "N should be at least " << maxMatrixSize + 5 << endl;
        exit(1);
    }

    ll n = maps.size();
    // forn(i, n) outputMatrix(maps[i]);
    forn(i, n-1){
        if (maps[i].r == 0 || maps[i].c == 0 || maps[i+1].r == 0 || maps[i+1].c == 0) continue;
        Matrix mat = matrixMult(maps[i], maps[i+1]);
        // outputMatrix(mat);
        for (auto x : mat.mat){
            for (auto y : x){
                assert(!y);
            }
        }
    }
    vector<vn> matrices(n+2, vn(0));
    forn(i, n){
        matrices[i+1] = vn(maps[i].size(), num(0));
        forn(j, maps[i].size()){
            forn(k, maps[i][j].size())
                matrices[i+1][j][k] = maps[i][j][k];
        }
    }
    if (maps[n-1].size() == 0) matrices[n+1] = {};
    else matrices[n+1] = vn(maps[n-1][0].size());
    vector<vn> matrixTransposes(n+2, vn(0));
    forn(i, n){
        maps[i] = takeTranspose(Matrix(maps[i]));
    }

    if (maps[0].size() == 0) matrixTransposes[0] = vn(0);
    else matrixTransposes[0] = vn(maps[0][0].size());
    forn(i, n){
        matrixTransposes[i+1] = vn(maps[i].size(), num(0));
        forn(j, maps[i].size()){
            forn(k, maps[i][j].size())
                matrixTransposes[i+1][j][k] = maps[i][j][k];
        }
    }
    if (outputLengths){
        cout << "Lengths:" << endl;
        forn(i, n+1) cout << minDist(matrices[i], matrices[i+1], 1, 0, 0) << ' ';
        cout << endl;
        forn(i, n+1) cout << minDist(matrixTransposes[i+1], matrixTransposes[i], 1, 0, 0) << ' ';
        cout << endl;
    }
    if (outputHomologyDimension){
        cout << "Homology:" << endl;
        forn(i, n+1) cout << minDist(matrices[i], matrices[i+1], 0, 1, 0) << ' ';
        cout << endl;
        forn(i, n+1) cout << minDist(matrixTransposes[i+1], matrixTransposes[i], 0, 1, 0) << ' ';
        cout << endl;
    }
    if (outputDistance){
        cout << "Distances:" << endl;
        forn(i, n+1) cout << minDist(matrices[i], matrices[i+1], 0, 0, 0) << ' ';
        cout << endl;
        forn(i, n+1) cout << minDist(matrixTransposes[i+1], matrixTransposes[i], 0, 0, 0) << ' ';
        cout << endl;
    }
    if (outputCounts){
        cout << "Number of Minimially Weighted Elements:" << endl;
        forn(i, n+1) cout << minDist(matrices[i], matrices[i+1], 0, 0, 1) << ' ';
        cout << endl;
        forn(i, n+1) cout << minDist(matrixTransposes[i+1], matrixTransposes[i], 0, 0, 1) << ' ';
        cout << endl;
    }
}

int main(){
    bool takeAnnular = 1;
    bool restrictAnnularGrading = 1;
    vector<Matrix> maps;
    if (takeAnnular) maps = annular::planarDiagramToMaps(restrictAnnularGrading);
    else maps = getMaps(getPlanarDiagram(), 1); // always takes reduced homology

    // maps = getMatrices();
    getAllDistances(maps, 1, 1, 1, 1);
}