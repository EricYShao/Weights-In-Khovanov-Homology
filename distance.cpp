// program to find distances using planar diagram notation directly from input.txt
// should modify N to be the maximum matrix size + 10 to be safe
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
#pragma GCC optimize("O2")
#include "differentialMaps.hpp"
#include "matrices.hpp"

using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back

const int N = 400;
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

int minDist(vn oldMap, vn newMap){
    vn basis1(N), basis2(N);
    ll rank1, rank2;
    rank1 = rank2 = 0;
    for (auto x : oldMap) insertVector(rank1, basis1, x);
    for (auto x : newMap) insertVector(rank2, basis2, x);
    // cerr << "Dimension: " << newMap.size() << ", Rank: " << rank2 << endl;
    if (rank1 == newMap.size() - rank2){
        return 0;
    }
    // cerr << endl << "Homology dimension: " << newMap.size() - rank2 - rank1 << endl;
    vn vectors = newMap;
    ll n = vectors.size();
    ld tic = clock();

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
            forn(j, n){
                if (w[j]){
                    mask ^= vectors[j];
                    good = 1;
                }
            }
            if (!good) break;
            if (mask == num(0)){
                if (isLinearlyIndependent(basis1, w)){
                    return k;
                }
            }
        }
    }
    
    assert(0); // shouldn't reach here
    return 0;
}

void getAllDistances(vector<vector<vector<bool>>> &maps){
    ll n = maps.size();
    vector<vn> matrices(n+2, vn(0));
    forn(i, n){
        matrices[i+1] = vn(maps[i].size(), num(0));
        forn(j, maps[i].size()){
            forn(k, maps[i][j].size())
                matrices[i+1][j][k] = maps[i][j][k];
        }
    }
    ll maxMatrixSize = 0;
    forn(i, n){
        maxMatrixSize = max(maxMatrixSize, (ll) matrices[i+1].size());
    }
    if (maxMatrixSize > N){
        cerr << "N should be at least " << maxMatrixSize + 5 << endl;
        exit(1);
    }
    matrices[n+1] = vn(maps[n-1][0].size());
    forn(i, n+1) cout << minDist(matrices[i], matrices[i+1]) << ' ';
    cout << endl;

    matrices.clear();
    matrices = vector<vn>(n+2, vn(0));

    forn(i, n){
        maps[i] = takeTranspose(maps[i]);
    }

    matrices[0] = vn(maps[0][0].size());
    forn(i, n){
        matrices[i+1] = vn(maps[i].size(), num(0));
        forn(j, maps[i].size()){
            forn(k, maps[i][j].size())
                matrices[i+1][j][k] = maps[i][j][k];
        }
    }
    forn(i, n+1) cout << minDist(matrices[i+1], matrices[i]) << ' ';
}

int main(){
    bool takeAnnular = 1;
    vector<vector<vector<bool>>> maps;
    if (takeAnnular) maps = annular::planarDiagramToMaps();
    else maps = getMaps(getPlanarDiagram(), 1); // always takes reduced homology
    getAllDistances(maps);
}