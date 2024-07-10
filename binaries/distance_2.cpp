
// program to find distances using planar diagram notation directly from input.txt
// should modify N to be the maximum matrix size + 10 to be safe

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <bitset>
#include "differentialMaps.hpp"

using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back

const int N = 310;

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;
using pll = pair<ll, ll>;
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

    for (int k = 1; k <= n; ++k) {
        num w(0);
        forn(i, k) w[i] = 1;
        for (; !w[n]; w = nextPerm(w)) {
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
                // cerr << " good " << k << ' ';
                if (1 || isLinearlyIndependent(basis1, w)){
                    // if (k == 15){
                    //     forn(i, N) if (mask[i]) assert(0);
                    //     cerr << endl;
                    //     forn(i, N){
                    //         if (w[i]) cerr << i << ' ';
                    //     }
                    //     cerr << endl;
                    // }
                    return k;
                }
            }
        }
    }
    
    assert(0); // shouldn't reach here
    return 0;
}

signed main()
{

    // cout << "This program finds the minimum distance. Please enter the number of basis elements for the image: " << endl;
    // ll n; cin >> n;
    // while (n >= 120){
    //     cout << "Please enter a number less than 120" << endl;
    //     cin >> n;
    // }
    // cout << "Input " << n << " integers, each representing a basis element in (Z/2Z)^" << n << ':' << endl;
    vector<vector<vector<bool>>> maps = planarDiagramToMaps(1);

    // freopen("output.txt", "r", stdin);
    // int numCases; cin >> numCases;
    ll n = maps.size();
    vector<vn> matrices(n+2, vn(0));
    forn(i, n){
        matrices[i+1] = vn(maps[i].size(), num(0));
        forn(j, maps[i].size()){
            forn(k, maps[i][j].size())
                matrices[i+1][j][k] = maps[i][j][k];
        }
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

    return 0;
}
