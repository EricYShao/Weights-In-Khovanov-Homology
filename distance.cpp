// program to find distances when fed matrices from outputDifferentialMaps.cpp
// in integer form

// todo:
// use bitsets and take in matrices as normal

#include <iostream>
#include <vector>
#include <map>
#include <cassert>

using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back

using ll = __int128_t;
using vi = vector<int>;
using vll = vector<ll>;
using pll = pair<ll, ll>;

#define forn(i, n) for (int i = 0; i < int(n); i++)

const ll INF = 1e18;

inline ll next_bit_perm(ll v) { // doesn't work for big v
    ll t = v | (v - 1);
    return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1));
}

void insertVector(ll &rank, vll &basis, ll mask) {
	for (ll i = 0; i < basis.size(); i++) {
		if ((mask & 1ll << i) == 0) continue;
		if (!basis[i]) {
			basis[i] = mask;		
            rank++;
		}
		mask ^= basis[i];
	}
}

bool isLinearlyIndependent(vll &basis, ll mask){
    for (ll i = 0; i < basis.size(); i++){
        if ((mask & 1ll << i) == 0) continue;
        if (!basis[i]) return 1;
        mask ^= basis[i];
    }
    return 0;
}

void outputMinDist(vll &oldMap, vll &newMap){
    vll basis1(127), basis2(127);
    ll rank1, rank2;
    rank1 = rank2 = 0;
    for (auto x : oldMap) insertVector(rank1, basis1, x);
    for (auto x : newMap) insertVector(rank2, basis2, x);
    if (rank1 == newMap.size() - rank2){
        cout << 0 << ' ';
        return;
    }
    vll vectors = newMap;
    ll n = vectors.size();

    for (int k = 1; k <= n; ++k) {
        for (ll w = ((ll)(1) <<k)-1; w < ((ll)(1) << n); w = next_bit_perm(w)) {
            ll mask = 0;
            forn(j, n){
                if (w & ((ll)(1) << j)){
                    mask ^= vectors[j];
                }
            }
            if (!mask){
                if (isLinearlyIndependent(basis1, w)){
                    cout << k << ' ';
                    return;
                }
            }
        }
    }
    
    assert(0);
    return;
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


    freopen("output.txt", "r", stdin);
    int numCases; cin >> numCases;
    vector<vll> matrices(numCases+1);
    forn(i, numCases){
        int n; cin >> n;
        forn(j, n){
            long long x; cin >> x;
            matrices[i+1].push_back(x);
        }
    }
    forn(i, numCases) outputMinDist(matrices[i], matrices[i+1]);
    cout << 1 << endl << "1 ";
    matrices.clear();
    matrices.resize(numCases);
    forn(i, numCases){
        int n; cin >> n;
        forn(j, n){
            long long x; cin >> x;
            matrices[i].push_back(x);
        }
    }
    forn(i, numCases) outputMinDist(matrices[i+1], matrices[i]);

    return 0;
}