#include <bits/stdc++.h>
#include "differentialMaps.hpp"
using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back
#define forn(i, n) for (int i = 0; i < int(n); i++)

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;


int main(){
    freopen("input.txt", "r", stdin);
    // input format as such:
    // first line contains number of crossings
    // each of the next lines should be a crossing in planar diagram notation
    // then, line containing number of faces
    // then, for each face, line containing the number of edges that bound a face
    // then, line containing the edges that bound a face

    int n; std::cin >> n;
    // assume edges are always 1-indexed
    PD D = readPlanarDiagram(n);
    int f; std::cin >> f;
    vector<vi> faces(f);
    forn(i, f){
        ll numEdges;
        cin >> numEdges;
        forn(j, numEdges){
            ll x; cin >> x;
            faces[i].push_back(x);
        }
    }

    auto differentialMap = annularDifferentialMap(D, faces);

    for (ll i = 0; i < n; i++){ // i = number of 1 resolutions
        for (ll columnIndex = 0; columnIndex < differentialMap[i][0].size(); columnIndex++){
            for (ll rowIndex = 0; rowIndex < differentialMap[i].size(); rowIndex++){
                std::cout << differentialMap[i][rowIndex][columnIndex] << ' ';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}