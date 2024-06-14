#include <iostream>
#include <set>
#include <vector>
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

int main(){
    int n;
    cout << "Input the number of crossings (no more than 60):" << endl;
    cin >> n;
    cout << "Input " << 4 * n << " space-separated numbers in planar diagram notation." << endl;

    PD D(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 4; j++){
            cin >> D.crossings[i][j];
        }
    }

    // testing if resolutionCircles produces correct circles
    // set<set<int>> res = resolutionCircles(D, 2);
    // cout << res.size() << endl;

    // construct resolution cube
    vector<set<set<int>>> resolutionCube(1ll << n);
    for (ll i = 0; i < (1 << n); i++){
        resolutionCube[i] = resolutionCircles(D, i);
        cout << i << ": " << resolutionCube[i].size() << endl;
    }

    return 0;
}