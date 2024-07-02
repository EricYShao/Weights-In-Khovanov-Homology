#include <bits/stdc++.h>
using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back
#define forn(i, n) for (int i = 0; i < int(n); i++)

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;


vector<vi> perm3 = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};

vector<bool> intToBools(int a, int len){
    vector<bool> ret(len);
    forn(i, len){
        if (a & (1 << i)) ret[i] = 1;
    }
    return ret;
}

int boolsToInt(vector<bool> aRep){
    int ret = 0;
    int n = aRep.size();
    forn(i, n){
        if (aRep[i]) ret += (1 << i);
    }
    return ret;
}

// bool det2(vector<vector<bool>> mat){
//     return (mat[0][0] && mat[1][1]) != (mat[1][0] && mat[0][1]);
// }

bool det3(vector<vector<bool>> mat){
    int ans = 0;
    for (vi p : perm3){
        bool val = 1;
        forn(i, 3){
            val = val && mat[i][p[i]];
        }
        ans += val;
    }
    return (ans % 2);
}

vector<vector<bool>> takeTranspose(vector<vector<bool>> mat){
    vector<vector<bool>> ret(mat[0].size(), vector<bool>(mat.size()));
    for (ll i = 0; i < mat.size(); i++)
        for (ll j = 0; j < mat[0].size(); j++)
            ret[j][i] = mat[i][j];
    return ret;
}


int multPair(int a, int b){
    int ans = 0;
    forn(i, 3){
        forn(j, 3-i){
            if ((a & (1 << i)) && (b & (1 << j)))
                ans ^= (1 << (i + j));
        }
    }
    return ans;
}

int mult(int a){
    // a has 9 bits, represents a sum of elements in the A tensor A
    int ans = 0; // 3 bits
    forn(i, 9){
        if ((1 << i) & a){
            int first = i % 3;
            int second = i / 3;
            ans ^= multPair((1 << first), (1 << second));
        }
    }
    return ans;
}

int coMult(int a){ // big endian convention, A -> A x A
    // ordering goes 11, 1X, 1X^2, X1, XX, XX^2, X^2 1, X^2 X, X^2 X^2
    int ans = 0;
    if (a & 1) ans ^= (1 << 2) + (1 << 4) + (1 << 6); // has a 1 in representation
    if (a & 2) ans ^= (1 << 5) + (1 << 7); // has a X in representation
    if (a & 4) ans ^= (1 << 8); // has X^2 in representation
    return ans;
}

int multDual(int a){ // big endian, A* -> A* x A*
    // ordering goes 1* 1*, 1* X*, 1* X^2*, X* 1*, X* X*, X* X^2*, X^2* 1*, X^2* X*, X^2* X^2*, but with duals
    int ans = 0;
    if (a & 1) ans ^= (1 << 0); // has 1*
    if (a & 2) ans ^= (1 << 1) + (1 << 3); // has X*
    if (a & 4) ans ^= (1 << 2) + (1 << 4) + (1 << 6); // has X^2*
    return ans;
}

int coMultDualPair(int a, int b){
    int ans = 0;
    forn(i, 3){
        forn(j, 3){
            if (i + j >= 2 && (a & (1 << i)) && (b & (1 << j)))
                ans ^= (1 << (i + j - 2));
        }
    }
    return ans;
}

int coMultDual(int a){ // A* x A* -> A*
    // a has 9 bits, represents a sum of elements in the A* tensor A*
    int ans = 0; // 3 bits
    forn(i, 9){
        if ((1 << i) & a){
            int first = i % 3;
            int second = i / 3;
            ans ^= coMultDualPair((1 << first), (1 << second));
        }
    }
    return ans;
}

vector<bool> applyMat(const vector<vector<bool>> &mat, vector<bool> vec){
    int n = vec.size();
    vi ret(n);
    forn(i, n){
        forn(j, n){
            ret[i] += mat[i][j] && vec[j];
        }
    }
    vector<bool> retBool(n);
    forn(i, n) retBool[i] = ret[i] % 2;
    return retBool;
}

int applyMatInt(const vector<vector<bool>> &mat, int x){
    return boolsToInt(applyMat(mat, intToBools(x, mat.size())));
}

int dualA(vi p, int a){ // A_B -> A*_B
    int ans = 0;
    forn(i, 3){
        if (a & (1 << i)){
            ans += (1 << p[i]);
        }
    }
    return ans;
}

int dualAA(vi p, int a){ // A_B x A_B -> A*_B x A*_B
    int ans = 0;
    forn(i, 3){
        forn(j, 3){
            if (a & (1 << (3 * i + j))){
                ans += (1 << (3 * p[i] + p[j]));
            }
        }
    }
    return ans;
}

vi invMatGeneral(vi mat, int n){
    // mat should be n x n matrix
    // assumes invertible
    forn(i, n){
        mat[i] += (1 << (n + i));
    }
    forn(i, n){
        for(int j = i; j < n; j++){
            if (mat[j] & (1 << i)){
                swap(mat[i], mat[j]);
                break;
            }
        }
        forn(j, n){
            if (j == i) continue;
            if (mat[j] & (1 << i)){
                mat[j] ^= mat[i];
            }
        }
    }
    forn(i, n){
        mat[i] >>= n;
    }
    return mat;
}

vector<vector<bool>> invMatBools(vector<vector<bool>> mat){
    vector<vector<bool>> ret(mat.size());
    vi matIntForm(mat.size());
    forn(i, mat.size()) matIntForm[i] = boolsToInt(mat[i]);
    vi invMatInts = invMatGeneral(matIntForm, mat.size());
    forn(i, mat.size()){
        ret[i] = intToBools(invMatInts[i], mat.size());
    }
    return ret;
}

void check(vi basis){ // should be 3 integers each corresponding to a basis element
    // iterate over all permutations composed with taking the dual

    vector<vector<bool>> A_BtoA(3, vector<bool>(3));
    forn(i, 3){
        forn(j, 3){
            if ((1 << j) & basis[i]) A_BtoA[j][i] = 1;
        }
    }
    if (!det3(A_BtoA)) return; // not a basis, not invertible, det = 0
    auto AtoA_B = invMatBools(A_BtoA);

    // create matrix A_BA_BtoAA
    vector<vector<bool>> A_BA_BtoAA(9, vector<bool>(9));
    forn(i, 3){
        forn(j, 3){
            // considering where element basis[i] x basis[j] maps
            // turns on 1s in column 3*i + j
            forn(ii, 3){ // bit index in i
                forn(jj, 3){ // bit index in j
                    if ((basis[i] & (1 << ii)) && (basis[j] & (1 << jj))){
                        A_BA_BtoAA[3*ii+jj][3*i+j] = !A_BA_BtoAA[3*ii+jj][3*i+j];
                    }
                }
            }
        }
    }

    vector<vector<bool>> AAtoA_BA_B = invMatBools(A_BA_BtoAA);

    for (auto p : perm3){
        // check A -> A* tensor A*
        vi originalNums = {1, 2, 4};
        vi res1(3), res2(3);
        forn(i, 3){
            // debug stuff commented out
            // if (basis[0] == 1 && basis[1] == 2 && basis[2] == 5 && originalNums[i] == 5 && p[0] == 2 && p[1] == 1)
            //     cout << "here";
            // int x = originalNums[i];
            // x = applyMatInt(AtoA_B, x);
            // x = dualA(p, x);
            // x = applyMatInt(takeTranspose(AtoA_B), x);
            // x = multDual(x);
            // x = applyMatInt(takeTranspose(A_BA_BtoAA), x);
            res1[i] = multDual(applyMatInt(takeTranspose(AtoA_B), dualA(p, applyMatInt(AtoA_B, originalNums[i])))); // up first
            // debug stuff commented out            
            // int y = originalNums[i];
            // y = coMult(y);
            // y = applyMatInt(AAtoA_BA_B, y);
            // y = dualAA(p, y);
            // y = applyMatInt(takeTranspose(AAtoA_BA_B), y);
            res2[i] = applyMatInt(takeTranspose(AAtoA_BA_B), dualAA(p, applyMatInt(AAtoA_BA_B, coMult(originalNums[i])))); // right first
        }

        // check A x A -> A*
        originalNums = {1, 2, 4, 8, 16, 32, 64, 128, 256};
        vi res3(9), res4(9);
        forn(i, 9){
            res3[i] = coMultDual(applyMatInt(takeTranspose(AAtoA_BA_B), dualAA(p, applyMatInt(AAtoA_BA_B, originalNums[i])))); // up first
            // debug stuff commented out
            // int x = originalNums[i];
            // x = mult(x);
            // x = applyMatInt(AtoA_B, x);
            // x = dualA(p, x);
            // x = applyMatInt(A_BtoA, x);
            res4[i] = applyMatInt(takeTranspose(AtoA_B), dualA(p, applyMatInt(AtoA_B, mult(originalNums[i]))));
        }
        if ((res1 == res2) && (res3 == res4)) cout << "found: " << basis[0] << ' ' << basis[1] << ' ' << basis[2] << ": "
        << p[0] << p[1] << p[2] << endl;
    }
    return;
}


int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    // vector<vector<bool>> mat = {{1, 0, 0, 1}, {1, 0, 1, 0}, {0, 1, 1, 1}, {0, 0, 1, 0}};
    // auto x = invMatBools(mat);
    for (int i = 0; i < 8; i++){
        for(int j = i+1; j < 8; j++){
            for(int k = j+1; k < 8; k++){
                check({i, j, k});
            }
        }
    }
    // check({1, 2, 4});
    // cout << mult(1);
    // forn(i, 9){
    //      cout << (1 << i) << ": " << mult(1 << i) << endl;
    // }
    // cout << mult(1) << ' ' << mult(2) << ' ' << mult(6) << endl;
    // forn(i, 8){
    //     cout << i << ": " << coMult(i) << endl; // coMult is correct
    // }
    // cout << coMult(1) << ' ' << coMult(2) << ' ' << coMult(4) << endl;
    // forn(i, 8){
    //     cout << i << ": " << multDual(i) << endl; // multDual is correct
    // }
    // cout << multDual(1) << ' ' << multDual(2) << ' ' << multDual(4) << endl;
    // forn(i, 9){
    //     cout << (1 << i) << ": " << coMultDual(1 << i) << endl;
    // }
    // cout << coMultDual(1) << ' ' << coMultDual(2) << ' ' << coMultDual(4) << endl;
    return 0;
}