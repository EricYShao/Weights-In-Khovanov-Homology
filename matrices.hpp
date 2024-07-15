#ifndef MATRICES
#define MATRICES
#include <vector>
#include <cassert>
using namespace std;
using Matrix = vector<vector<bool>>;
Matrix takeTranspose(Matrix mat){
    Matrix ret(mat[0].size(), vector<bool>(mat.size()));
    for (int i = 0; i < mat.size(); i++)
        for (int j = 0; j < mat[0].size(); j++)
            ret[j][i] = mat[i][j];
    return ret;
}

vector<bool> intToBools(int a, int len){
    vector<bool> ret(len);
    for (int i = 0; i < len; i++){
        if (a & (1 << i)) ret[i] = 1;
    }
    return ret;
}

int boolsToInt(vector<bool> aRep){
    int ret = 0;
    int n = aRep.size();
    for(int i=0; i<n; i++){
        if (aRep[i]) ret += (1 << i);
    }
    return ret;
}

Matrix matrixMult(Matrix A, Matrix B){
    // returns AB
    assert(A[0].size() == B.size());
    Matrix C(A.size(), vector<bool>(B[0].size()));
    for(int i=0; i<A.size(); i++){
        for(int j=0; j<B[0].size(); j++){
            for(int k=0; k<A[0].size(); k++){
                C[i][j] = (C[i][j] != (A[i][k] && B[k][j]));
            }
        }
    }
    return C;
}

pair<bool, vector<int>> invMatGeneral(vector<int> mat, int n){
    // mat should be n x n matrix
    for(int i=0; i<n; i++){
        mat[i] += (1 << (n + i));
    }
    for(int i=0; i<n; i++){
        bool good = 0;
        for(int j = i; j < n; j++){
            if (mat[j] & (1 << i)){
                good = 1;
                swap(mat[i], mat[j]);
                break;
            }
        }
        if (!good){
            return {0, {}};
        }
        for(int j=0; j<n; j++){
            if (j == i) continue;
            if (mat[j] & (1 << i)){
                mat[j] ^= mat[i];
            }
        }
    }
    for(int i=0; i<n; i++){
        mat[i] >>= n;
    }
    return {1, mat};
}

pair<bool, Matrix> invMatBools(Matrix mat){
    Matrix ret(mat.size());
    vector<int> matIntForm(mat.size());
    for(int i=0; i<mat.size(); i++) matIntForm[i] = boolsToInt(mat[i]);
    auto res = invMatGeneral(matIntForm, mat.size());
    if (!(res.first)){
        return {0, {}};
    }
    auto invMatInts = res.second;
    for(int i=0; i<mat.size(); i++){
        ret[i] = intToBools(invMatInts[i], mat.size());
    }
    return {1, ret};
}
#endif