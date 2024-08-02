#include <bits/stdc++.h>
#include "matrices.hpp"
using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back
#define forn(i, n) for (int i = 0; i < int(n); i++)

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;

// need a function for evaluating bubble bursting

bool bubbleBurst2(vi &dots){ // will modify dots array
// evaluates closed foam with dots[i] dots on the ith facet
    while (dots.size() > 1){
        if (dots[dots.size()-1] >= 3 || dots[dots.size()-2] >= 3 || dots[dots.size()-3] >= 3) return 0;
        if (!dots[dots.size()-1] && !dots[dots.size()-2]) return 0;
        if (dots[dots.size()-1] == 1 && dots[dots.size()-2] == 1) return 0;
        if (dots[dots.size()-1] && dots[dots.size()-2] && dots[dots.size()-3]) return 0;
        if (dots[dots.size()-1] == 2 && dots[dots.size()-2] == 2) return 0;
        dots[dots.size()-3] += dots[dots.size()-1] + dots[dots.size()-2] - 1;
        dots.pop_back(); dots.pop_back();   
    }
    return (dots[0] == 2);
}

bool bubbleBurst(vi dots){
    forn(i, dots.size()) if (dots[i] >= 3) return 0;
    ll sum = 0;
    forn(i, dots.size()) sum += dots[i];
    assert(dots.size() % 2);
    if (sum != (dots.size()/2) + 2) return 0;
    return bubbleBurst2(dots);
}

bool bubbleBurst3(vi dotMasks){ // dots given as integers, 2^x bit toggled -> include sum with x dots
    forn(i, dotMasks.size()) dotMasks[i] %= 8; // restricts to stuff with at most 2 dots
    assert(dotMasks.size() % 2);
    while (dotMasks.size() > 1){
        int countLast = dotMasks[dotMasks.size() - 1];
        int countSecondLast = dotMasks[dotMasks.size() - 2];
        // use (X, 1, 2) -> X, (X, 2, 1) -> X, (X, 1, 4) -> 2X (mod 8), (X, 4, 1) -> 2X (mod 8),
        // (X, 4, 2) -> 4X (mod 8)
        int nextLast = dotMasks[dotMasks.size() - 3];
        int curNextVal = 0;
        forn(i, 3){
            forn(j, 3){
                if (i != j && ((1 << i) & countLast) && ((1 << j) & countSecondLast)){
                    curNextVal ^= ((nextLast) * (1 << (i + j - 1))) % 8;
                }
            }
        }
        dotMasks.pop_back();
        dotMasks.pop_back();
        dotMasks[dotMasks.size()-1] = curNextVal;
    }
    return dotMasks[0] & (1 << 2);
}

vector<vi> getBasis(int length){
    vector<vi> basis;
    for (ll i = 0; i < 3 * (1ll << length); i++){
        vi dots;
        bool flip = 0;
        ll cpy = i;
        for (int j = 0; j < length; j++){
            if (!flip){
                dots.pb(cpy % 2);
                dots.pb(0);
            }
            else{
                dots.pb(0);
                dots.pb(cpy % 2);
            }
            cpy /= 2;
            flip = !flip;
        }
        dots.pb(cpy);
        reverse(all(dots));
        basis.pb(dots);
    }
    return basis;
}

vector<vi> getDualBasis(int length){
    vector<vi> dualBasis;
    for (ll i = 0; i < 3 * (1ll << length); i++){
        vi dualDots;
        bool flip = 0;
        ll cpy = i;
        for (int j = 0; j < length; j++){
            if (!flip){
                dualDots.pb(0);
                dualDots.pb(1 - (cpy % 2));
            }
            else{
                dualDots.pb(1 - (cpy % 2));
                dualDots.pb(0);
            }
            cpy /= 2;
            flip = !flip;
        }
        dualDots.pb(2-cpy);
        reverse(all(dualDots));
        dualBasis.pb(dualDots);
    }
    return dualBasis;
}

Matrix changeBasis(vector<vi> newBasis, int thetaLength){
    // get map from theta_newbasis to theta_standardBasis
    vector<vi> thetaDualBasis = getDualBasis(thetaLength);
    assert(thetaDualBasis.size() == newBasis.size());
    Matrix ret(newBasis.size(), vector<bool>(newBasis.size()));
    forn(i, thetaDualBasis.size()){
        forn(j, newBasis.size()){
            vi dotMasks = newBasis[j];
            forn(k, thetaDualBasis[i].size()){
                dotMasks[k] *= (1 << thetaDualBasis[i][k]);
                dotMasks[k] %= 8;
            }
            ret[i][j] = bubbleBurst3(dotMasks);
        }
    }
    return ret;
}

Matrix thetaSplit(int origLength, int rightLength){ // assumes rightmost is ccw
    // origLength is number of initial "rungs" on the ladder
    // final length is number of final "rungs" on the rightmost ladder
    // two theta lengths are finalLength, origLength - 1 - finalLength
    vector<vi> basis = getBasis(origLength);
    vector<vi> rightDual = getDualBasis(rightLength);
    vector<vi> leftDual = getDualBasis(origLength - 1 - rightLength);
    Matrix retMat(leftDual.size() * rightDual.size(), vector<bool>(basis.size()));

    for(ll i = 0; i < basis.size(); i++){
        for (ll j = 0; j < leftDual.size() * rightDual.size(); j++){
            ll rightIndex = j % rightDual.size();
            ll leftIndex = j / rightDual.size();

            int leftLength = origLength - 1 - rightLength;
            vi dots = basis[i];
            forn(ii, 2 * leftLength+1) dots[ii] += leftDual[leftIndex][ii];
            forn(ii, 2 * rightLength+1) dots[2 * leftLength + 2 + ii] += rightDual[rightIndex][ii];
            retMat[j][i] = bubbleBurst(dots);
        }
    }
    return retMat;
}

Matrix thetaMerge(int leftLength, int rightLength){
    int finalLength = leftLength + rightLength + 1;
    vector<vi> leftBasis = getBasis(leftLength);
    vector<vi> rightBasis = getBasis(rightLength);
    vector<vi> dualBasis = getDualBasis(finalLength);

    Matrix retMat(dualBasis.size(), vector<bool>(leftBasis.size() * rightBasis.size()));
    
    for (ll i = 0; i < leftBasis.size() * rightBasis.size(); i++){
        ll leftIndex = i / rightBasis.size();
        ll rightIndex = i % rightBasis.size();

        for (ll j = 0; j < dualBasis.size(); j++){
            vi dots = dualBasis[j];
            forn(ii, 2 * leftLength+1) dots[ii] += leftBasis[leftIndex][ii];
            forn(ii, 2 * rightLength+1) dots[2 * leftLength + 2 + ii] += rightBasis[rightIndex][ii];
            retMat[j][i] = bubbleBurst(dots);
        }
    }
    return retMat;
}

Matrix kroneckerProduct(Matrix A, Matrix B){
    Matrix C(A.r * B.r, vector<bool>(A.c * B.c));
 
    for (int i = 0; i < A.r; i++)
        for (int j = 0; j < A.c; j++)
            for (int k = 0; k < B.r; k++)
                for (int l = 0; l < B.c; l++)
                    C[i * B.r + k][j * B.c + l] = A[i][j] && B[k][l];
 
    return C;
}

Matrix stackVertical(Matrix top, Matrix bot){
    Matrix C(top.r + bot.r, vector<bool>(top.c));
    forn(i, top.r){
        forn(j, top.c){
            C[i][j] = top[i][j];
        }
    }
    forn(i, bot.r){
        forn(j, bot.c){
            C[top.size() + i][j] = bot[i][j];
        }
    }
    return C;
}

Matrix stackVertical(vector<Matrix> matrices){
    // stacks first at the top, last at the bottom
    assert(matrices.size());
    if (matrices.size() == 1) return matrices[0];
    Matrix last = matrices.back();
    matrices.pop_back();
    return stackVertical(stackVertical(matrices), last);
}

Matrix stackHorizontal(Matrix left, Matrix right){
    Matrix C(left.r, vector<bool>(left.c + right.c));
    forn(i, left.r){
        forn(j, left.c){
            C[i][j] = left[i][j];
        }
        forn(k, right.c){
            C[i][left.c + k] = right[i][k];
        }
    }
    return C;
}

Matrix stackHorizontal(vector<Matrix> matrices){
    // stacks first at the left, last at the right
    assert(matrices.size());
    if (matrices.size() == 1) return matrices[0];
    Matrix last = matrices.back();
    matrices.pop_back();
    return stackHorizontal(stackHorizontal(matrices), last);
}

vector<Matrix> getMatrices()
// int main()
{
//     vector<vi> test = {{2, 1, 0, 1, 0}, {2, 0, 1, 0, 1}, {2, 1, 0, 0, 1}, {2, 0, 1, 1, 0}};
//     forn(i, 4) assert(bubbleBurst(test[i]));
//     vector<vi> test2 = {{4, 2, 1, 2, 1}, {4, 1, 2, 1, 2}, {4, 2, 1, 1, 2}, {4, 1, 2, 2, 1}, {4, 3, 3}};
//     assert(bubbleBurst3(test2[0]));
//     assert(bubbleBurst3(test2[1]));
//     assert(bubbleBurst3(test2[2]));
//     assert(bubbleBurst3(test2[3]));
//     assert(bubbleBurst3(test2[4]) == 0);

//     cout << "Asserts passed";

    /*
    auto merge0 = thetaMerge(0, 0);
    auto split1 = thetaSplit(1, 0);

    

   
    auto topLeftMap = kroneckerProduct(I3, split1);

    auto botLeftMap = thetaMerge(0, 1);

    auto leftMap = stackVertical(topLeftMap, botLeftMap);

    //outputMatrix(merge0);

    auto topRightMap = kroneckerProduct(merge0, I3); // sketchy
    //outputMatrix(topRightMap); cout << endl;
    auto botRightMap = thetaSplit(2, 0);
    //outputMatrix(botRightMap); cout << endl;

    auto rightMap = stackHorizontal(topRightMap, botRightMap);
    //outputMatrix(rightMap);

    // auto res = matrixMult(rightMap, leftMap);
    // outputMatrix(res);
    // outputMatrix(rightMap);
    */

    // Matrix firstBasisChange = changeBasis({{1}, {2}, {5}}, 0);
    // Matrix left = kroneckerProduct(firstBasisChange, firstBasisChange);
    // Matrix middle = thetaMerge(0, 0);
    // Matrix right = invMatBools(changeBasis({{1, 1, 1}, {1, 2, 1}, {2, 1, 1}, {2, 2, 1}, {5, 1, 1}, {5, 2, 1}}, 1)).second;
    // Matrix diff = matrixMult(right, matrixMult(middle, left));

    // Matrix left = changeBasis({{1, 1, 1}, {1, 2, 1}, {2, 1, 1}, {2, 2, 1}, {5, 1, 1}, {5, 2, 1}}, 1);
    // Matrix middle = thetaSplit(1, 0);
    // Matrix changeBack = invMatBools(changeBasis({{1}, {2}, {5}}, 0)).second;
    // Matrix right = kroneckerProduct(changeBack, changeBack);
    // vector<Matrix> maps = {matrixMult(right, matrixMult(middle, left))};
    
    
    Matrix I3({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    Matrix I6(6, 6);
    forn(i, 6) I6[i][i] = 1;
    
    vector<vi> basis1 = {{1}, {2}, {5}};
    vector<vi> basis2 = {{1, 1, 1}, {1, 1, 2}, {2, 1, 1}, {2, 1, 2}, {4, 1, 1}, {4, 1, 2}};
    vector<vi> basis3 = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 2}, {1, 1, 2, 1, 1}, {1, 1, 2, 1, 2},
                        {2, 1, 1, 1, 1}, {2, 1, 1, 1, 2}, {2, 1, 2, 1, 1}, {2, 1, 2, 1, 2},
                        {4, 1, 1, 1, 1}, {4, 1, 1, 1, 2}, {4, 1, 2, 1, 1}, {4, 1, 2, 1, 2}};
    vector<vi> basis4 = {{1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 2}, {1, 1, 1, 1, 2, 1, 1}, {1, 1, 1, 1, 2, 1, 2},
                        {1, 1, 2, 1, 1, 1, 1}, {1, 1, 2, 1, 1, 1, 2}, {1, 1, 2, 1, 2, 1, 1}, {1, 1, 2, 1, 2, 1, 2},
                        {2, 1, 1, 1, 1, 1, 1}, {2, 1, 1, 1, 1, 1, 2}, {2, 1, 1, 1, 2, 1, 1}, {2, 1, 1, 1, 2, 1, 2},
                        {2, 1, 2, 1, 1, 1, 1}, {2, 1, 2, 1, 1, 1, 2}, {2, 1, 2, 1, 2, 1, 1}, {2, 1, 2, 1, 2, 1, 2},
                        {5, 1, 1, 1, 1, 1, 1}, {5, 1, 1, 1, 1, 1, 2}, {5, 1, 1, 1, 2, 1, 1}, {5, 1, 1, 1, 2, 1, 2},
                        {5, 1, 2, 1, 1, 1, 1}, {5, 1, 2, 1, 1, 1, 2}, {5, 1, 2, 1, 2, 1, 1}, {5, 1, 2, 1, 2, 1, 2}};
    vector<vi> basis5 = {{1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 2}, {1, 1, 1, 1, 1, 1, 2, 1, 1}, {1, 1, 1, 1, 1, 1, 2, 1, 2},
                        {1, 1, 1, 1, 2, 1, 1, 1, 1}, {1, 1, 1, 1, 2, 1, 1, 1, 2}, {1, 1, 1, 1, 2, 1, 2, 1, 1}, {1, 1, 1, 1, 2, 1, 2, 1, 2},
                        {1, 1, 2, 1, 1, 1, 1, 1, 1}, {1, 1, 2, 1, 1, 1, 1, 1, 2}, {1, 1, 2, 1, 1, 1, 2, 1, 1}, {1, 1, 2, 1, 1, 1, 2, 1, 2},
                        {1, 1, 2, 1, 2, 1, 1, 1, 1}, {1, 1, 2, 1, 2, 1, 1, 1, 2}, {1, 1, 2, 1, 2, 1, 2, 1, 1}, {1, 1, 2, 1, 2, 1, 2, 1, 2},
                        {2, 1, 1, 1, 1, 1, 1, 1, 1}, {2, 1, 1, 1, 1, 1, 1, 1, 2}, {2, 1, 1, 1, 1, 1, 2, 1, 1}, {2, 1, 1, 1, 1, 1, 2, 1, 2},
                        {2, 1, 1, 1, 2, 1, 1, 1, 1}, {2, 1, 1, 1, 2, 1, 1, 1, 2}, {2, 1, 1, 1, 2, 1, 2, 1, 1}, {2, 1, 1, 1, 2, 1, 2, 1, 2},
                        {2, 1, 2, 1, 1, 1, 1, 1, 1}, {2, 1, 2, 1, 1, 1, 1, 1, 2}, {2, 1, 2, 1, 1, 1, 2, 1, 1}, {2, 1, 2, 1, 1, 1, 2, 1, 2},
                        {2, 1, 2, 1, 2, 1, 1, 1, 1}, {2, 1, 2, 1, 2, 1, 1, 1, 2}, {2, 1, 2, 1, 2, 1, 2, 1, 1}, {2, 1, 2, 1, 2, 1, 2, 1, 2},
                        {5, 1, 1, 1, 1, 1, 1, 1, 1}, {5, 1, 1, 1, 1, 1, 1, 1, 2}, {5, 1, 1, 1, 1, 1, 2, 1, 1}, {5, 1, 1, 1, 1, 1, 2, 1, 2},
                        {5, 1, 1, 1, 2, 1, 1, 1, 1}, {5, 1, 1, 1, 2, 1, 1, 1, 2}, {5, 1, 1, 1, 2, 1, 2, 1, 1}, {5, 1, 1, 1, 2, 1, 2, 1, 2},
                        {5, 1, 2, 1, 1, 1, 1, 1, 1}, {5, 1, 2, 1, 1, 1, 1, 1, 2}, {5, 1, 2, 1, 1, 1, 2, 1, 1}, {5, 1, 2, 1, 1, 1, 2, 1, 2},
                        {5, 1, 2, 1, 2, 1, 1, 1, 1}, {5, 1, 2, 1, 2, 1, 1, 1, 2}, {5, 1, 2, 1, 2, 1, 2, 1, 1}, {5, 1, 2, 1, 2, 1, 2, 1, 2},};
    
    Matrix changeBasis1 = changeBasis(basis1, 0);
    Matrix changeBasis2 = changeBasis(basis2, 1);
    Matrix changeBasis3 = changeBasis(basis3, 2);
    Matrix changeBasis4 = changeBasis(basis4, 3);
    Matrix changeBasis5 = changeBasis(basis5, 4);
    auto inv = invMatBools(changeBasis1);
    assert(inv.first);
    Matrix invChangeBasis1 = inv.second.mat;
    inv = invMatBools(changeBasis2);
    assert(inv.first);
    Matrix invChangeBasis2 = inv.second.mat;
    inv = invMatBools(changeBasis3);
    assert(inv.first);
    Matrix invChangeBasis3 = inv.second.mat;
    inv = invMatBools(changeBasis4);
    assert(inv.first);
    Matrix invChangeBasis4 = inv.second;
    inv = invMatBools(changeBasis5);
    assert(inv.first);
    Matrix invChangeBasis5 = inv.second;

    /*
    // for computing the distance in the unknot code with one -, one + crossing
    
    Matrix leftBasisChange = kroneckerProduct(changeBasis2, changeBasis1);
    Matrix topLeftMat = kroneckerProduct(thetaSplit(1, 0), I3);
    Matrix botLeftMat = thetaMerge(1, 0);
    Matrix newTopLeftMat = matrixMult(kroneckerProduct(kroneckerProduct(invChangeBasis1, invChangeBasis1), invChangeBasis1), topLeftMat);
    Matrix newBotLeftMat = matrixMult(invChangeBasis3, botLeftMat);
    Matrix leftDiff = stackVertical(newTopLeftMat, newBotLeftMat);

    Matrix topRightBasisChange = kroneckerProduct(kroneckerProduct(changeBasis1, changeBasis1), changeBasis1);
    Matrix topRightMap = kroneckerProduct(I3, thetaMerge(0, 0));
    Matrix basisChangeBack = kroneckerProduct(invChangeBasis1, invChangeBasis2);
    Matrix newTopRightMap = matrixMult(basisChangeBack, matrixMult(topRightMap, topRightBasisChange));

    Matrix botRightBasisChange = changeBasis3;
    Matrix botRightMap = thetaSplit(2, 1);
    Matrix newBotRightMap = matrixMult(basisChangeBack, matrixMult(botRightMap, botRightBasisChange));
    Matrix rightDiff = stackHorizontal(newTopRightMap, newBotRightMap);
    vector<Matrix> maps = {leftDiff, rightDiff};

    // outputMatrix(matrixMult(rightDiff, leftDiff));
    */

    /*
    // code for 2 negative crossings 
    Matrix leftBasisChange = kroneckerProduct(kroneckerProduct(changeBasis1, changeBasis1), changeBasis1);
    Matrix topLeftMap = matrixMult(kroneckerProduct(I3, thetaMerge(0, 0)), leftBasisChange);
    Matrix basisChangeBackTop = kroneckerProduct(invChangeBasis1, invChangeBasis2);
    Matrix newTopLeftMap = matrixMult(basisChangeBackTop, topLeftMap);

    Matrix botLeftMap = matrixMult(kroneckerProduct(thetaMerge(0, 0), I3), leftBasisChange);
    Matrix basisChangeBackBot = kroneckerProduct(invChangeBasis2, invChangeBasis1);
    Matrix newBotLeftMap = matrixMult(basisChangeBackBot, botLeftMap);

    Matrix leftDiff = stackVertical(newTopLeftMap, newBotLeftMap);

    Matrix topRightMap = matrixMult(thetaMerge(0, 1), kroneckerProduct(changeBasis1, changeBasis2));
    Matrix botRightMap = matrixMult(thetaMerge(1, 0), kroneckerProduct(changeBasis2, changeBasis1));
    Matrix rightDiff = stackHorizontal(topRightMap, botRightMap);

    Matrix squared = matrixMult(rightDiff, leftDiff);
    */
    
    /*
    // code for 2 positive crossings, homology at rightmost degree
    Matrix topLeftMap = matrixMult(kroneckerProduct(invChangeBasis2, invChangeBasis1), thetaSplit(2, 0)).mat;
    Matrix botLeftMap = matrixMult(kroneckerProduct(invChangeBasis1, invChangeBasis2), thetaSplit(2, 1)).mat;
    Matrix leftDiff = stackVertical(topLeftMap, botLeftMap);

    Matrix topRightMap = matrixMult(kroneckerProduct(thetaSplit(1, 0), I3), kroneckerProduct(changeBasis2, changeBasis1)).mat;
    Matrix changeBasisFinal = kroneckerProduct(kroneckerProduct(invChangeBasis1, invChangeBasis1), invChangeBasis1);
    Matrix newTopRightMap = matrixMult(changeBasisFinal, topRightMap).mat;
    Matrix botRightMap = matrixMult(kroneckerProduct(I3, thetaSplit(1, 0)), kroneckerProduct(changeBasis1, changeBasis2)).mat;
    Matrix newBotRightMap = matrixMult(changeBasisFinal, botRightMap).mat;
    Matrix rightDiff = stackHorizontal(newTopRightMap, newBotRightMap);

    Matrix squared = matrixMult(rightDiff, leftDiff).mat;
    */
    

    /* to try to calculate unknot code, 4 crossings
    Matrix left11 = matrixMult(invChangeBasis5, matrixMult(thetaMerge(2, 1), kroneckerProduct(changeBasis3, changeBasis2)));
    Matrix left12 = matrixMult(kroneckerProduct(kroneckerProduct(invChangeBasis2, invChangeBasis1), invChangeBasis2),
                    kroneckerProduct(matrixMult(thetaSplit(2, 0), changeBasis3), changeBasis2));
    Matrix left13 = matrixMult(kroneckerProduct(kroneckerProduct(invChangeBasis1, invChangeBasis2), invChangeBasis2),
                    kroneckerProduct(matrixMult(thetaSplit(2, 1), changeBasis3), changeBasis2));

    Matrix left21 = matrixMult(invChangeBasis5, matrixMult(thetaMerge(3, 0), kroneckerProduct(changeBasis4, changeBasis1)));
    Matrix left22 = matrixMult(kroneckerProduct(kroneckerProduct(invChangeBasis2, invChangeBasis2), invChangeBasis1),
                    kroneckerProduct(matrixMult(thetaSplit(3, 1), changeBasis4), changeBasis1));
    Matrix left23 = kroneckerProduct(matrixMult(kroneckerProduct(invChangeBasis1, invChangeBasis3), matrixMult(thetaSplit(3, 2), changeBasis4)), I3);

    Matrix left31 = kroneckerProduct(I6, kroneckerProduct(I3, matrixMult(invChangeBasis2, matrixMult(thetaMerge(0, 0), kroneckerProduct(changeBasis1, changeBasis1)))));
    Matrix left32 = kroneckerProduct(I6, kroneckerProduct(matrixMult(invChangeBasis2, matrixMult(thetaMerge(0, 0), kroneckerProduct(changeBasis1, changeBasis1))), I3));
    Matrix left33 = kroneckerProduct(matrixMult(kroneckerProduct(invChangeBasis1, invChangeBasis1), matrixMult(thetaSplit(1, 0), changeBasis2)), kroneckerProduct(I3, kroneckerProduct(I3, I3)));

    Matrix left41 = kroneckerProduct(kroneckerProduct(I3, I6), matrixMult(invChangeBasis2, matrixMult(thetaMerge(0, 0), kroneckerProduct(changeBasis1, changeBasis1))));
    Matrix left42 = kroneckerProduct(I3, kroneckerProduct(matrixMult(invChangeBasis3, matrixMult(thetaMerge(1, 0), kroneckerProduct(changeBasis2, changeBasis1))), I3));
    Matrix left43 = kroneckerProduct(kroneckerProduct(I3, matrixMult(kroneckerProduct(invChangeBasis1, invChangeBasis1), matrixMult(thetaSplit(1, 0), changeBasis2))), kroneckerProduct(I3, I3));
    
    Matrix left1 = stackVertical({left11, left12, left13, Matrix(108+108+243, 72)});    
    Matrix left2 = stackVertical({left21, Matrix(108+108, 72), left22, left23, Matrix(243, 72)});
    Matrix left3 = stackVertical({Matrix(48, 162), left31, Matrix(108, 162), left32, Matrix(108, 162), left33});
    Matrix left4 = stackVertical({Matrix(48+108, 162), left41, Matrix(108, 162), left42, left43});

    Matrix leftDiff = stackHorizontal({left1, left2, left3, left4});


    Matrix right11 = matrixMult(kroneckerProduct(invChangeBasis2, invChangeBasis3), matrixMult(thetaSplit(4, 2), changeBasis5));
    Matrix right12 = matrixMult(kroneckerProduct(invChangeBasis1, invChangeBasis4), matrixMult(thetaSplit(4, 3), changeBasis5));

    Matrix right21 = kroneckerProduct(I6, matrixMult(invChangeBasis3, matrixMult(thetaMerge(0, 1), kroneckerProduct(changeBasis1, changeBasis2))));
    */


    // Matrix squared = matrixMult(rightDiff, leftDiff);
    // outputMatrix(squared);
    Matrix rightDiff;
    Matrix leftDiff;

    Matrix squared;
    forn(i, squared.size()){
        forn(j, squared[i].size()){
            assert(squared[i][j] == 0);
        }
    }

    return {leftDiff, rightDiff};
}