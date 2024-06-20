#include <iostream>
#include <vector>
#include <map>

using namespace std;

#define sz(x) ((int)x.size())
#define all(a) (a).begin(), (a).end()
#define pb push_back

using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;

#define forn(i, n) for (int i = 0; i < int(n); i++)

const ll INF = 1e18;

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // cout << "This program finds the minimum distance. Please enter the number of basis elements for the image: " << endl;
    // ll n; cin >> n;
    // while (n >= 120){
    //     cout << "Please enter a number less than 120" << endl;
    //     cin >> n;
    // }
    // cout << "Input " << n << " integers, each representing a basis element in (Z/2Z)^" << n << ':' << endl;
    
    vll vectors;

    // forn(i, n){
    //     ll x; cin >> x;
    //     vectors.pb(x);
    // }
    ll x;
    while (cin >> x) vectors.pb(x);
    ll n = vectors.size();
    cout << n << " column vectors" << endl;

    ll minDist = INF;

    map<ll, ll> maskCount;

    for (ll i = 1; i < (1ll << (n/2)); i++){
        ll mask = 0;
        if ((ll)__builtin_popcountll(i) >= minDist) continue;
        forn(j, n/2){
            if (i & (1ll << j)){
                mask ^= vectors[j];
            }
        }
        if (!maskCount.count(mask)) maskCount[mask] = (ll)__builtin_popcountll(i);
        else maskCount[mask] = min(maskCount[mask], (ll)__builtin_popcountll(i));
        if (!mask) minDist = min(minDist, (ll)__builtin_popcountll(i));
    }

    for (ll i = 1; i < (1ll << ((n+1)/2)); i++){
        if ((ll)__builtin_popcountll(i) >= minDist) continue;
        ll mask = 0;
        forn(j, (n+1)/2){
            if (i & (1ll << j)){
                mask ^= vectors[j + n/2];
            }
        }
        if (maskCount.count(mask)){
            minDist = min(minDist, __builtin_popcountll(i) + maskCount[mask]);
        }
        if (!mask) minDist = min(minDist, (ll)__builtin_popcountll(i));
    }

    if (minDist != INF) cout << "The minimum distance is " << minDist << '.' << endl;
    else cout << "The kernel is empty." << endl;

    return 0;
}