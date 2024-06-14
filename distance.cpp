#include <bits/stdc++.h>
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

    cout << "This program finds the minimum distance. Please enter the number of basis elements for the image: " << endl;
    ll n; cin >> n;
    while (n >= 60){
        cout << "Please enter a number less than 60" << endl;
        cin >> n;
    }
    cout << "Input " << n << " integers, each representing a basis element in (Z/2Z)^" << n << ':' << endl;
    
    vll vectors;

    forn(i, n){
        ll x; cin >> x;
        vectors.pb(x);
    }

    ll minDist = INF;

    for (ll i = 1; i < (1ll << n); i++){
        ll mask = 0;
        forn(j, n){
            if (i & 1ll << j){
                mask ^= vectors[j];
            }
        }
        if (!mask) minDist = min(minDist, (ll)__builtin_popcountll(i));
    }

    if (minDist != INF) cout << "The minimum distance is " << minDist << '.' << endl;
    else cout << "The kernel is empty." << endl;

    return 0;
}