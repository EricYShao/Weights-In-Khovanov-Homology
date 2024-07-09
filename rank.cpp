// This code only outputs the dimension of the image and
// the dimension of the kernel. To find an explicit basis,
// consider this source https://codeforces.com/blog/entry/98376

#include "bits/stdc++.h"
 
using namespace std;
using ll = long long;
 
const ll N = 1e5 + 10, LOG_A = 60;

ll basis[LOG_A];

ll sz;

void insertVector(ll mask) {
	for (ll i = 0; i < LOG_A; i++) {
		if ((mask & 1ll << i) == 0) continue;

		if (!basis[i]) {
			basis[i] = mask;
			++sz;
			
			return;
		}

		mask ^= basis[i];
	}
}

int main() {
	ll n;
	cin >> n;

	for (ll i = 0; i < n; i++) {
		ll a;
		scanf("%d", &a);

		insertVector(a);
	}

	cout << "Rank: " << sz << endl;
    cout << "Nullity: " << n - sz << endl;

	return 0;
}
