#include <bits/stdc++.h>
 
using namespace std;
 
const int N = 1e5 + 10, LOG_A = 20;

int basis[LOG_A];

int sz;

void insertVector(int mask) {
	for (int i = 0; i < LOG_A; i++) {
		if ((mask & 1 << i) == 0) continue;

		if (!basis[i]) {
			basis[i] = mask;
			++sz;
			
			return;
		}

		mask ^= basis[i];
	}
}

int main() {
	int n;
	cin >> n;

	for (int i = 0; i < n; i++) {
		int a;
		scanf("%d", &a);

		insertVector(a);
	}

	cout << "Rank: " << sz << endl;
    cout << "Nullity: " << n - sz << endl;

	return 0;
}