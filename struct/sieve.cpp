int tot;
int p[N], pr[N], pe[N];
// p[x]为x最小的质因子, pr[x]为第x个质数, pe[x]为x的最小质因子个数幂

void sieve(int n) {
	for (int i = 2; i <= n; i++) {
		if (pe[i] == 0) p[i] = i, pe[i] = i, pr[++tot] = i; 
		for (int j = 1; j <= tot && i * pr[j] <= n; j++) {
			p[i * pr[j]] = pr[j];
			if (p[i] == pr[j]) {
				pe[i * pr[j]] = pe[i] * pr[j];
				break;
			} else {
				pe[i * pr[j]] = pr[j];
			}
		}
	}
}

void compute(int f[], int n, std::function<int(int)> calcpe) {
	f[1] = 1;
	for (int i = 2; i <= n; i++) {
		if (i == pe[i]) f[i] = calcpe(i);
		else f[i] = f[pe[i]] * f[i / pe[i]];
	}
}