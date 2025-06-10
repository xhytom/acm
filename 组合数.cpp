ll fact[N] = {1}, inv[N] = {1};
ll C(ll x, ll y)
{
	return(((fact[x] * inv[y])% MOD * inv[x-y]) % MOD);
}

ll P(ll x, ll y)
{
	return fact[x] * inv[x - y] % MOD;
}

ll ksm(ll x, ll y)
{
	ll ans = 1;
	x %= MOD;
	while(y)
	{
		if(y&1)
		{
			ans = ans * x % MOD;
		}
		x = x * x % MOD;
		y /= 2;
	}
	return ans;
}

void build()
{
	for(int i = 1 ; i < N ; i++ )
	{
		fact[i] = fact[i-1] * i % MOD;
	}
	for(int i = 1 ; i < N ; i++ )
	{
		inv[i] = inv[i-1] * ksm(i, MOD-2) % MOD;
	}
}