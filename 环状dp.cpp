#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
const ll MOD_N = 998244353;
int dp1[405][405],dp2[405][405];
int num[405],sum[405];
int main()
{
	memset(dp2,0x3f,sizeof dp2);
	int n;
	cin >> n;
	rep(i,1,n+1)
	{
		cin >> num[i];
		num[i+n] = num[i];
	}
	rep(i,1,n+n+1)
	{
		sum[i] = sum[i-1] + num[i];
		dp1[i][i] = 0;
		dp2[i][i] = 0;
	}
	int maxans = 0,minans = INT_MAX;
	rep(p,1,n+1)
	{
		rep(len,2,n+1) 
		{

				for(int l = p ; l+len-1 <= p+n-1 ; l++)
				{
					int r = l+len-1;
					rep(k,l,r)
					{
						dp1[l][r] = max(dp1[l][r],dp1[l][k]+dp1[k+1][r]+sum[r]-sum[l-1]);
						dp2[l][r] = min(dp2[l][r],dp2[l][k]+dp2[k+1][r]+sum[r]-sum[l-1]);
					//	printf("dp1[%d][%d] = %d\n",l,r,dp1[l][r]);
					//	printf("dp2[%d][%d] = %d\n",l,r,dp2[l][r]);
					}
				}
		}
		maxans = max(maxans,dp1[p][p+n-1]);
		minans = min(minans,dp2[p][p+n-1]);
		memset(dp2,0x3f,sizeof dp2);
		memset(dp1,0,sizeof dp1);
		rep(i,1,2*n+1)
		{
			dp1[i][i] = 0;
			dp2[i][i] = 0;
		}
	}
	cout << minans << endl << maxans;
	return 0;
}
