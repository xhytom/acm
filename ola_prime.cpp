#include<bits/stdc++.h>
using namespace std;
const int PN = 1000; 
int isprime[PN+5]={1,1,0,0},idx_prime=1;
int prime[PN]={0};
void build()
{
	for(int i = 2 ; i <= PN ; i++ )
	{
		if(isprime[i]==0) prime[idx_prime++] = i; 
		for(int j = 1 ; j <= idx_prime &&i*prime[j]<=PN; j++ )
		{
			isprime[i*prime[j]] = 1;
			if(i%prime[j]==0) break;
		}
	}
	// for(int i = 1 ; i <= idx_prime ; i++ )
	// {
	// 	cout << prime[i] << "\n";
	// }
}



int main()
{
	build();

}