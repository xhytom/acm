int mod=998244353;
 
int qpow(int a,int b){
    int ans=1;
    for(;b;b>>=1){
        if(b&1) ans=1LL*ans*a%mod;
        a=1LL*a*a%mod;
    }
    return ans;
}
 
vector<int>roots{0,1};
vector<int>rev;
 
void dft(vector<int>&a){
    int n=a.size();
    if(rev.size()!=n){
        rev.resize(n);
        int k=__builtin_ctz(n)-1;
        for(int i=0;i<n;++i) rev[i]=rev[i>>1]>>1|(i&1)<<k;
    }
    for(int i=0;i<n;++i) if(i<rev[i]) swap(a[i],a[rev[i]]);
    if(roots.size()<n){
        int k=__builtin_ctz(roots.size());
        roots.resize(n);
        while((1<<k)<n){
            int e=qpow(3,(mod-1)>>(k+1));
            for(int i=(1<<(k-1));i<(1<<k);++i){
                roots[2*i]=roots[i];
                roots[2*i+1]=1LL*roots[i]*e%mod;
            }
            k++;
        }
    }
    for(int k=1;k<n;k<<=1){
        for(int i=0;i<n;i+=2*k){
            for(int j=0;j<k;++j){
                int u=a[i+j],v=1LL*a[i+j+k]*roots[k+j]%mod;
                a[i+j]=(u+v)%mod;
                a[i+j+k]=(u-v+mod)%mod;
            }
        }
    }
}
void idft(vector<int>&a){
    reverse(a.begin()+1,a.end());
    dft(a);
    int n=a.size(),inv=(1-mod)/n+mod;
    for(int i=0;i<n;++i) a[i]=1LL*a[i]*inv%mod;
}
 
struct Poly{
    vector<int>a;
    friend Poly operator*(Poly a,Poly b){
        int sz=1,tot=a.a.size()+b.a.size()-1;
        while(sz<tot) sz<<=1;
        a.a.resize(sz);b.a.resize(sz);
        dft(a.a);dft(b.a);
        for(int i=0;i<sz;++i) a.a[i]=1LL*a.a[i]*b.a[i]%mod;
        idft(a.a);
        a.a.resize(tot);
        return a;
    }
};
