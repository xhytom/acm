
std::mt19937 rng(std::random_device{}());
int A[]={2,325,9375,28178,450775,9780504,1795265022};

int qpow(int a,int b,int p){
    int ans=1;
    for(;b;b>>=1){
        if(b&1) ans=(__int128)ans*a%p;
        a=(__int128)a*a%p;
    }
    return ans;
}

bool isprime(int x){
    if(x<3) return x==2;
    if(x%2==0) return false;
    int d=x-1,r=__builtin_ctz(d);
    d>>=r;
    for(auto a:A){
        int v=qpow(a,d,x);
        if(v<=1 || v==x-1) continue;
        for(int i=0;i<r;++i){
            v=(__int128)v*v%x;
            if(v==x-1 && i!=r-1){
                v=1;
                break;
            }
            if(v==1) return false;
        }
        if(v!=1) return false;
    }
    return true;
}

int Pollard_Rho(int x) {
    int s=0,t=0,val=1;
    int c=rng()%(x-1)+1;
    for(int i=1;;i*=2,s=t,val=1){
        for(int j=1;j<=i;++j){
            t=((__int128)t*t+c)%x;
            val=(__int128)val*abs(t-s)%x;
            if(j%127==0){
                int d=__gcd(val,x);
                if(d>1) return d;
            }
        }
        int d=__gcd(val,x);
        if(d>1) return d;
    }
}

void getfac(int x, map<int, int> &c) {
    if(x<2) return;
    if(isprime(x)){
        c[x]++;
        return;
    }
    int p=x;
    while(p>=x) p=Pollard_Rho(x);
    if(x%p==0) x/=p;
    getfac(x, c),getfac(p, c);
}
map<int, int> getprimes(int x) {
    map<int, int> c;
    getfac(x, c);
    return c;
}