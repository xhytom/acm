void fft(vector<complex<double>>&a){
    int n=a.size(),L=31-__builtin_clz(n);
    vector<complex<long double>>R(2,1);
    vector<complex<double>>rt(2,1);
    for(int k=2;k<n;k*=2){
        R.resize(n);
        rt.resize(n);
        auto x=polar(1.0L, acos(-1.0L)/k);
        for(int i=k;i<2*k;++i) rt[i]=R[i]=i&1?R[i/2]*x:R[i/2];
    }
    vector<int>rev(n);
    for(int i=0;i<n;++i) rev[i]=(rev[i/2]|(i&1)<<L)/2;
    for(int i=0;i<n;++i) if(i<rev[i]) swap(a[i],a[rev[i]]);
    for (int k=1;k<n;k*=2)
        for (int i=0;i<n;i+=2*k)
            for(int j=0;j<k;++j){
                complex<double>z=rt[j+k]*a[i+j+k];
                a[i+j+k]=a[i+j]-z;
                a[i+j]+=z;
            }
}

vector<double>mul(const vector<double>&a,const vector<double>&b){
    if(a.empty() || b.empty()) return {};
    vector<double>res(a.size()+b.size()-1);
    int L=32-__builtin_clz(res.size()),n=1<<L;
    vector<complex<double>>in(n),out(n);
    copy(a.begin(),a.end(),in.begin());
    for(int i=0;i<b.size();++i) in[i].imag(b[i]);
    fft(in);
    for(auto &x:in) x*=x;
    for(int i=0;i<n;++i) out[i]=in[-i&(n-1)]-conj(in[i]);
    fft(out);
    for(int i=0;i<res.size();++i) res[i]=imag(out[i])/(4 * n);
    return res;
}