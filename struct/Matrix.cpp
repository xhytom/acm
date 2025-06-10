struct Matrix{
    int n, m ;
    vector<vector<Z>> s;
     
    Matrix(int _n, int _m) : n(_n), m(_m), s(_n, vector<Z>(_m, 0)){}
     
    friend Matrix operator * (Matrix a, Matrix b){
        assert(a.m == b.n);    
        Matrix res(a.n , b.m);
        for(int k = 0; k < a.m; k++)
            for(int i = 0; i < a.n; i++)
                for(int j = 0; j < b.m; j++)
                    res.s[i][j] = (res.s[i][j] + a.s[i][k] * b.s[k][j]);
        return res;
    }
     
    Matrix qmi(int b){
        assert(n == m);
        Matrix res(n, n);
        for(int i = 0; i < n; i ++ )
            res.s[i][i] = 1;
        while(b){
            if(b & 1) res = ((*this) * res );
            b >>= 1;
            *this = (*this) * (*this);
        }
        return (*this) = res;
    };
     
    Matrix inv() {
        std::vector<std::vector<Z>> a(n, std::vector<Z> (2 * n));
         
        for (int i = 0; i < n; i++) {
            a[i][i + n] = 1;
            for (int j = 0; j < n; j++) {
                a[i][j] = s[i][j];
            }
        }
         
        for (int i = 0; i < n; i++) {
            Z tmp = Z(1) / a[i][i];
            for (int j = i; j < 2 * n; j++) {
                a[i][j] *= tmp;
            }
             
            for (int k = 0; k < n; k++) {
                Z tmp = a[k][i];
                if (k == i) continue;
                for (int j = i; j < 2 * n; j++) {
                    a[k][j] -= tmp * a[i][j];
                }
            }
        }
         
        Matrix ans(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                ans.s[i][j] = a[i][j + n];
            }
        }
        return ans;
    }
     
};