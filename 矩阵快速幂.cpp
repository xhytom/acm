template<typename T>
struct Matrix{
    int n , m ;
    vector<vector<T>> s;
    
    Matrix(int n ,int m):n(n) ,m(m), s(n, vector<T>(m, T())){}
    
    friend Matrix operator * (Matrix a , Matrix b){
        assert(a.m == b.n);    
        Matrix res(a.n , b.m);
        for(int k = 0 ; k < a.m ; k ++ )
            for(int i = 0 ; i < a.n ; i ++ )
                for(int j = 0 ; j < b.m ; j ++ )
                    res.s[i][j] = res.s[i][j] + a.s[i][k] * b.s[k][j];
        return res;
    }
    
    Matrix qmi(T b){
        assert(n == m);
        Matrix res(n , n);
        for(int i = 0 ; i < n ; i ++ )
            res.s[i][i] = 1;
        while(b){
            if(b & 1)res = ((*this) * res );
            b >>= 1;
            *this = (*this) * (*this);
        }
        return (*this) = res;
    };
    
};

