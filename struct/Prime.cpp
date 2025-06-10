struct Prime{
    vector<int> prime, is;
    int n;
    Prime(int _n) : is(_n + 1, 1) 
    {
        n = _n;
        for(int i = 2 ; i <= n ; i++ )
        {
            if(is[i] == 1) prime.push_back(i);
            for(int j = 0 ; i * prime[j] <= n && j < prime.size() ; j++ )
            {
                is[i * prime[j]] = 0;
                if(i % prime[j] == 0) break;
            }
        }
    }
     
};