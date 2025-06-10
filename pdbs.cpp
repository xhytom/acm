#include<ext/pb_ds/tree_policy.hpp>
#include<ext/pb_ds/assoc_container.hpp>

using namespace __gnu_pbds;
__gnu_pbds::tree<ll, null_type, less<ll>, rb_tree_tag, tree_order_statistics_node_update> T;

if(op == 1)
{
    T.insert({x, i});
}else if (op == 2)
{
    T.erase(T.lower_bound({x, 0}));
}else if (op == 3)
{
    cout << T.order_of_key({x, 0}) + 1 << "\n";
}else if (op == 4)
{
    cout << T.find_by_order(x - 1)->first << "\n";
}else if (op == 5)
{
    cout << prev(T.lower_bound({x, 0}))->first << "\n";
}else if (op == 6)
{
    cout << T.lower_bound({x + 1, 0})->first << "\n";
}