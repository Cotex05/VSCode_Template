#include <bits/stdc++.h>
#define ll unsigned long long
#define tc    \
    ll t;     \
    cin >> t; \
    while (t--)
#define all(x) x.begin(), x.end()
#define deb(x) cout << #x << " " << x << "\n";
const int INF = INT_MAX;
#define INFL LLONG_MAX
#define sum(a) (accumulate((a).begin(), (a).end(), 0ll))
#define pb push_back
#define eb emplace_back
#define F(left, a, n) for (int left = a; left < n; left++)
#define mod 1000000007
#define vi vector<int>
#define vll vector<long long int>
#define vvi vector<vector<int>>
#define vvl vector<vector<ll>>
#define vip vector<pair<int, int>>
#define NL endl
#define __ " "
#define fi first
#define se second
#define f(i, a, b) for (long long i = a; i < b; i++)
#define rf(i, a, b) for (long long i = a; i >= b; i--)
using namespace std;

clock_t time_p = clock();
void clk()
{
    time_p = clock() - time_p;
    cerr << "Execution Time: " << 100 * (float)(time_p) / CLOCKS_PER_SEC << "ms\n";
}

inline void multiply(ll mat1[2][2], ll mat2[2][2])
{
    ll new_mat[2][2] = {
        {(((mat1[0][0])%mod * (mat2[0][0])%mod)%mod + ((mat1[0][1])%mod * (mat2[1][0])%mod)%mod)%mod,
         (((mat1[0][0])%mod * (mat2[0][1])%mod)%mod + ((mat1[0][1])%mod * (mat2[1][1])%mod)%mod)%mod},
        {(((mat1[1][0])%mod * (mat2[0][0])%mod)%mod + ((mat1[1][1])%mod * (mat2[1][0])%mod)%mod)%mod,
         (((mat1[1][0])%mod * (mat2[0][1])%mod)%mod + ((mat1[1][1])%mod * (mat2[1][1])%mod)%mod)%mod}};
    mat1[0][0] = new_mat[0][0]%mod;
    mat1[0][1] = new_mat[0][1]%mod;
    mat1[1][0] = new_mat[1][0]%mod;
    mat1[1][1] = new_mat[1][1]%mod;
}
inline void power(ll mat[2][2], ll n)
{
    if (n == 0 || n == 1)
    {
        return;
    }
    power(mat, n / 2);
    multiply(mat, mat);
    if (n % 2 != 0)
    {
        ll temp[2][2] = {{1, 1}, {1, 0}};
        multiply(mat, temp);
    }
}
inline ll fib(ll n)
{
    ll mat[2][2] = {{1, 1}, {1, 0}};
    if (n == 0)
    {
        return 0;
    }
    power(mat, n - 1);
    return mat[0][0]%mod;
}
inline unsigned long long fiboSum(unsigned long long m,unsigned long long n)
{
    return (fib(n + 2)%mod - fib(m + 1)%mod+mod)%mod;
}


void faisal()
{
    ll m, n;
    cin >> m >> n;
    cout << fiboSum(m, n) << "\n";
}

int main()
{

#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    tc
    {
        faisal();
    }
    // clk();

    return 0;
}