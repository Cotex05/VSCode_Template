#include <bits/stdc++.h>
#define ll long long int
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

template <typename T>
std::ostream &operator<<(std::ostream &output, std::vector<T> const &values)
{
    for (auto const &value : values)
    {
        output << value << " ";
    }
    return output;
}

template <typename T>
std::istream &operator>>(std::istream &input, std::vector<T> &values)
{
    for (auto &value : values)
    {
        input >> value;
    }
    return input;
}

ll power(ll x, ll y)
{
    ll result = 1;
    while (y > 0)
    {
        if (y % 2 == 0) // y is even
        {
            x = x * x;
            y = y / 2;
        }
        else // y isn't even
        {
            result = result * x;
            y = y - 1;
        }
    }
    return result;
}

char x, y;

struct Comp
{
    bool operator()(const pair<ll, ll> &a, const pair<ll, ll> &b)
    {
        return (abs(a.first - x) > abs(b.first - x));
    }
};

void solveC()
{
    ll n;
    cin >> n;
    vll arr(n), brr(n);
    cin >> arr >> brr;
    
}

int main()
{

#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // tc
    // {
    //     solveC();
    // }
    solveC();
    // clk();
    return 0;
}