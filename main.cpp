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

const int borne = (int)(1e6) + 5;
const int WAIT = 0, ENTERED = 1, LEFT = 2;
int state[borne];

long long Ceiling(long long a, long long b)
{
    return (a + b - 1) / b;
}

clock_t time_p = clock();
void clk()
{
    time_p = clock() - time_p;
    cerr << "Execution Time: " << 100 * (float)(time_p) / CLOCKS_PER_SEC << "ms\n";
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

ll nearestSquare(ll n)
{
    return power(ceil(sqrt(n)), 2);
}

void faisal()
{
    ll n;
    cin >> n;
    vector<ll> arr(n, -1);
    ll val = n - 1;
    ll c = 0;
    while (c < n)
    {
        ll ns = nearestSquare(val);
        ll k = ns - val;
        for (ll i = k; i < n; i++)
        {
            if (arr[i] == -1)
            {
                arr[i] = val;
                val--;
                c++;
            }
            else
            {
                break;
            }
        }
    }
    for (auto e : arr)
    {
        cout << e << " ";
    }
    cout << "\n";
}

bool isPrime(int x)
{
    x = abs(x);
    if(x<2){
        return false;
    }
    for (int i = 2; i*i <= x; i++)
    {
        if (x % i == 0)
        {
            return false;
        }
    }
    return true;
}

int sumPrime(int l, int r){
    int ans = 0;
    for(int i=l; i<=r; i++){
        if(isPrime(i)){
            ans = i;
            break;
        }
    }
    for(int i=r; i>=l; i--){
        if(isPrime(i)){
            ans+=i;
            break;
        }
    }
    return ans;
}

vector<int> funcArrange(vector<int> inputArr){
    int n = inputArr.size();
    vector<int> ans;
    for(int i=0; i<n; i++){
        if(inputArr[i]%2==0){
            ans.push_back(inputArr[i]);
        }
    }
    for(int i=0; i<n; i++){
        if(inputArr[i]%2==1){
            ans.push_back(inputArr[i]);
        }
    }
    return ans;
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
        int n;
        cin >> n;
        vector<int> arr(n);
        for(auto &e: arr){
            cin >> e;
        }
        auto v = funcArrange(arr);
        for(auto e: v){
            cout << e << " ";
        }
    }
    // clk();

    return 0;
}