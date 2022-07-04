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
using namespace std;

void solve1()
{
    ll n;
    cin >> n;
    vector<pair<int, int>> vp(n);
    for (int i = 0; i < n; i++)
    {
        cin >> vp[i].first >> vp[i].second;
    }
    vector<int> hash(n + 1);
    for (auto p : vp)
    {
        for (int i = p.first; i <= p.second; i++)
        {
            hash[i]++;
        }
    }
    int mx = 0;
    for (int i = 1; i <= n; i++)
    {
        mx = max(hash[i], mx);
    }
    int c = 0;
    for (int i = 1; i < n; i++)
    {
        if (hash[i] == mx)
        {
            c++;
        }
    }
    cout << c << "\n";
    for (int i = 1; i <= n; i++)
    {
        if (hash[i] == mx)
        {
            cout << i << "\n";
        }
    }
}

int solve(int n, string s)
{
    int result = 0;
 
    int count[26] = {0};
    for (int i=0; i<n; i++)
        count[s[i]-'a']++;
 
    for (int i=0; i<26; i++)
        result += (count[i]*(count[i]+1)/2);
 
    int ans = n*(n+1)/2 - result;
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
        string s;
        cin >> s;
        cout << solve(n, s);
    }

    return 0;
}