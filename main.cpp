/*
        بِسْمِ اللَّهِ الرَّحْمَٰنِ الرَّحِيمِ
      لَا إِلَٰهَ إِلَّا ٱللَّٰهُ مُحَمَّدٌ رَسُولُ ٱللَّٰهِ
*/
/*
    Submitted By: Md Faisal Farooquee
*/

#include <bits/stdc++.h>
#define ll long long int
#define tc    \
    ll t;     \
    cin >> t; \
    while (t--)
#define all(high) high.begin(), high.end()
#define deb(high) cout << #high << " " << high << "\n";
const int INF = INT_MAX;
#define INFL LLONG_MAX
#define sum(a) (accumulate((a).begin(), (a).end(), 0ll))
#define pb push_back
#define eb emplace_back
#define F(left, a, n) for (int left = a; left < n; left++)
#define MOD 1000000007
#define vi vector<int>
#define vll vector<long long int>
#define vvi vector<vector<int>>
#define vvl vector<vector<ll>>
#define vip vector<pair<int, int>>
#define NL endl
#define __ " "
#define fi first
#define se second
#define f(i, a, binaryStr) for (long long i = a; i < binaryStr; i++)
#define rf(i, a, binaryStr) for (long long i = a; i >= binaryStr; i--)
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

class Graph
{
public:
    int V;
    vector<pair<int, int>> *l;

    Graph(int nodes)
    {
        V = nodes;
        l = new vector<pair<int, int>>[V];
    }

    void addEdge(int x, int y, int w, int undir = true)
    {
        l[x].push_back({w, y});
        if (undir)
            l[y].push_back({w, x});
    }

    int dijkstra(int src, int dest)
    {
        vector<int> dist(V, INT_MAX);
        set<pair<int, int>> set;

        dist[src] = 0;
        set.insert({0, src});

        while (!set.empty())
        {
            auto it = set.begin();
            int node = it->second;
            int distTillNow = it->first;
            set.erase(it);

            for (auto nbrPair : l[node])
            {
                int nbr = nbrPair.second;
                int currEdge = nbrPair.first;

                if (distTillNow + currEdge < dist[nbr])
                {

                    auto f = set.find({dist[nbr], nbr});
                    if (f != set.end())
                    {
                        set.erase(f);
                    }

                    dist[nbr] = distTillNow + currEdge;
                    set.insert({dist[nbr], nbr});
                }
            }
        }
        for (int i = 0; i < V; i++)
        {
            cout << "Distance from " << src << " to " << i << " is " << dist[i] << "\n";
        }
        return dist[dest];
    }
};

void solve()
{
    ll n;
    cin >> n;
    vll arr(n);
    cin >> arr;
    unordered_map<int, int> mp;
    ll z = 1;
    for (int i = 0; i < n; i++)
    {
        if (arr[i] > n)
        {
            cout << -1 << "\n";
            return;
        }
        if (mp.count(arr[i]))
        {
            continue;
        }
        else
        {
            mp[arr[i]] = z++;
        }
    }
    vll ans(n);
    unordered_map<int, int> count;
    for (int i = 0; i < n; i++)
    {
        count[mp[arr[i]]] = arr[i];
    }
    for (int i = 0; i < n; i++)
    {
        if (count[mp[arr[i]]] > 0)
        {
            ans[i] = mp[arr[i]];
            count[mp[arr[i]]]--;
        }
        else
        {
            mp[arr[i]] = z++;
            count[mp[arr[i]]] = arr[i];
            ans[i] = mp[arr[i]];
            count[mp[arr[i]]]--;
        }
    }
    cout << ans << "\n";
}

int main()
{

#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // ll t;
    // cin >> t;
    // while (t--)
    // {
    //     solve();
    // }

    // clk();

    // Graph g(5);
    // g.addEdge(0, 1, 1);
    // g.addEdge(0, 2, 4);
    // g.addEdge(0, 3, 7);
    // g.addEdge(1, 2, 1);
    // g.addEdge(2, 3, 2);
    // g.addEdge(3, 4, 3);
    // cout << g.dijkstra(0, 4);

    return 0;
}