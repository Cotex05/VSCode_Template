
void BFS(ll s, vector<int> adj[], ll n)
{
    bool *visited = new bool[n];
    for (ll i = 0; i < n; i++)
        visited[i] = false;

    list<ll> queue;
    visited[s] = true;
    queue.push_back(s);

    while (!queue.empty())
    {
        s = queue.front();
        cout << s << " ";
        queue.pop_front();

        for (auto i : adj[s])
        {
            if (!visited[i])
            {
                visited[i] = true;
                queue.push_back(i);
            }
        }
    }
}

#include <bits/stdc++.h>
using namespace std;
#define ll long long

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    ll kitne_cases_hain;
    kitne_cases_hain = 1;
    cin >> kitne_cases_hain;
    while (kitne_cases_hain--)
    {
        ll n;
        cin >> n;
        string s;
        cin >> s;
        if (s[0] == '0')
        {
            for (int i = 0; i < 2 * n; i++)
            {
                if (s[i] == '1')
                {
                    s[i] = '0';
                }
                else
                {
                    s[i] = '1';
                }
            }
        }
        vector<ll> v;
        ll cnt;
        ll flag = 0;
        ll p;
        for (int i = 0; i < 2 * n; i++)
        {
            if (s[i] == '0')
            {
                v.push_back(i);
            }
            if (s[i] != s[2 * n - i - 1])
            {
                flag = 1;
                break;
            }
        }
        if (flag)
        {
            cout << "1\n";
            cout << 2 * n << "\n";
            continue;
        }
        if (v.size() == 0)
        {
            cout << "-1\n";
        }
        else
        {
            cnt = 2 * n - 1 - v.back();
            for (int i = v.size() - 1; i >= 1; i--)
            {
                if (v[i] - v[i - 1] - 1 != cnt)
                {
                    flag = 1;
                    cout << "2\n";
                    cout << v[i - 1] + 1 << " " << 2 * n - 1 - v[i - 1] << "\n";
                    break;
                }
            }
            if (flag == 0)
            {
                cout << "2\n";
                cout << v[v.size() - 2] + 2 << " " << 2 * n - 2 - v[v.size() - 2] << "\n";
            }
        }
    }
    return 0;
}

if (n == 1)
    {
        cout << arr[0] << "\n";
        return;
    }
    unordered_map<ll, ll> mp;
    ll mx = 0;
    for (auto &e : arr)
    {
        mp[e]++;
        mx = max(e, mx);
    }
    ll ans = 0;
    for (ll i = mx; i >= 1; i--)
    {
        ll c = i;
        auto hash = mp;
        for (ll j = i; j > 0;)
        {
            ll f = 0;
            if (hash.find(j) != hash.end())
            {
                hash[j]--;
                if (hash[j] <= 0)
                {
                    hash.erase(j);
                }
                ll k = 1;
                while (1)
                {
                    if (c - 2 * k <= 0)
                    {
                        f = 1;
                        break;
                    }
                    if (hash.find(c - 2 * k) != hash.end())
                    {
                        hash[c - 2 * k]--;
                        if (hash[c - 2 * k] <= 0)
                        {
                            hash.erase(c - 2 * k);
                        }
                        hash[c - 2 * k + j - 1]++;
                        break;
                    }
                    k++;
                }
                j -= 2;
            }
            else
            {
                c = 0;
                break;
            }
            if (f)
            {
                c = 0;
                break;
            }
        }
        ans = max(c, ans);
    }
    cout << ans << "\n";

#include <bits/stdc++.h>
#define int long long
using namespace std;
int n, m, s, t;
string a[500001];
map<pair<int, int>, bool> vis;
void dfs(int x, int y)
{
    if (x < 1 || y < 1 || x > n || y > m)
        return;
    if (vis[{x, y}] || a[x][y - 1] == '#')
        return;
    t = max(t, x);
    vis[{x, y}] = true;
    dfs(x - 1, y);
    dfs(x + 1, y);
    dfs(x, y - 1);
    dfs(x, y + 1);
}
void solve()
{
    vis.clear();
    t = s = 0;
    cin >> n >> m;
    for (int i = 1; i <= n; ++i)
        cin >> a[i];
    if (n == 1)
    {
        for (int i = 0; i < m; ++i)
            if (a[1][i] == '#')
            {
                puts("1");
                return;
            }
        puts("0");
        return;
    }
    dfs(1, 1);
    while (!vis[{n, m}])
    {
        s++;
        if (t == n)
            break;
        for (int i = 1; i <= m; ++i)
            a[t + 1][i - 1] = '.';
        for (int i = 1; i <= m; ++i)
            vis[{t + 1, i}] = false;
        dfs(t + 1, 1);
    }
    cout << s << endl;
}
signed main()
{
    int T;
    cin >> T;
    while (T--)
        solve();
}

void DFS(int s, vector<int> adj[], int n, vector<bool> &vis)
{
    vis[s] = true;
    cout << s << " ";
    for (auto it : adj[s])
        if (!vis[it])
            DFS(it, adj, n, vis);
}

class Graph
{

    list<pair<int, int>> edge_list;
    int V;

public:
    Graph(int V)
    {
        this->V = V;
    }

    void addEdge(int u, int v)
    {
        edge_list.push_back(make_pair(u, v));
    }

    auto getEdgeList()
    {
        return edge_list;
    }

    int getSize()
    {
        return V;
    }
};

int maxWeightNode(vector<int> &arr)
{
    unordered_map<int, int> mp;
    int n = arr.size();
    for (int i = 0; i < n; i++)
    {
        if (arr[i] != -1)
            mp[arr[i]] += i;
    }
    int mx = INT_MIN, node = -1;
    for (int i = 0; i < n; i++)
    {
        if (mp[i] > mx)
        {
            mx = mp[i];
            node = i;
        }
    }
    cout << node << " ";
}

int maxSumCycle(vector<int> &arr)
{
    vector<int> sm;
    int n = arr.size();
    for (int i = 0; i < n; i++)
    {
        vector<int> path;
        int j = i;
        int temp = 0;
        while (arr[j] < n && arr[j] != i && arr[j] != -1 && find(path.begin(), path.end(), j) == path.end())
        {
            path.push_back(j);
            temp += j;
            j = arr[j];
            if (arr[j] == i)
            {
                temp += j;
                break;
            }
        }
        if (j < n && i == arr[j])
            sm.push_back(temp);
    }
    if (sm.size() < 1)
        return -1;
    return *max_element(sm.begin(), sm.end());
}

// largest cycle
// res stores result
int res = 0;
// visit to check in before visiting the node, to stop repeat visiting
unordered_map<int, bool> visit;

void dfs(vector<int> &a, unordered_map<int, int> &mp, int i, int k)
{
    if (visit.find(i) != visit.end()) // already visited
        return;
    if (a[i] == -1)
    {
        visit[i] = true;
        return;
    }
    if (mp.find(i) != mp.end())
    {
        res = max(res, k - mp[i]);
        visit[i] = true;
        return;
    }
    mp[i] = k;
    dfs(a, mp, a[i], k + 1);
    visit[i] = true;
}

// then find for each vertex
for (int i = 0; i < n; i++)
{
    if (visit.find(i) == visit.end())
    {
        unordered_map<int, int> mp;
        dfs(arr, mp, i, 0);
    }
}
cout << res << endl;

int main()
{

#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(0, 4);
    g.addEdge(4, 5);
    g.addEdge(3, 5);
    g.addEdge(3, 2);
    g.addEdge(2, 1);
    g.addEdge(3, 4);

    int n = g.getSize();

    auto edgeList = g.getEdgeList();

    vector<int> adj[n];
    for (auto e : edgeList)
    {
        adj[e.first].push_back(e.second);
    }

    int mx = INT_MIN;

    unordered_map<int, vector<int>> mp;

    for (auto e : edgeList)
    {
        mp[e.second].push_back(e.first);
    }

    for (auto it : mp)
    {
        int sm = 0;
        for (auto e : it.second)
        {
            sm += e;
        }
        mx = max(mx, sm);
    }

    cout << mx << "\n";

    return 0;
}

// tree of space
#include <bits/stdc++.h>
using namespace std;

unordered_map<string, int> mpNodes; // {name, position}

vector<pair<pair<string, int>, string>> findInitial(vector<vector<string>> &queries, unordered_map<string, pair<string, string>> &status)
{
    vector<pair<pair<string, int>, string>> initial;
    for (auto e : queries)
    {
        stringstream ss;
        ss << e[0];
        int id;
        ss >> id;
        initial.push_back({{e[1], id}, e[2]});
    }
    for (auto val : initial)
    {
        status[val.first.first] = {"0", "-1"};
    }
    return initial;
}

// Lock implementation
string lock(string name, vector<string> &nodes, unordered_map<string, pair<string, string>> &status, string lockedBy)
{
    // auto pos = find(nodes.begin(), nodes.end(), name);
    // int index = pos - nodes.begin();
    int index = mpNodes[name];
    index++;
    int n = nodes.size();
    int child1 = index * 2;
    int child2 = index * 2 + 1;
    if ((child1 > 0 && child1 < n) && (child2 > 0 && child2 < n))
    {
        status[nodes[child1 - 1]].first = "fail";
        status[nodes[child2 - 1]].first = "fail";
    }
    if (status[name].first == "lock" || status[name].first == "fail")
    {
        return "false";
    }
    else
    {
        int half = index / 2;
        status[nodes[half - 1]].first = "fail";
        status[name].first = "lock";
        status[name].second = lockedBy;
        return "true";
    }
}

// Unlock implementation
string unlock(string name, vector<string> &nodes, unordered_map<string, pair<string, string>> &status, string lockedBy)
{
    int index = mpNodes[name];
    index++;
    if (status[name].first == "lock" && status[name].second == lockedBy)
    {
        status[name].first = "unlock";
        status[name].second = "-1";
        int n = nodes.size();
        int child1 = index * 2;
        int child2 = index * 2 + 1;
        if ((child1 > 0 && child1 < n) && (child2 > 0 && child2 < n))
        {
            status[nodes[child1 - 1]].first = "unlock";
            status[nodes[child1 - 1]].second = "-1";
            status[nodes[child2 - 1]].first = "unlock";
            status[nodes[child2 - 1]].second = "-1";
            status[nodes[index - 1]].first = "unlock";
            status[nodes[index - 1]].second = "-1";
        }
        return "true";
    }
    else
    {
        return "false";
    }
}

// upgradeLock implementation
string upgradeLock(string name, vector<string> &nodes, unordered_map<string, pair<string, string>> &status, string lockedBy)
{
    // auto pos = find(nodes.begin(), nodes.end(), name);
    // int index = pos - nodes.begin();
    int index = mpNodes[name];
    index++;
    int child1 = index * 2;
    int child2 = index * 2 + 1;
    int n = nodes.size();
    if ((child1 > 0 && child1 < n) && (child2 > 0 && child2 < n))
    { // boundary check
        if ((status[nodes[child1 - 1]].first == "lock" && status[nodes[child1 - 1]].second == lockedBy) && (status[nodes[child2 - 1]].first == "lock" && status[nodes[child2 - 1]].second == lockedBy))
        {
            status[nodes[child1 - 1]].first = "unlock";
            status[nodes[child1 - 1]].second = "-1";
            status[nodes[child2 - 1]].first = "unlock";
            status[nodes[child2 - 1]].second = "-1";
            status[nodes[index - 1]].first = "lock";
            status[nodes[index - 1]].second = lockedBy;
            return "true";
        }
        else
        {
            status[nodes[index - 1]].first = "fail";
            return "false";
        }
    }
}

string operationType(string name, int code, vector<string> &nodes, unordered_map<string, pair<string, string>> &status, string lockedBy)
{
    string result = "false";
    if (code == 1)
    {
        result = lock(name, nodes, status, lockedBy);
    }
    else if (code == 2)
    {
        result = unlock(name, nodes, status, lockedBy);
    }
    else if (code == 3)
    {
        result = upgradeLock(name, nodes, status, lockedBy);
    }
    return result;
}

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m, q;

    cin >> n >> m >> q;

    unordered_map<string, pair<string, string>> status; // {name, status}
    vector<string> nodes(n);
    vector<vector<string>> queries(q, vector<string>(3)); // {type, name, uid}

    for (int i = 0; i < n; i++)
    {
        cin >> nodes[i];
        mpNodes[nodes[i]] = i;
    }
    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cin >> queries[i][j];
        }
    }

    vector<pair<pair<string, int>, string>> initial = findInitial(queries, status);

    for (int i = 0; i < initial.size(); i++)
    {
        stringstream ss;
        ss << initial[i].first.second;
        int id;
        ss >> id;
        string lockedBy = initial[i].second;
        string ans = operationType(initial[i].first.first, id, nodes, status, lockedBy);
        cout << ans << "\n";
    }
}

// Codevita
int check(char c, vector<vector<char>> &board)
{
    int count = 0;
    // 1 row
    if (board[0][0] == c)
    {
        if (board[0][1] == c)
        {
            if (board[0][2] == c)
            {
                count++;
            }
        }
    }
    // 2 row
    if (board[1][0] == c)
    {
        if (board[1][1] == c)
        {
            if (board[1][2] == c)
            {
                count++;
            }
        }
    }
    // 3 row
    if (board[2][0] == c)
    {
        if (board[2][1] == c)
        {
            if (board[2][2] == c)
            {
                count++;
            }
        }
    }
    // 1 col
    if (board[0][0] == c)
    {
        if (board[1][0] == c)
        {
            if (board[2][0] == c)
            {
                count++;
            }
        }
    }
    // 2 col
    if (board[0][1] == c)
    {
        if (board[1][1] == c)
        {
            if (board[2][1] == c)
            {
                count++;
            }
        }
    }
    // 3 col
    if (board[0][2] == c)
    {
        if (board[1][2] == c)
        {
            if (board[2][2] == c)
            {
                count++;
            }
        }
    }
    // leading diag
    if (board[0][0] == c)
    {
        if (board[1][1] == c)
        {
            if (board[2][2] == c)
            {
                count++;
            }
        }
    }
    // other diag
    if (board[0][2] == c)
    {
        if (board[1][1] == c)
        {
            if (board[2][0] == c)
            {
                count++;
            }
        }
    }

    return count;
}

void solve2()
{
    vector<string> mat;
    for (int i = 0; i < 3; i++)
    {
        string s;
        cin >> s;
        mat.push_back(s);
    }
    int x = 0, y = 0;
    vector<vector<char>> board(3, vector<char>('_', 3));
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            if (mat[j][i] == 'X')
            {
                x++;
            }
            else if (mat[j][i] == 'O')
            {
                y++;
            }
            board[j][i] = mat[j][i];
        }
    }
    int flagX = 0, flagO = 0;
    flagX = check('X', board);
    flagO = check('O', board);
    if (x == y)
    {
        if (flagX == 0 && flagO == 1)
        {
            cout << "YES";
            return;
        }
        else
        {
            cout << "NO";
            return;
        }
    }
    if (x == y + 1)
    {
        if ((flagX == 2 || flagX == 1) && flagO == 0)
        {
            cout << "YES";
            return;
        }
        if (flagX == 0 && flagO == 0)
        {
            cout << "YES";
            return;
        }
    }
    cout << "NO";
}

bool isOValid(vector<vector<char>> &b)
{
    if (b[0][0] == 'O' && b[1][1] == 'O' && b[2][2] == 'O')
        return true;
    else if (b[2][0] == 'O' && b[1][1] == 'O' && b[0][2] == 'O')
        return true;
    else if (b[0][0] == 'O' && b[0][1] == 'O' && b[0][2] == 'O')
        return true;
    else if (b[1][0] == 'O' && b[1][1] == 'O' && b[1][2] == 'O')
        return true;
    else if (b[2][0] == 'O' && b[2][1] == 'O' && b[2][2] == 'O')
        return true;
    else if (b[0][0] == 'O' && b[1][0] == 'O' && b[2][0] == 'O')
        return true;
    else if (b[0][1] == 'O' && b[1][1] == 'O' && b[2][1] == 'O')
        return true;
    else if (b[0][2] == 'O' && b[1][2] == 'O' && b[2][2] == 'O')
        return true;
    else
        return false;
}
bool isXValid(vector<vector<char>> &b)
{
    if (b[0][0] == 'X' && b[1][1] == 'X' && b[2][2] == 'X')
        return true;
    else if (b[2][0] == 'X' && b[1][1] == 'X' && b[0][2] == 'X')
        return true;
    else if (b[0][0] == 'X' && b[0][1] == 'X' && b[0][2] == 'X')
        return true;
    else if (b[1][0] == 'X' && b[1][1] == 'X' && b[1][2] == 'X')
        return true;
    else if (b[2][0] == 'X' && b[2][1] == 'X' && b[2][2] == 'X')
        return true;
    else if (b[0][0] == 'X' && b[1][0] == 'X' && b[2][0] == 'X')
        return true;
    else if (b[0][1] == 'X' && b[1][1] == 'X' && b[2][1] == 'X')
        return true;
    else if (b[0][2] == 'X' && b[1][2] == 'X' && b[2][2] == 'X')
        return true;
    else
        return false;
}
int countRow(vector<vector<char>> &b)
{
    if (b[0][0] == 'X' && b[0][1] == 'X' && b[0][2] == 'X')
    {
        if (b[1][0] == 'O' && b[1][1] == 'O' && b[1][2] == 'O')
        {
            return 0;
        }
        else if (b[2][0] == 'O' && b[2][1] == 'O' && b[2][2] == 'O')
            return 0;
        return 1;
    }
    else if (b[1][0] == 'X' && b[1][1] == 'X' && b[1][2] == 'X')
    {
        if (b[0][0] == 'O' && b[0][1] == 'O' && b[0][2] == 'O')
            return 0;
        else if (b[2][0] == 'O' && b[2][1] == 'O' && b[2][2] == 'O')
            return 0;
        return 1;
    }
    else if (b[2][0] == 'X' && b[2][1] == 'X' && b[2][2] == 'X')
    {
        if (b[1][0] == 'O' && b[1][1] == 'O' && b[1][2] == 'O')
            return 0;
        else if (b[0][0] == 'O' && b[0][1] == 'O' && b[0][2] == 'O')
            return 0;
        return 1;
    }
    if (isOValid(b))
        return -1;
    return 1;
}
int countColumn(vector<vector<char>> &b)
{
    if (b[0][0] == 'X' && b[1][0] == 'X' && b[2][0] == 'X')
    {
        if (b[0][1] == 'O' && b[1][1] == 'O' && b[2][1] == 'O')
            return 0;
        else if (b[0][2] == 'O' && b[1][2] == 'O' && b[2][2] == 'O')
            return 0;
        return 1;
    }
    else if (b[0][1] == 'X' && b[1][1] == 'X' && b[2][1] == 'X')
    {
        if (b[0][0] == 'O' && b[1][0] == 'O' && b[2][0] == 'O')
            return 0;
        else if (b[0][2] == 'O' && b[1][2] == 'O' && b[2][2] == 'O')
            return 0;
        return 1;
    }
    else if (b[0][2] == 'X' && b[1][2] == 'X' && b[2][2] == 'X')
    {
        if (b[0][1] == 'O' && b[1][1] == 'O' && b[2][1] == 'O')
            return 0;
        else if (b[0][0] == 'O' && b[1][0] == 'O' && b[2][0] == 'O')
            return 0;
        return 1;
    }
    if (isOValid(b))
        return -1;
    else
        return 1;
}

bool isValid(vector<vector<char>> &b)
{
    int countx = 0, counto = 0, count = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (b[i][j] == 'X')
                countx++;
            else if (b[i][j] == 'O')
                counto++;
            else if (b[i][j] == '_')
                count++;
        }
    }
    if (counto > countx || countx - counto > 1 || (countx + counto + count) != 9)
        return false;
    else if (isXValid(b))
    {
        if (countx == counto + 1)
            return true;
        return false;
    }
    else if (isOValid(b))
    {
        if (countx != counto)
            return false;
        else
            return true;
    }
    int numberOfRows = countRow(b);
    int numberOfColumns = countColumn(b);
    if (numberOfRows == 0 || numberOfColumns == 0)
        return false;
    else if (numberOfRows == -1 && numberOfColumns == -1)
        return false;
    return true;
}

void solve()
{
    vector<string> mat;
    for (int i = 0; i < 3; i++)
    {
        string s;
        cin >> s;
        mat.push_back(s);
    }
    int x = 0, y = 0;
    vector<vector<char>> board(3, vector<char>('_', 3));
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            if (mat[j][i] == 'X')
            {
                x++;
            }
            else if (mat[j][i] == 'O')
            {
                y++;
            }
            board[j][i] = mat[j][i];
        }
    }
    bool flag = isValid(board);
    if (flag)
    {
        cout << "YES";
        return;
    }
    cout << "NO";
}

void Mario()
{
    int n, m;
    cin >> n >> m;
    vector<vector<char>> arr(n, vector<char>(m, '0'));
    for (int i = 0; i < n; i++)
    {
        string s;
        cin >> s;
        for (int j = 0; j < m; j++)
        {
            arr[i][j] = s[j];
        }
    }
    int cost = 0, coins = 0;
    for (int j = 0; j < m; j++)
    {
        if (arr[n - 1][j] == 'H')
        {
            int k = n - 1;
            while (arr[k][j] != '0')
            {
                cost++;
                k--;
            }
        }
        else
        {
            int k = 0;
            bool flag = true;
            while (k <= n - 1)
            {
                if (arr[k][j] == 'C')
                {
                    coins++;
                    if (flag)
                        cost += n - 1 - k;
                    flag = false;
                }
                k++;
            }
        }
    }
    cout << coins << " " << cost * 2;
}