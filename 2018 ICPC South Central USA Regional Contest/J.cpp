/***
    This algorithm uses ford-fulkerson algorithm combing SPFA to find the
    minimum cost of the flow while making sure that total flow capacity is
    above certain threshold

    The reason to use ford-fulkerson algorithm is because possibly
    we want to undo the previous choice and preserves the total flow
    capacity

    The SPFA is just a random choice to find out the current minimum cost
    from s to t

    The way to do is this is to
      1. always first SPFA find the shortest path from s to t as S_PATH
      2. find the flow capacity of S_PATH as F_PATH
      3. update the residual graph by augmenting path S_PATH with F_PATH
      4. update the total_cost by sending F_PATH with S_PATH, resulting to
         total cost of F_PATH*dist(S_PATH)
      5. repeat this process until no SPFA can be found or require_flow met
***/


#pragma GCC diagnostic ignored "-Wunused-parameter"
#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#define INP "solve"
#define OUT "solve"

const long long INF_LL = 1e17;
const int INF = 1e9 + 100;
const int MOD = 1e9 + 7;
const int Base = 311;
const double EPS = 1e-9;
const int BLOCK = 1000;
const int dx[4] = {0, 0, 1, -1};
const int dy[4] = {1, -1, 0, 0};
mt19937 rnd(chrono::system_clock::now().time_since_epoch().count());

void open_file() {
    #ifdef THEMIS
        freopen(INP ".inp","r",stdin);
        freopen(OUT ".out","w",stdout);
    #endif // THEMIS
}

const int maxN = 2e5 + 100;

#define MCMF MincostMaxflow
#define flow_t int
#define cost_t int
const flow_t foo = (flow_t) 1e9;
const cost_t coo = (cost_t) 1e9;
namespace MincostMaxflow {
    const int maxv = 1e5 + 5;
    const int maxe = 1e6 + 5;
    int n, s, t, E;
    int adj[maxe], nxt[maxe], lst[maxv], frm[maxv], vis[maxv];
    flow_t cap[maxe], flw[maxe], totalFlow;
    cost_t cst[maxe], dst[maxv], totalCost;

    void init(int nn, int ss, int tt) {
        n = nn, s = ss, t = tt;
        // E: current index of add function
        fill_n(lst, n + 7, -1), E = 0;
    }
    void add(int u, int v, flow_t ca, cost_t co) {
        // u: starting node
        // v: destination node
        // E <-> u <-> v 
        // purpose of E: u and v are not in order, but E is; so from E
        //               we can transfrom back to u and v and vice versa
        //               using adj
        // adj: transform from E to u and v
        // cap: capacity forward and backward
        // flw: flow forward and backward initialized from 0
        // cst: weight forward and backward
        // nxt: the next node to tranverse
        adj[E] = v, cap[E] = ca, flw[E] = 0, cst[E] = +co, nxt[E] = lst[u], lst[u] = E++;
        adj[E] = u, cap[E] =  0, flw[E] = 0, cst[E] = -co, nxt[E] = lst[v], lst[v] = E++;
    }
    int spfa() {
        fill_n(dst, n + 7, coo), dst[s] = 0;

        // spfa 
        queue<int> que; que.push(s);
        while (que.size()) {
            int u = que.front(); que.pop();
            for (int e = lst[u]; e != -1; e = nxt[e]) 
                if (flw[e] < cap[e]) {
                    int v = adj[e];
                    if (dst[v] > dst[u] + cst[e]) {
                        dst[v] = dst[u] + cst[e];
                        frm[v] = e;
                        if (!vis[v]) {
                            vis[v] = 1;
                            que.push(v);
                        }
                    }
                }
            vis[u] = 0;
        }
        // find out if the distance(from s to t) < 1e9
        return dst[t] < coo;
    }
    pair<flow_t, cost_t> mincost(int P) {
        totalCost = 0, totalFlow = 0;
        while (1) {
            // if distance(from s to t) >= 1e9 -> break
            if (!spfa()) break;

            // mn: the minimum capacity along the edge
            flow_t mn = foo;
            // frm: parent of the spfa node traversal trajectory
            for (int v = t, e = frm[v]; v != s; v = adj[e ^ 1], e = frm[v]) mn = min(mn, cap[e] - flw[e]);
            
            // update capcities of all edges based on mn
            for (int v = t, e = frm[v]; v != s; v = adj[e ^ 1], e = frm[v]) {
                // flw[e] the forward edge
                flw[e] += mn;
                // flw[e^1] the backward edge 
                flw[e ^ 1] -= mn;
            }

            // consider if totalFlow already sufficiently meets
            // the flow needed (P) no need to proceed the loop
            if(totalFlow + mn >= P){
                totalFlow = P;
                mn = P - totalFlow;
                totalCost += mn * dst[t];
                break;
            }

            // update the maximum flow capacity
            totalFlow += mn;
            
            // update the total cost
            totalCost += mn * dst[t];
        }
        return {totalFlow, totalCost};
    }
}

void sol() {
    int P, R, L;
    cin >> P >> R >> L;
    // s: the starting node
    // t: the terminal node
    int s = R + 2, t = R + 1;
    MCMF::init(R + 3, s, t);
    for (int i = 0; i < L; i++) {
        int u, v;
        cin >> u >> v;
        if (u == -2) u = R;
        if (u == -1) u = R + 1;
        if (v == -2) v = R;
        if (v == -1) v = R + 1;
        MCMF::add(u, v, 1, 1);
        MCMF::add(v, u, 1, 1);
    }
    MCMF::add(s, R, P, 0);
    pair<flow_t, cost_t> ans = MCMF::mincost(P);
    if (ans.first == P) cout << ans.second;
    else cout << (P - ans.first) << " " << "people left behind";
}

void solve() {
    int T;
    //cin >> T;
    T = 1;
    int TestCase = 0;
    while (T--) {
        TestCase += 1;
        //cout << "Case " << TestCase << ":" << ' ';
        ///cout << "Case #" << TestCase << '\n';
        sol();
    }
}
int main(int argc,char *argv[]) {
    //srand(time(nullptr));
    ios_base::sync_with_stdio(0); cin.tie(nullptr); cout.tie(nullptr);
    open_file();
    solve();
}
