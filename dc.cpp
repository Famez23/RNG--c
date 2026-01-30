#include <bits/stdc++.h>
using namespace std;

using ll = long long;

inline bool eq(ll a, ll b) { return a == b; }
inline bool lt(ll a, ll b) { return a < b; }
inline int sgn(ll a) { return (a > 0) - (a < 0); }

struct pt {
    ll x, y;
    pt() {} pt(ll X, ll Y) : x(X), y(Y) {}
    pt operator-(const pt& p) const { return pt(x - p.x, y - p.y); }
    ll cross(const pt& p) const { return x * p.y - y * p.x; }
    ll cross(const pt& a, const pt& b) const { return (a - *this).cross(b - *this); }
    ll dot(const pt& p) const { return x * p.x + y * p.y; }
    ll dot(const pt& a, const pt& b) const { return (a - *this).dot(b - *this); }
    ll sqrLength() const { return this->dot(*this); }
    bool operator==(const pt& p) const { return x == p.x && y == p.y; }
    bool operator<(const pt& p) const { return x < p.x || (x == p.x && y < p.y); }
};
static const pt inf_pt = pt((ll)1e18, (ll)1e18);

struct PtHash {
    size_t operator()(const pt& p) const noexcept {
        uint64_t x = (uint64_t)p.x, y = (uint64_t)p.y;
        x ^= y + 0x9e3779b97f4a7c15ull + (x << 6) + (x >> 2);
        return (size_t)x;
    }
};
struct PtEq {
    bool operator()(const pt& a, const pt& b) const noexcept {
        return a.x == b.x && a.y == b.y;
    }
};

struct QuadEdge;
static std::vector<QuadEdge*> ALL_EDGES;

struct QuadEdge {
    pt origin;
    QuadEdge* rot = nullptr;
    QuadEdge* onext = nullptr;
    bool alive = true; // logical lifetime flag
    QuadEdge* rev() const { return rot->rot; }
    QuadEdge* lnext() const { return rot->rev()->onext->rot; }
    QuadEdge* oprev() const { return rot->onext->rot; }
    QuadEdge* rnext() const { return rot->onext->rot->rot; }
    QuadEdge* rprev() const { return rev()->onext; }
    QuadEdge* dnext() const { return rev()->onext->rev(); }
    QuadEdge* dprev() const { return onext->rev(); }
    pt dest() const { return rev()->origin; }
};

QuadEdge* make_edge(pt from, pt to) {
    QuadEdge* e1 = new QuadEdge;
    QuadEdge* e2 = new QuadEdge;
    QuadEdge* e3 = new QuadEdge;
    QuadEdge* e4 = new QuadEdge;
    e1->origin = from; e2->origin = to; e3->origin = e4->origin = inf_pt;
    e1->rot = e3; e2->rot = e4; e3->rot = e2; e4->rot = e1;
    e1->onext = e1; e2->onext = e2; e3->onext = e4; e4->onext = e3;

    ALL_EDGES.push_back(e1);
    ALL_EDGES.push_back(e2);
    ALL_EDGES.push_back(e3);
    ALL_EDGES.push_back(e4);
    return e1;
}

void splice(QuadEdge* a, QuadEdge* b) {
    swap(a->onext->rot->onext, b->onext->rot->onext);
    swap(a->onext, b->onext);
}

void delete_edge(QuadEdge* e) {
    splice(e, e->oprev());
    splice(e->rev(), e->rev()->oprev());
    QuadEdge* a = e;
    QuadEdge* b = e->rev();
    QuadEdge* ar = a->rot;
    QuadEdge* br = b->rot;
    if (a) a->alive = false;
    if (b) b->alive = false;
    if (ar) ar->alive = false;
    if (br) br->alive = false;
}

QuadEdge* connect(QuadEdge* a, QuadEdge* b) {
    QuadEdge* e = make_edge(a->dest(), b->origin);
    splice(e, a->lnext());
    splice(e->rev(), b);
    return e;
}

inline bool left_of(pt p, QuadEdge* e) { return p.cross(e->origin, e->dest()) > 0; }
inline bool right_of(pt p, QuadEdge* e) { return p.cross(e->origin, e->dest()) < 0; }

template<class T>
T det3(T a1, T a2, T a3, T b1, T b2, T b3, T c1, T c2, T c3) {
    return a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) + a3 * (b1 * c2 - c1 * b2);
}

static inline int orient_sign(const pt& a, const pt& b, const pt& c) {
    __int128 v = (__int128)(b.x - a.x) * (c.y - a.y) - (__int128)(b.y - a.y) * (c.x - a.x);
    return (v > 0) - (v < 0);
}

bool in_circle(pt a, pt b, pt c, pt d) {
    __int128 det = -det3<__int128>(b.x, b.y, b.sqrLength(), c.x, c.y, c.sqrLength(), d.x, d.y, d.sqrLength());
    det += det3<__int128>(a.x, a.y, a.sqrLength(), c.x, c.y, c.sqrLength(), d.x, d.y, d.sqrLength());
    det -= det3<__int128>(a.x, a.y, a.sqrLength(), b.x, b.y, b.sqrLength(), d.x, d.y, d.sqrLength());
    det += det3<__int128>(a.x, a.y, a.sqrLength(), b.x, b.y, b.sqrLength(), c.x, c.y, c.sqrLength());
    if (orient_sign(a, b, c) < 0) det = -det;
    return det > 0;
}

pair<QuadEdge*, QuadEdge*> build_tr(int l, int r, vector<pt>& p) {
    if (r - l + 1 == 2) {
        QuadEdge* res = make_edge(p[l], p[r]);
        return {res, res->rev()};
    }
    if (r - l + 1 == 3) {
        QuadEdge* a = make_edge(p[l], p[l + 1]), *b = make_edge(p[l + 1], p[r]);
        splice(a->rev(), b);
        int sg = sgn(p[l].cross(p[l + 1], p[r]));
        if (sg == 0) return {a, b->rev()};
        QuadEdge* c = connect(b, a);
        return sg == 1 ? make_pair(a, b->rev()) : make_pair(c->rev(), c);
    }
    int mid = (l + r) / 2;
    QuadEdge *ldo, *ldi, *rdo, *rdi;
    tie(ldo, ldi) = build_tr(l, mid, p);
    tie(rdi, rdo) = build_tr(mid + 1, r, p);
    while (true) {
        if (left_of(rdi->origin, ldi)) { ldi = ldi->lnext(); continue; }
        if (right_of(ldi->origin, rdi)) { rdi = rdi->rev()->onext; continue; }
        break;
    }
    QuadEdge* basel = connect(rdi->rev(), ldi);
    auto valid = [&](QuadEdge* e) { return right_of(e->dest(), basel); };
    if (ldi->origin == ldo->origin) ldo = basel->rev();
    if (rdi->origin == rdo->origin) rdo = basel;
    while (true) {
        QuadEdge* lcand = basel->rev()->onext;
        if (valid(lcand)) {
            while (in_circle(basel->dest(), basel->origin, lcand->dest(), lcand->onext->dest())) {
                QuadEdge* t = lcand->onext;
                delete_edge(lcand);
                lcand = t;
            }
        }
        QuadEdge* rcand = basel->oprev();
        if (valid(rcand)) {
            while (in_circle(basel->dest(), basel->origin, rcand->dest(), rcand->oprev()->dest())) {
                QuadEdge* t = rcand->oprev();
                delete_edge(rcand);
                rcand = t;
            }
        }
        if (!valid(lcand) && !valid(rcand)) break;
        if (!valid(lcand) || (valid(rcand) && in_circle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())))
            basel = connect(rcand, basel->rev());
        else
            basel = connect(basel->rev(), lcand->rev());
    }
    return {ldo, rdo};
}

vector<tuple<pt, pt, pt>> delaunay(vector<pt> p) {
    if (p.size() < 3) {
        throw runtime_error("Delaunay triangulation requires at least 3 points");
    }
    sort(p.begin(), p.end()); // uses pt::operator<

    ALL_EDGES.clear();

    QuadEdge *res_first = nullptr, *res_second = nullptr;
    tie(res_first, res_second) = build_tr(0, (int)p.size() - 1, p);

    struct TrioHash {
        size_t operator()(const array<pair<ll,ll>,3>& a) const noexcept {
            auto h = std::hash<ll>{}(a[0].first) ^ (std::hash<ll>{}(a[0].second)<<1);
            h ^= (std::hash<ll>{}(a[1].first)<<1) ^ (std::hash<ll>{}(a[1].second)<<2);
            h ^= (std::hash<ll>{}(a[2].first)<<2) ^ (std::hash<ll>{}(a[2].second)<<3);
            h ^= (h << 13); h ^= (h >> 7); h ^= (h << 17);
            return h;
        }
    };
    struct TrioEq {
        bool operator()(const array<pair<ll,ll>,3>& A, const array<pair<ll,ll>,3>& B) const noexcept {
            return A[0]==B[0] && A[1]==B[1] && A[2]==B[2];
        }
    };
    auto tri_key = [&](const pt& a, const pt& b, const pt& c){
        array<pair<ll,ll>,3> v = {{{a.x,a.y},{b.x,b.y},{c.x,c.y}}};
        sort(v.begin(), v.end());
        return v;
    };

    unordered_set<array<pair<ll,ll>,3>, TrioHash, TrioEq> seen;
    vector<tuple<pt,pt,pt>> result;

    vector<QuadEdge*> cands; cands.reserve(ALL_EDGES.size()*2);
    for (auto e : ALL_EDGES) {
        if (!e || !e->alive) continue;
        cands.push_back(e);
        QuadEdge* r = e->rev();
        if (r && r->alive) cands.push_back(r);
    }

    for (auto e1 : cands) {
        if (!e1) continue;
        if (!e1->alive) continue;
        if (e1->origin==inf_pt || e1->dest()==inf_pt) continue;

        QuadEdge* e2 = e1->lnext(); if (!e2 || !e2->alive) continue;
        if (e2->origin==inf_pt || e2->dest()==inf_pt) continue;

        QuadEdge* e3 = e2->lnext(); if (!e3 || !e3->alive) continue;
        if (e3->origin==inf_pt || e3->dest()==inf_pt) continue;

        if (e3->lnext() != e1) continue;

        const pt& A = e1->origin;
        const pt& B = e2->origin;
        const pt& C = e3->origin;

        auto key = tri_key(A,B,C);
        if (seen.insert(key).second) {
            result.emplace_back(A,B,C);
        }
    }

    for (auto e : ALL_EDGES) {
        delete e;
    }
    ALL_EDGES.clear();

    return result;
}

// ---------- KD-tree ----------
struct KDNode {
    int idx = -1, left = -1, right = -1, axis = 0;
    ll minx, maxx, miny, maxy;
};

struct KDTree {
    const vector<pt>* P = nullptr;
    vector<int> order;
    vector<KDNode> nodes;

    int build_rec(int l, int r) {
        if (l >= r) return -1;
        int m = (l + r) / 2;
        ll minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
        for (int i = l; i < r; ++i) {
            const pt& a = (*P)[order[i]];
            minx = min(minx, a.x);
            maxx = max(maxx, a.x);
            miny = min(miny, a.y);
            maxy = max(maxy, a.y);
        }
        int axis = (maxx - minx >= maxy - miny) ? 0 : 1;
        nth_element(order.begin() + l, order.begin() + m, order.begin() + r, [&](int ia, int ib) {
            return axis == 0 ? ((*P)[ia].x < (*P)[ib].x) : ((*P)[ia].y < (*P)[ib].y);
        });
        KDNode node;
        node.idx = order[m];
        node.axis = axis;
        node.minx = minx;
        node.maxx = maxx;
        node.miny = miny;
        node.maxy = maxy;
        int id = (int)nodes.size();
        nodes.push_back(node);
        nodes[id].left = build_rec(l, m);
        nodes[id].right = build_rec(m + 1, r);
        auto tighten = [&](int child) {
            if (child < 0) return;
            nodes[id].minx = min(nodes[id].minx, nodes[child].minx);
            nodes[id].maxx = max(nodes[id].maxx, nodes[child].maxx);
            nodes[id].miny = min(nodes[id].miny, nodes[child].miny);
            nodes[id].maxy = max(nodes[id].maxy, nodes[child].maxy);
        };
        tighten(nodes[id].left);
        tighten(nodes[id].right);
        return id;
    }

    int root = -1;
    void build(const vector<pt>& pts) {
        P = &pts;
        order.resize(pts.size());
        iota(order.begin(), order.end(), 0);
        nodes.reserve(pts.size() * 2);
        root = build_rec(0, (int)order.size());
    }

    static inline __int128 rect_dist2(const pt& a, const KDNode& n) {
        __int128 dx = 0, dy = 0;
        if (a.x < n.minx) dx = (__int128)n.minx - a.x; 
        else if (a.x > n.maxx) dx = (__int128)a.x - n.maxx;
        if (a.y < n.miny) dy = (__int128)n.miny - a.y; 
        else if (a.y > n.maxy) dy = (__int128)a.y - n.maxy;
        return dx * dx + dy * dy;
    }

    template<typename F>
    void radius_enumerate(int node_id, const pt& center, __int128 r2, F&& visit) const {
        if (node_id < 0) return;
        const KDNode& n = nodes[node_id];
        if (rect_dist2(center, n) >= r2) return;
        const pt& p = (*P)[n.idx];
        __int128 dx = p.x - center.x, dy = p.y - center.y;
        if (dx * dx + dy * dy < r2) visit(n.idx);
        radius_enumerate(n.left, center, r2, visit);
        radius_enumerate(n.right, center, r2, visit);
    }
};

// ---------- RNG extraction ----------
set<pair<pt, pt>> extract_rng(const vector<tuple<pt, pt, pt>>& triangles, 
                               const vector<pt>& all_points) {
    unordered_map<pt, unordered_set<pt, PtHash, PtEq>, PtHash, PtEq> neighbors;
    neighbors.reserve(triangles.size() * 2 + 1);

    for (const auto& t : triangles) {
        pt a, b, c;
        tie(a, b, c) = t;
        neighbors[a].insert(b); neighbors[b].insert(a);
        neighbors[b].insert(c); neighbors[c].insert(b);
        neighbors[c].insert(a); neighbors[a].insert(c);
    }

    KDTree kdt;
    kdt.build(all_points);

    vector<pair<pt, pt>> candidates;
    candidates.reserve(neighbors.size() * 3);
    for (const auto& kv : neighbors) {
        const pt& u = kv.first;
        for (const pt& v : kv.second) {
            if (u < v) candidates.emplace_back(u, v);
        }
    }

    set<pair<pt, pt>> rng_edges;

    for (const auto& uv : candidates) {
        const pt& u = uv.first;
        const pt& v = uv.second;
        __int128 dx = u.x - v.x, dy = u.y - v.y;
        __int128 r2 = dx * dx + dy * dy;
        bool keep = true;

        auto check_witnesses = [&](const pt& center) {
            if (!keep) return;
            kdt.radius_enumerate(kdt.root, center, r2, [&](int k) {
                if (!keep) return;
                const pt& w = all_points[k];
                if (w == u || w == v) return;
                __int128 dx_u = w.x - u.x, dy_u = w.y - u.y;
                __int128 d2_u = dx_u * dx_u + dy_u * dy_u;
                __int128 dx_v = w.x - v.x, dy_v = w.y - v.y;
                __int128 d2_v = dx_v * dx_v + dy_v * dy_v;
                if (d2_u < r2 && d2_v < r2) keep = false;
            });
        };
        
        check_witnesses(u);
        check_witnesses(v);

        if (keep) rng_edges.insert(uv);
    }
    
    return rng_edges;
}

static inline uint64_t mix64(uint64_t x,uint64_t y){ x^= y + 0x9e3779b97f4a7c15ull + (x<<6)+(x>>2); return x; }

vector<pt> load_points_from_file(const string& filename){
    ifstream in(filename); if(!in) throw runtime_error("Failed to open file: "+filename);
    unordered_set<uint64_t> seen; vector<pt> pts; string line; size_t bad=0;
    while(getline(in,line)){
        if(line.empty()) continue; stringstream ss(line); long long x,y; char c;
        if((ss>>x>>c>>y) && c==','){ uint64_t k=mix64((uint64_t)x,(uint64_t)y); if(seen.insert(k).second) pts.emplace_back(x,y); }
        else ++bad;
    }
    if(bad) cerr<<"Warning: skipped "<<bad<<" malformed line(s)\n";
    return pts;
}
static size_t convex_hull_size(vector<pt> v){
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    int n=(int)v.size(); if(n<=1) return n;
    vector<pt> H(2*n);
    int k=0;
    auto cross = [](const pt& O, const pt& A, const pt& B){
        __int128 x1=A.x-O.x, y1=A.y-O.y, x2=B.x-O.x, y2=B.y-O.y;
        __int128 c = x1*y2 - y1*x2;
        return (c>0)-(c<0);
    };
    for(int i=0;i<n;i++){
        while(k>=2 && cross(H[k-2],H[k-1],v[i])<=0) k--;
        H[k++]=v[i];
    }
    for(int i=n-2, t=k+1;i>=0;i--){
        while(k>=t && cross(H[k-2],H[k-1],v[i])<=0) k--;
        H[k++]=v[i];
    }
    return (size_t)(k-1);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        auto start = chrono::high_resolution_clock::now();

        vector<pt> points = load_points_from_file("points.csv");
        cout << "Loaded " << points.size() << " points\n";

        auto start_DT = chrono::high_resolution_clock::now();
        cout << "Execution time: " << chrono::duration<double>(start_DT - start).count() << " seconds\n";
        auto triangles = delaunay(points);
        auto end_DT = chrono::high_resolution_clock::now();
        cout << "Number of triangles: " << triangles.size() << "\n";
        cout << "Execution time: " << chrono::duration<double>(end_DT - start_DT).count() << " seconds\n";

        auto start_rng = chrono::high_resolution_clock::now();
        auto rng_edges = extract_rng(triangles, points);
        auto end_rng = chrono::high_resolution_clock::now();
        cout << "Number of RNG edges: " << rng_edges.size() << "\n";
        cout << "Execution time: " << chrono::duration<double>(end_rng - start_rng).count() << " seconds\n";
        

        auto end = chrono::high_resolution_clock::now();
        cout << "Execution time: " << chrono::duration<double>(end - start).count() << " seconds\n";
        cout << "--------------------------------\n";
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
