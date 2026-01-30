#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <array>
#include <chrono>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/scatter.h>
#include <thrust/gather.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cmath>
#include <cuda_profiler_api.h>

#ifdef _OPENMP
  #include <omp.h>
#endif
#ifndef BASE_N
#define BASE_N            256       
#endif
#ifndef TASK_CUTOFF
#define TASK_CUTOFF       50000   
#endif
#ifndef BOTH_TASK_CUTOFF
#define BOTH_TASK_CUTOFF  400000   
#endif 

using namespace std;

using ll = long long;

inline bool eq(ll a, ll b) { return a == b; }
inline bool lt(ll a, ll b) { return a < b; }
inline int sgn(ll a) { return (a > 0) - (a < 0); }

struct pt {
    ll x, y;
    __host__ __device__ pt() : x(0), y(0) {}
    __host__ __device__ pt(ll X, ll Y) : x(X), y(Y) {}
    pt operator-(const pt& p) const { return pt(x - p.x, y - p.y); }
    ll cross(const pt& p) const { return x * p.y - y * p.x; }
    ll cross(const pt& a, const pt& b) const { return (a - *this).cross(b - *this); }
    ll dot(const pt& p) const { return x * p.x + y * p.y; }
    ll dot(const pt& a, const pt& b) const { return (a - *this).dot(b - *this); }
    ll sqrLength() const { return this->dot(*this); }
    __host__ __device__ bool operator==(const pt& p) const { return x == p.x && y == p.y; }
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
static inline uint64_t pack_pair_ids(int a, int b) {
    uint32_t lo = (uint32_t)std::min(a,b);
    uint32_t hi = (uint32_t)std::max(a,b);
    return ( (uint64_t)hi << 32 ) | (uint64_t)lo;
}
static inline void unpack_pair_ids(uint64_t k, int &a, int &b) {
    uint32_t lo = (uint32_t)(k & 0xFFFFFFFFull);
    uint32_t hi = (uint32_t)(k >> 32);
    a = (int)lo; b = (int)hi;
}
struct Edge {
    pt u, v;
    __host__ __device__ Edge() {}
    __host__ __device__ Edge(pt a, pt b) : u(a), v(b) {}
};
struct TriI { int a, b, c; };
int find_idx(const std::vector<std::pair<pt,int>>& vp, const pt& p)
{
    auto it = std::lower_bound(
        vp.begin(), vp.end(), p,
        [](auto const& pair, pt const& val){ return pair.first < val; }
    );

    if (it != vp.end() && it->first == p)
        return it->second;

    return -1; 
}

std::vector<TriI>
triangles_to_indices(const std::vector<std::tuple<pt,pt,pt>>& tris,
                     const std::vector<pt>& points)
{
    std::vector<std::pair<pt,int>> vp;
    vp.reserve(points.size());
    for (int i = 0; i < (int)points.size(); ++i)
        vp.emplace_back(points[i], i);

    std::sort(vp.begin(), vp.end(),
              [](auto const& A, auto const& B){ return A.first < B.first; });

    std::vector<TriI> out(tris.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)tris.size(); ++i) {
        auto& t = tris[i];
        pt a,b,c;
        std::tie(a,b,c) = t;

        int ia = find_idx(vp, a);
        int ib = find_idx(vp, b);
        int ic = find_idx(vp, c);

        out[i] = { ia, ib, ic };
    }
    return out;
}



struct QuadEdge;
static std::vector<QuadEdge*> ALL_EDGES;
thread_local std::vector<QuadEdge*> TLS_EDGES;

struct QuadEdge {
    pt origin;
    QuadEdge* rot = nullptr;
    QuadEdge* onext = nullptr;
    bool alive = true;
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

    e1->origin = from;  e2->origin = to;   e3->origin = inf_pt; e4->origin = inf_pt;

    e1->rot = e3;  e2->rot = e4;  e3->rot = e2;  e4->rot = e1;
    e1->onext = e1; e2->onext = e2; e3->onext = e4; e4->onext = e3;

    e1->alive = true; e2->alive = true; e3->alive = true; e4->alive = true;

    TLS_EDGES.push_back(e1);
    TLS_EDGES.push_back(e2);
    TLS_EDGES.push_back(e3);
    TLS_EDGES.push_back(e4);

    return e1;
}

void splice(QuadEdge* a, QuadEdge* b) {
    using std::swap;
    swap(a->onext->rot->onext, b->onext->rot->onext);
    swap(a->onext, b->onext);
}

inline void delete_edge(QuadEdge* e) {
    if (!e || !e->alive) return;
    splice(e, e->oprev());
    splice(e->rev(), e->rev()->oprev());
    QuadEdge* a  = e;
    QuadEdge* b  = e->rev();
    QuadEdge* ar = a->rot;
    QuadEdge* br = b->rot;
    if (a)  a->alive  = false;
    if (b)  b->alive  = false;
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
static std::atomic<long> g_task_cnt{0};
static std::atomic<long> g_big_spawns{0};
static std::atomic<long> g_mid_spawns{0};

void merge_tls_edges_into_all(std::vector<QuadEdge*>& ALL_EDGES) {
    int T = 1;
    #pragma omp parallel
    {
        #pragma omp single
        T = omp_get_num_threads();
    }

    std::vector<size_t> sizes(T, 0), offs(T+1, 0);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        sizes[tid] = TLS_EDGES.size();
    }
    for (int t = 0; t < T; ++t) offs[t+1] = offs[t] + sizes[t];

    ALL_EDGES.resize(offs[T]);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::copy(TLS_EDGES.begin(), TLS_EDGES.end(),
                  ALL_EDGES.begin() + offs[tid]);
        std::vector<QuadEdge*>().swap(TLS_EDGES); 
    }
}

static inline std::pair<QuadEdge*, QuadEdge*>
build_tr_serial(int l, int r, std::vector<pt>& p)
{
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
        return sg == 1 ? std::make_pair(a, b->rev()) : std::make_pair(c->rev(), c);
    }

    int mid = (l + r) / 2;
    QuadEdge *ldo, *ldi, *rdo, *rdi;
    std::tie(ldo, ldi) = build_tr_serial(l, mid, p);
    std::tie(rdi, rdo) = build_tr_serial(mid + 1, r, p);

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
            while (in_circle(basel->dest(), basel->origin,
                             lcand->dest(), lcand->onext->dest())) {
                QuadEdge* t = lcand->onext;
                delete_edge(lcand);   // lazy unlink
                lcand = t;
            }
        }
        QuadEdge* rcand = basel->oprev();
        if (valid(rcand)) {
            while (in_circle(basel->dest(), basel->origin,
                             rcand->dest(), rcand->oprev()->dest())) {
                QuadEdge* t = rcand->oprev();
                delete_edge(rcand);   // lazy unlink
                rcand = t;
            }
        }
        if (!valid(lcand) && !valid(rcand)) break;

        if (!valid(lcand) ||
            (valid(rcand) &&
             in_circle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())))
            basel = connect(rcand, basel->rev());
        else
            basel = connect(basel->rev(), lcand->rev());
    }
    return {ldo, rdo};
}

static inline std::pair<QuadEdge*, QuadEdge*>
build_tr_task(int l, int r, std::vector<pt>& p)
{
    int n = r - l + 1;
    if (n <= BASE_N) return build_tr_serial(l, r, p);

    int mid = (l + r) / 2;
    std::pair<QuadEdge*, QuadEdge*> L, R;

    #pragma omp taskgroup
    {
        if (n >= BOTH_TASK_CUTOFF) {
            #pragma omp task shared(L,p) firstprivate(l,mid)
            { L = build_tr_task(l, mid, p); }
            g_big_spawns.fetch_add(1, std::memory_order_relaxed);
            g_task_cnt.fetch_add(1, std::memory_order_relaxed);

            #pragma omp task shared(R,p) firstprivate(mid,r)
            { R = build_tr_task(mid+1, r, p); }
            g_big_spawns.fetch_add(1, std::memory_order_relaxed);
            g_task_cnt.fetch_add(1, std::memory_order_relaxed);

            #pragma omp taskwait
        } else if (n >= TASK_CUTOFF) {
            #pragma omp task shared(L,p) firstprivate(l,mid)
            { L = build_tr_task(l, mid, p); }
            g_mid_spawns.fetch_add(1, std::memory_order_relaxed);
            g_task_cnt.fetch_add(1, std::memory_order_relaxed);

            R = build_tr_task(mid+1, r, p);
            #pragma omp taskwait
        } else {
            L = build_tr_task(l, mid, p);
            R = build_tr_task(mid+1, r, p);
        }
    }

    // ---- merge ----
    QuadEdge *ldo = L.first, *ldi = L.second;
    QuadEdge *rdi = R.first, *rdo = R.second;

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
            while (in_circle(basel->dest(), basel->origin,
                             lcand->dest(), lcand->onext->dest())) {
                QuadEdge* t = lcand->onext;
                delete_edge(lcand);
                lcand = t;
            }
        }
        QuadEdge* rcand = basel->oprev();
        if (valid(rcand)) {
            while (in_circle(basel->dest(), basel->origin,
                             rcand->dest(), rcand->oprev()->dest())) {
                QuadEdge* t = rcand->oprev();
                delete_edge(rcand);
                rcand = t;
            }
        }
        if (!valid(lcand) && !valid(rcand)) break;

        if (!valid(lcand) ||
            (valid(rcand) &&
             in_circle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())))
            basel = connect(rcand, basel->rev());
        else
            basel = connect(basel->rev(), lcand->rev());
    }
    return {ldo, rdo};
}

inline std::pair<QuadEdge*, QuadEdge*>
build_tr(int l, int r, std::vector<pt>& p)
{
    std::pair<QuadEdge*, QuadEdge*> res;
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            res = build_tr_task(l, r, p);
        }
    }
    return res;
}



vector<tuple<pt, pt, pt>> delaunay(vector<pt> p) {
    if (p.size() < 3) {
        throw std::runtime_error("Delaunay triangulation requires at least 3 points");
    }
std::cerr << "OMP threads: " << omp_get_max_threads() << "\n";
    auto t00 = std::chrono::high_resolution_clock::now();
    std::sort(p.begin(), p.end());
    auto t002 = std::chrono::high_resolution_clock::now();
    cout << "Time taken for sorting points: "
         << std::chrono::duration<double,std::milli>(t002 - t00).count()
         << " milliseconds" << endl;

    ALL_EDGES.clear();

    std::pair<QuadEdge*, QuadEdge*> res;
    auto t01 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            res = build_tr_task(0, (int)p.size() - 1, p);
        }
    }   
    std::cerr << "tasks created: " << g_task_cnt.load() << "\n";
 std::cerr << "tasks created: " << g_task_cnt.load()
          << "  big:" << g_big_spawns.load()
          << "  mid:" << g_mid_spawns.load() << "\n";

    merge_tls_edges_into_all(ALL_EDGES);
    auto t02 = std::chrono::high_resolution_clock::now();
    cout << "Time taken for build_tr: "
         << std::chrono::duration<double>(t02 - t01).count()
         << " seconds" << endl;

    auto t1 = std::chrono::high_resolution_clock::now();

struct Tri { pt a, b, c; };
auto lex_lt = [](const pt& p, const pt& q){
    return (p.x < q.x) || (p.x == q.x && p.y < q.y);
};

auto floop = std::chrono::high_resolution_clock::now();

std::vector<size_t> local_counts(omp_get_max_threads(), 0);

#pragma omp parallel
{
    int tid = omp_get_thread_num();
    size_t count = 0;

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < ALL_EDGES.size(); ++i) {
        QuadEdge* e = ALL_EDGES[i];
        if (!e || !e->alive) continue;

        ++count;                            
        if (QuadEdge* r = e->rev(); r && r->alive)
            ++count;                        
    }
    local_counts[tid] = count;
}

std::vector<size_t> offsets = local_counts;
for (int i = 1; i < offsets.size(); ++i)
    offsets[i] += offsets[i-1];

size_t total_count = offsets.back();

std::vector<QuadEdge*> cands(total_count);

#pragma omp parallel
{
    int tid = omp_get_thread_num();
    size_t start = (tid == 0) ? 0 : offsets[tid-1];
    size_t pos   = start;

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < ALL_EDGES.size(); ++i) {
        QuadEdge* e = ALL_EDGES[i];
        if (!e || !e->alive) continue;

        cands[pos++] = e;
        if (QuadEdge* r = e->rev(); r && r->alive)
            cands[pos++] = r;
    }
}

auto tloop = std::chrono::high_resolution_clock::now();

std::cout << "Time taken for collecting candidate edges (OpenMP, 2-pass): "
          << std::chrono::duration<double>(tloop - floop).count() << " s\n";
size_t count = cands.size();

auto tloop2 = std::chrono::high_resolution_clock::now();

const int T = omp_get_max_threads();

std::vector<size_t> tcounts(T, 0);

#pragma omp parallel
{
    const int tid = omp_get_thread_num();
    size_t local = 0;

    #pragma omp for schedule(static)
    for (ptrdiff_t i = 0; i < (ptrdiff_t)count; ++i) {
        QuadEdge* e1 = cands[i];
        if (!e1 || !e1->alive) continue;
        if (e1->origin==inf_pt || e1->dest()==inf_pt) continue;

        QuadEdge* e2 = e1->lnext(); if (!e2 || !e2->alive) continue;
        if (e2->origin==inf_pt || e2->dest()==inf_pt) continue;

        QuadEdge* e3 = e2->lnext(); if (!e3 || !e3->alive) continue;
        if (e3->origin==inf_pt || e3->dest()==inf_pt) continue;

        if (e3->lnext() != e1) continue;  

        const pt& A = e1->origin;
        const pt& B = e2->origin;
        const pt& C = e3->origin;

        if (!lex_lt(B,A) && !lex_lt(C,A)) {
            ++local; 
        }
    }
    tcounts[tid] = local;
}

std::vector<size_t> offs(T+1, 0);
for (int t = 0; t < T; ++t) offs[t+1] = offs[t] + tcounts[t];

std::vector<std::tuple<pt,pt,pt>> result;
result.resize(offs[T]);

#pragma omp parallel
{
    const int tid = omp_get_thread_num();
    size_t out_idx = offs[tid];

    #pragma omp for schedule(static)
    for (ptrdiff_t i = 0; i < (ptrdiff_t)count; ++i) {
        QuadEdge* e1 = cands[i];
        if (!e1 || !e1->alive) continue;
        if (e1->origin==inf_pt || e1->dest()==inf_pt) continue;

        QuadEdge* e2 = e1->lnext(); if (!e2 || !e2->alive) continue;
        if (e2->origin==inf_pt || e2->dest()==inf_pt) continue;

        QuadEdge* e3 = e2->lnext(); if (!e3 || !e3->alive) continue;
        if (e3->origin==inf_pt || e3->dest()==inf_pt) continue;

        if (e3->lnext() != e1) continue;

        const pt& A = e1->origin;
        const pt& B = e2->origin;
        const pt& C = e3->origin;

        if (!lex_lt(B,A) && !lex_lt(C,A)) {
            result[out_idx++] = std::make_tuple(A,B,C);
        }
    }
}

auto tloop3 = std::chrono::high_resolution_clock::now();
std::cout << "Time taken for processing candidate edges (parallel, single-emission): "
          << std::chrono::duration<double>(tloop3 - tloop2).count()
          << " seconds\n";

#pragma omp parallel for schedule(static)
for (ptrdiff_t i = 0; i < (ptrdiff_t)ALL_EDGES.size(); ++i) {
    QuadEdge* e = ALL_EDGES[i];
    if (e) delete e;          
}
ALL_EDGES.clear();

return result;

}


#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " \
                      << cudaGetErrorString(err) << "\n"; \
            exit(EXIT_FAILURE); \
        } \
    } while(0)
__device__ __forceinline__ pt load_pt_ldg(const pt* p) {
    pt tmp;
    tmp.x = __ldg(&p->x);
    tmp.y = __ldg(&p->y);
    return tmp;
}
__global__ void verify_rng_edges_spatial_kernel_idx(
    const pt* __restrict__ points,           
    const int n_points,
    const int* __restrict__ cand_u,
    const int* __restrict__ cand_v,
    const int n_candidates,
    uint8_t* __restrict__ keep_flags,
    const int* __restrict__ cell_starts,
    const int* __restrict__ cell_counts,
    const int* __restrict__ point_indices,
    long long min_x, long long min_y, long long cell_size,
    const int grid_width, const int grid_height)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_candidates) return;

    const int u_idx = __ldg(&cand_u[idx]);
    const int v_idx = __ldg(&cand_v[idx]);
    const pt u = load_pt_ldg(&points[u_idx]);
    const pt v = load_pt_ldg(&points[v_idx]);

    // Use exact integer arithmetic for distance comparisons
    const long long dx = u.x - v.x;
    const long long dy = u.y - v.y;
    const long long r2 = dx*dx + dy*dy;

    if (r2 <= 0) {
        keep_flags[idx] = 1;
        return;
    }

    // Use float only for approximate cell radius calculation
    const float r = sqrtf((float)r2);
    const float inv_cell = 1.0f / (float)cell_size;
    const int cell_radius = (int)(r * inv_cell) + 1;

    const int cu_x = (int)floorf(((float)((long long)u.x - min_x)) * inv_cell);
    const int cu_y = (int)floorf(((float)((long long)u.y - min_y)) * inv_cell);

    const int max_rx = max(cu_x, grid_width  - 1 - cu_x);
    const int max_ry = max(cu_y, grid_height - 1 - cu_y);
    const int max_r  = min(cell_radius, max(max_rx, max_ry));

    uint8_t keep = 1;

    for (int rring = 0; rring <= max_r && keep; ++rring) {
        for (int dy_cell = -rring; dy_cell <= rring && keep; ++dy_cell) {
            int abs_dy = (dy_cell < 0) ? -dy_cell : dy_cell;
            for (int dx_cell = -rring; dx_cell <= rring && keep; ++dx_cell) {
                int abs_dx = (dx_cell < 0) ? -dx_cell : dx_cell;
                if (abs_dx != rring && abs_dy != rring) continue; 

                int cell_x = cu_x + dx_cell;
                int cell_y = cu_y + dy_cell;
                if ((unsigned)cell_x >= (unsigned)grid_width ||
                    (unsigned)cell_y >= (unsigned)grid_height) continue;

                int cell_id = cell_y * grid_width + cell_x;
                int count = __ldg(&cell_counts[cell_id]);
                if (count <= 0) continue;

                int start = __ldg(&cell_starts[cell_id]);
                for (int ii = 0; ii < count && keep; ++ii) {
                    int sorted_idx = start + ii;
                    
                    // Compare in SORTED space
                    if (sorted_idx == u_idx || sorted_idx == v_idx) continue;
                    
                    const pt w = load_pt_ldg(&points[sorted_idx]);
                    
                    // Exact integer distance comparisons
                    const long long dxu = w.x - u.x;
                    const long long dyu = w.y - u.y;
                    const long long d2u = dxu*dxu + dyu*dyu;
                    if (d2u >= r2) continue;
                    
                    const long long dxv = w.x - v.x;
                    const long long dyv = w.y - v.y;
                    const long long d2v = dxv*dxv + dyv*dyv;
                    if (d2v < r2) {
                        keep = 0;
                        break; 
                    }
                }
            }
        }
    }

    keep_flags[idx] = keep;
}
__global__
void convert_candidates(const int* cu, const int* cv,
                        int* cu_new, int* cv_new,
                        const int* old2new,
                        int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    cu_new[i] = old2new[ cu[i] ];
    cv_new[i] = old2new[ cv[i] ];
}

__global__
void compute_cell_id_kernel(const pt* points, int* cell_id,
                            int n_points, long long min_x, long long min_y,
                            long long cell_size, int grid_width)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_points) return;

    long long x = points[i].x;
    long long y = points[i].y;

    int cx = (x - min_x) / cell_size;
    int cy = (y - min_y) / cell_size;

    cell_id[i] = cy * grid_width + cx;
}

struct SpatialGrid {
    int* cell_starts;  
    int* cell_counts;  
    int* point_indices;
    pt*  sorted_points; 
    ll min_x, min_y;
    ll cell_size;
    int grid_width, grid_height;
    int total_cells;
};


SpatialGrid build_spatial_grid(const std::vector<pt>& points,
                               thrust::device_vector<pt>& d_sorted_points_out,
                               thrust::device_vector<int>& d_point_indices_out,
                               thrust::device_vector<int>& d_old2new_out)
{
    SpatialGrid grid{};
    const int n = (int)points.size();
    if (n <= 0) return grid;

    ll min_x = points[0].x, max_x = points[0].x;
    ll min_y = points[0].y, max_y = points[0].y;
    for (int i = 1; i < n; ++i) {
        ll x = points[i].x, y = points[i].y;
        if (x < min_x) min_x = x; if (x > max_x) max_x = x;
        if (y < min_y) min_y = y; if (y > max_y) max_y = y;
    }

    ll range_x = max_x - min_x + 1;
    ll range_y = max_y - min_y + 1;

    const double P = 32.0;
    double cells_target = std::max(1.0, double(n)/P);
    int side = int(std::ceil(std::sqrt(cells_target)));
    ll cell_size = std::max<ll>(1, (std::max(range_x, range_y) + side - 1)/side);

    grid.min_x = min_x;
    grid.min_y = min_y;
    grid.cell_size = cell_size;

    grid.grid_width  = (range_x + cell_size - 1) / cell_size;
    grid.grid_height = (range_y + cell_size - 1) / cell_size;
    grid.total_cells = grid.grid_width * grid.grid_height;

    thrust::device_vector<pt> d_points(points.begin(), points.end());

    thrust::device_vector<int> d_cell_id(n);
    thrust::device_vector<int> d_point_id(n);
    thrust::sequence(d_point_id.begin(), d_point_id.end());

    int TPB = 256;
    int BLK = (n + TPB - 1) / TPB;

    compute_cell_id_kernel<<<BLK, TPB>>>(
        thrust::raw_pointer_cast(d_points.data()),
        thrust::raw_pointer_cast(d_cell_id.data()),
        n,
        min_x, min_y,
        cell_size,
        grid.grid_width
    );
    cudaDeviceSynchronize();

    thrust::sort_by_key(d_cell_id.begin(), d_cell_id.end(), d_point_id.begin());

    d_old2new_out.resize(n);
    thrust::device_vector<int> tmp(n);
    thrust::sequence(tmp.begin(), tmp.end());

    thrust::scatter(
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(n),
        d_point_id.begin(),
        d_old2new_out.begin()
    );

    d_sorted_points_out.resize(n);
    thrust::gather(
        d_point_id.begin(),
        d_point_id.end(),
        d_points.begin(),
        d_sorted_points_out.begin()
    );

    d_point_indices_out = d_point_id;

    thrust::host_vector<int> h_cell_id = d_cell_id;

    std::vector<int> cell_starts(grid.total_cells, -1);
    std::vector<int> cell_counts(grid.total_cells, 0);

    for (int i = 0; i < n; ++i) {
        int cid = h_cell_id[i];
        if (cell_starts[cid] < 0) cell_starts[cid] = i;
        cell_counts[cid]++;
    }

    cudaMalloc(&grid.cell_starts, grid.total_cells*sizeof(int));
    cudaMalloc(&grid.cell_counts, grid.total_cells*sizeof(int));
    cudaMalloc(&grid.point_indices, n*sizeof(int));
    cudaMalloc(&grid.sorted_points, n*sizeof(pt));

    cudaMemcpy(grid.cell_starts, cell_starts.data(),
               grid.total_cells*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(grid.cell_counts, cell_counts.data(),
               grid.total_cells*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(grid.point_indices,
               thrust::raw_pointer_cast(d_point_indices_out.data()),
               n*sizeof(int),
               cudaMemcpyDeviceToDevice);

    cudaMemcpy(grid.sorted_points,
               thrust::raw_pointer_cast(d_sorted_points_out.data()),
               n*sizeof(pt),
               cudaMemcpyDeviceToDevice);

    return grid;
}

void free_spatial_grid(SpatialGrid& grid) {
    cudaFree(grid.cell_starts);
    cudaFree(grid.cell_counts);
    cudaFree(grid.point_indices);
    cudaFree(grid.sorted_points);
}

std::vector<std::pair<pt,pt>>
extract_rng_cuda_vec_idx(const std::vector<TriI>& triangles_i,
                         const std::vector<pt>& all_points)
{
    using clk = std::chrono::high_resolution_clock;
    auto t0 = clk::now();

    const int n_points = (int)all_points.size();
    const int n_tris   = (int)triangles_i.size();

    auto tB0 = clk::now();
    std::vector<uint64_t> keys((size_t)n_tris * 3);

    #pragma omp parallel for schedule(static)
    for (int t = 0; t < n_tris; ++t) {
        const auto &tr = triangles_i[t];
        size_t base = (size_t)t * 3;
        keys[base + 0] = pack_pair_ids(tr.a, tr.b);
        keys[base + 1] = pack_pair_ids(tr.b, tr.c);
        keys[base + 2] = pack_pair_ids(tr.c, tr.a);
    }
    auto tB1 = clk::now();


    auto tC0 = clk::now();
    thrust::device_vector<uint64_t> d_keys(keys.begin(), keys.end());
    thrust::sort(d_keys.begin(), d_keys.end());
    d_keys.erase(thrust::unique(d_keys.begin(), d_keys.end()), d_keys.end());
    std::vector<uint64_t> uniq_keys(d_keys.size());
    thrust::copy(d_keys.begin(), d_keys.end(), uniq_keys.begin());
    auto tC1 = clk::now();


    auto tD0 = clk::now();
    const int n_candidates = (int)uniq_keys.size();
    std::vector<int> cand_u(n_candidates), cand_v(n_candidates);

    for (int i = 0; i < n_candidates; ++i) {
        int u,v; unpack_pair_ids(uniq_keys[i], u, v);
        cand_u[i] = u;
        cand_v[i] = v;
    }
    auto tD1 = clk::now();

    std::cout << "GPU: Checking " << n_candidates
              << " candidates with " << n_points << " points\n";


    auto tE0 = clk::now();

    thrust::device_vector<pt> d_sorted_points;
    thrust::device_vector<int> d_point_indices;
    thrust::device_vector<int> d_old2new;

    SpatialGrid grid = build_spatial_grid(
        all_points,
        d_sorted_points,
        d_point_indices,
        d_old2new
    );

    auto tE1 = clk::now();

    auto tF0 = clk::now();

    thrust::device_vector<int> d_cu(cand_u.begin(), cand_u.end());
    thrust::device_vector<int> d_cv(cand_v.begin(), cand_v.end());

    thrust::device_vector<int> d_cu_new(n_candidates);
    thrust::device_vector<int> d_cv_new(n_candidates);

    int TPB = 256;
    int BLK = (n_candidates + TPB - 1) / TPB;

    convert_candidates<<<BLK, TPB>>>(
        thrust::raw_pointer_cast(d_cu.data()),
        thrust::raw_pointer_cast(d_cv.data()),
        thrust::raw_pointer_cast(d_cu_new.data()),
        thrust::raw_pointer_cast(d_cv_new.data()),
        thrust::raw_pointer_cast(d_old2new.data()),
        n_candidates
    );
    cudaDeviceSynchronize();

    thrust::device_vector<uint8_t> d_keep(n_candidates);
    auto tF1 = clk::now();


    verify_rng_edges_spatial_kernel_idx<<<BLK, TPB>>>(
        grid.sorted_points,      
        n_points,
        thrust::raw_pointer_cast(d_cu_new.data()),
        thrust::raw_pointer_cast(d_cv_new.data()),
        n_candidates,
        thrust::raw_pointer_cast(d_keep.data()),
        grid.cell_starts,
        grid.cell_counts,
        grid.point_indices,
        grid.min_x, grid.min_y,
        grid.cell_size,
        grid.grid_width, grid.grid_height
    );

    auto tG0 = clk::now();
    std::vector<uint8_t> keep(n_candidates);
    thrust::copy(d_keep.begin(), d_keep.end(), keep.begin());

    free_spatial_grid(grid);
    auto tG1 = clk::now();


    auto tH0 = clk::now();
    std::vector<std::pair<pt,pt>> out;
    out.reserve(n_candidates);

    for (int i = 0; i < n_candidates; ++i) {
        if (keep[i]) {
            out.emplace_back(
                all_points[cand_u[i]],
                all_points[cand_v[i]]
            );
        }
    }
    auto tH1 = clk::now();

    auto ms = [](auto a, auto b){
        return std::chrono::duration<double,std::milli>(b-a).count();
    };

    std::cout << "prep.emit_edges(idx) = " << ms(tB0,tB1) << " ms\n";
    std::cout << "prep.sort_unique(GPU) = " << ms(tC0,tC1) << " ms\n";
    std::cout << "prep.candidates(idx)  = " << ms(tD0,tD1) << " ms\n";
    std::cout << "prep.grid(build)      = " << ms(tE0,tE1) << " ms\n";
    std::cout << "prep.H2D+convert      = " << ms(tF0,tF1) << " ms\n";
    std::cout << "post.D2H+cleanup      = " << ms(tG0,tG1) << " ms\n";
    std::cout << "post.build_vec        = " << ms(tH0,tH1) << " ms\n";

    std::cout << "Total extract_rng_cuda_vec_idx time: "
              << std::chrono::duration<double>(clk::now()-t0).count()
              << " s\n";

    return std::move(out);
}


vector<pt> load_points_from_file(const string& filename) {
    std::ifstream in(filename);
    if (!in.is_open())
        throw std::runtime_error("Failed to open file: " + filename);

    std::vector<std::string> lines;
    lines.reserve(1000);

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty())
            lines.push_back(std::move(line));
    }

    std::vector<pt> pts(lines.size());

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string& s = lines[i];
        size_t c = s.find(',');
        if (c != std::string::npos) {
            long long x = std::stoll(s.substr(0, c));
            long long y = std::stoll(s.substr(c + 1));
            pts[i] = pt(x, y);
        }
    }

    return pts;
}


// static size_t convex_hull_size(vector<pt> v){
//     sort(v.begin(), v.end());
//     v.erase(unique(v.begin(), v.end()), v.end());
//     int n=(int)v.size(); if(n<=1) return n;
//     vector<pt> H(2*n);
//     int k=0;
//     auto cross = [](const pt& O, const pt& A, const pt& B){
//         __int128 x1=A.x-O.x, y1=A.y-O.y, x2=B.x-O.x, y2=B.y-O.y;
//         __int128 c = x1*y2 - y1*x2;
//         return (c>0)-(c<0);
//     };
//     for(int i=0;i<n;i++){
//         while(k>=2 && cross(H[k-2],H[k-1],v[i])<=0) k--;
//         H[k++]=v[i];
//     }
//     for(int i=n-2, t=k+1;i>=0;i--){
//         while(k>=t && cross(H[k-2],H[k-1],v[i])<=0) k--;
//         H[k++]=v[i];
//     }
//     return (size_t)(k-1);
// }

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        auto start = chrono::high_resolution_clock::now();
        
        auto start_load = chrono::high_resolution_clock::now();
        vector<pt> points = load_points_from_file("points.csv");
        auto end_load = chrono::high_resolution_clock::now();
        cout << "Loaded " << points.size() << " points\n";
        cout << "File loading time: " << chrono::duration<double>(end_load - start_load).count() << " seconds\n";
        
        std::cerr << "OMP threads: " << omp_get_max_threads()
          << " BASE_N=" << BASE_N
          << " TASK_CUTOFF=" << TASK_CUTOFF << "\n";
        
        cout << "--------------------------------\n";
        auto start_dt = chrono::high_resolution_clock::now();
        auto triangles = delaunay(points);
        cout << g_task_cnt.load() << " tasks will be created\n";
        cout << g_big_spawns.load() << " big tasks will be created\n";
        cout << g_mid_spawns.load() << " mid-size tasks will be created\n";
        auto start_convert = chrono::high_resolution_clock::now();
        auto tris_i = triangles_to_indices(triangles, points);
        auto end_convert = chrono::high_resolution_clock::now();
        cout << "triangles_to_indices time: " << chrono::duration<double>(end_convert - start_convert).count() << " seconds\n";
        
        auto end_dt = chrono::high_resolution_clock::now();
        cout << "Number of triangles: " << triangles.size() << "\n";
        cout << "delaunay execution time: " << chrono::duration<double>(end_dt - start_dt).count() << " seconds\n";
        
        auto start_rng = chrono::high_resolution_clock::now();
        size_t rng_count = 0;
            cudaProfilerStart();  // Start profiling
        auto vec = extract_rng_cuda_vec_idx(tris_i, points);
            cudaProfilerStop();   // Stop profiling
        rng_count = vec.size();
        
        auto end_rng = chrono::high_resolution_clock::now();
        cout << "Number of RNG edges: " << rng_count << "\n";
        cout << "extract execution time: " << chrono::duration<double>(end_rng - start_rng).count() << " seconds\n";

        auto end = chrono::high_resolution_clock::now();
        cout << "Total execution time: " << chrono::duration<double>(end - start).count() << " seconds\n";
        cout << "--------------------------------\n";
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
