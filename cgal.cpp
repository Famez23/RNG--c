#include <bits/stdc++.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using P   = K::Point_2;

using Vb  = CGAL::Triangulation_vertex_base_with_info_2<int, K>;
using Fb  = CGAL::Triangulation_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using DT  = CGAL::Delaunay_triangulation_2<K, Tds>;

using Traits = CGAL::Search_traits_2<K>;
using Tree   = CGAL::Kd_tree<Traits>;
using Sphere = CGAL::Fuzzy_sphere<Traits>;

struct I64PairHash {
    size_t operator()(const std::pair<long long,long long>& a) const noexcept {
        uint64_t x = (uint64_t)a.first, y = (uint64_t)a.second;
        x ^= y + 0x9e3779b97f4a7c15ull + (x<<6) + (x>>2);
        return (size_t)x;
    }
};

vector<P> load_points_csv_unique(const string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Failed to open file: " + filename);

    std::unordered_set<pair<long long,long long>, I64PairHash> uniq;
    uniq.reserve(1<<20);

    string line; size_t bad = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        long long x, y; char comma;
        if ((ss >> x >> comma >> y) && comma == ',') {
            uniq.emplace(x,y);
        } else {
            ++bad;
        }
    }
    if (bad) cerr << "Warning: skipped " << bad << " malformed line(s)\n";

    vector<P> pts; pts.reserve(uniq.size());
    for (const auto& q : uniq) pts.emplace_back((double)q.first, (double)q.second);
    return pts;
}

static inline bool is_rng_edge(const P& p1, const P& p2, const Tree& tree) {
    K::FT r2 = CGAL::squared_distance(p1, p2);     
    const double r = std::sqrt(CGAL::to_double(r2));

    vector<P> near_p1;
    near_p1.reserve(16);
    Sphere s1(p1, r, 0.0);
    tree.search(std::back_inserter(near_p1), s1);

    for (const P& q : near_p1) {
        if (q == p1 || q == p2) continue;
        if (CGAL::squared_distance(p1, q) < r2 &&
            CGAL::squared_distance(p2, q) < r2) {
            return false; 
        }
    }
    return true;
}

int main() {
    try {
        auto t_start1 = std::chrono::high_resolution_clock::now();
        vector<P> pts = load_points_csv_unique("points.csv");
        cout << "Loaded " << pts.size() << " unique points\n";
        if (pts.size() < 2) {
            cerr << "Need at least 2 points.\n";
            return 0;
        }
        auto t_start2 = std::chrono::high_resolution_clock::now();
        cout << "Point loading time: "
             << std::chrono::duration<double>(t_start2 - t_start1).count() << " s\n";
        auto t_start = std::chrono::high_resolution_clock::now();
        vector<std::pair<P,int>> with_info;
        with_info.reserve(pts.size());
        for (int i = 0; i < (int)pts.size(); ++i) with_info.emplace_back(pts[i], i);

        DT dt;
        dt.insert(with_info.begin(), with_info.end());
        cout << "Delaunay: " << dt.number_of_vertices() << " vertices, "
             << dt.number_of_faces() << " faces\n";
        auto t_delaunay = std::chrono::high_resolution_clock::now();
        cout << "Delaunay construction time: "
             << std::chrono::duration<double>(t_delaunay - t_start).count() << " s\n";
        Tree tree(pts.begin(), pts.end());

        vector<std::pair<int,int>> rng_edges;
        rng_edges.reserve(dt.number_of_vertices()); 

        for (auto e = dt.finite_edges_begin(); e != dt.finite_edges_end(); ++e) {
            auto fh = e->first;
            int  i  = e->second;
            auto v1 = fh->vertex((i+1)%3);
            auto v2 = fh->vertex((i+2)%3);
            const P& p1 = v1->point();
            const P& p2 = v2->point();
            int id1 = v1->info();
            int id2 = v2->info();
            if (id1 == id2) continue;

            if (is_rng_edge(p1, p2, tree)) {
                if (id2 < id1) std::swap(id1, id2);
                rng_edges.emplace_back(id1, id2);
            }
        }

        std::sort(rng_edges.begin(), rng_edges.end());
        rng_edges.erase(std::unique(rng_edges.begin(), rng_edges.end()), rng_edges.end());
        auto t_end = std::chrono::high_resolution_clock::now();
        cout << "RNG extraction time: "
             << std::chrono::duration<double>(t_end - t_delaunay).count() << " s\n";
        cout << "Total time: "
             << std::chrono::duration<double>(t_end - t_start1).count() << " s\n";

             cout << "RNG edges: " << rng_edges.size() << "\n";
        cout << "Number of Triangles: " << dt.number_of_faces() << "\n";

    } catch (const std::exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
