#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 30'500'000;
    std::mt19937_64 rng(12345); 
    std::uniform_int_distribution<long long> dist(0, 200'000'000); 

    ofstream out("points.csv");
    if (!out) {
        cerr << "Failed to open output file\n";
        return 1;
    }

    unordered_set<uint64_t> seen;
    auto mix64 = [](uint64_t x,uint64_t y){
        x ^= y + 0x9e3779b97f4a7c15ull + (x<<6)+(x>>2);
        return x;
    };

    int written = 0;
    while (written < N) {
        long long x = dist(rng);
        long long y = dist(rng);
        uint64_t key = mix64((uint64_t)x,(uint64_t)y);
        if (seen.insert(key).second) {
            out << x << "," << y << "\n";
            written++;
        }
    }

    cerr << "Wrote " << written << " unique points to points.csv\n";
    return 0;
}
