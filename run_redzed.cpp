#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <string>
#include <map>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <limits>
#include <memory>
#include <stdexcept>

using ValueT = double;
using IndexT = int64_t;

// ============================================================
// Distance Matrix Representation
// ============================================================
struct DistanceMatrix {
    int n;
    std::vector<ValueT> mat;
    DistanceMatrix(int n, const std::vector<std::vector<ValueT>>& m) : n(n), mat(n * n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                mat[i * n + j] = m[i][j];
            }
        }
    }
    inline ValueT get(int i, int j) const {
        return mat[i * n + j];
    }
};

// ============================================================
// Binomial Table / Combinadics Lookups
// ============================================================
struct BinomialTable {
    std::vector<IndexT> data;
    int nmax;
    int kmax;
    BinomialTable(int nmax, int kmax) : nmax(nmax), kmax(kmax), data((kmax + 1) * (nmax + 1), 0) {
        for (int n = 0; n <= nmax; ++n) {
            data[0 * (nmax + 1) + n] = 1;
            for (int k = 1; k <= std::min(n, kmax); ++k) {
                if (k == n) {
                    data[k * (nmax + 1) + n] = 1;
                } else {
                    data[k * (nmax + 1) + n] = data[(k - 1) * (nmax + 1) + (n - 1)] + data[k * (nmax + 1) + (n - 1)];
                }
            }
        }
    }
    inline IndexT choose(int n, int k) const {
        if (k > n || k < 0 || n < 0) return 0;
        return data[k * (nmax + 1) + n];
    }
};

int get_max_vertex(const BinomialTable& binom, IndexT idx, int k, int top_n) {
    int low = k - 1;
    int high = top_n;
    while (low < high) {
        int mid = (low + high + 1) >> 1;
        if (binom.choose(mid, k) <= idx) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }
    return low;
}

void decode_simplex(std::vector<int>& out, IndexT idx, int dim, int n, const BinomialTable& binom) {
    IndexT x = idx;
    int top_n = n - 1;
    for (int k = dim + 1; k >= 2; --k) {
        int v = get_max_vertex(binom, x, k, top_n);
        out[k - 1] = v;
        x -= binom.choose(v, k);
        top_n = v - 1;
    }
    out[0] = (int)x;
}

// ============================================================
// Boundary Enumeration
// ============================================================
struct BoundaryEnumerator {
    IndexT below;
    IndexT above;
    int j;
    int k;
    IndexT idx;
    int dim;

    BoundaryEnumerator() : below(0), above(0), j(0), k(0), idx(0), dim(0) {}

    void set_simplex(IndexT index, int dimension, int n) {
        below = index;
        above = 0;
        j = n - 1;
        k = dimension;
        idx = index;
        dim = dimension;
    }

    inline bool has_next_face() const { return k >= 0; }

    IndexT next_face_index(const BinomialTable& binom) {
        int j_val = get_max_vertex(binom, below, k + 1, j);
        IndexT c1 = binom.choose(j_val, k + 1);
        IndexT face_index = above - c1 + below;
        j = j_val;
        below -= c1;
        above += binom.choose(j_val, k);
        k -= 1;
        return face_index;
    }
};

// ============================================================
// Filtration Enumeration (Heap & Sibling Streams)
// ============================================================
struct Simplex {
    int dim;
    IndexT idx;
    ValueT diam;
};

inline bool simplex_less(const Simplex& a, const Simplex& b) {
    if (a.diam != b.diam) return a.diam < b.diam;
    if (a.dim != b.dim) return a.dim < b.dim;
    return a.idx < b.idx;
}

struct SiblingStream {
    std::vector<Simplex> items;
    size_t pos;
    SiblingStream(std::vector<Simplex>&& items) : items(std::move(items)), pos(0) {}
    inline bool is_empty() const { return pos >= items.size(); }
    inline const Simplex& head() const { return items[pos]; }
    inline void advance() { pos++; }
};

struct HeadEntry {
    Simplex simplex;
    int id;
};

struct HeadHeap {
    std::vector<HeadEntry> data;
    inline bool is_empty() const { return data.empty(); }
    void push_head(const HeadEntry& x) {
        data.push_back(x);
        std::push_heap(data.begin(), data.end(), [](const HeadEntry& a, const HeadEntry& b) {
            return simplex_less(b.simplex, a.simplex);
        });
    }
    HeadEntry pop_head() {
        std::pop_heap(data.begin(), data.end(), [](const HeadEntry& a, const HeadEntry& b) {
            return simplex_less(b.simplex, a.simplex);
        });
        HeadEntry x = data.back();
        data.pop_back();
        return x;
    }
};

inline ValueT child_diameter(const std::vector<int>& verts, int parent_dim, ValueT parent_diam, int w, const DistanceMatrix& dist) {
    ValueT diam = parent_diam;
    int m = parent_dim + 1;
    for (int i = 0; i < m; ++i) {
        ValueT d = dist.get(w, verts[i]);
        if (d > diam) diam = d;
    }
    return diam;
}

struct FiltrationEnumerator {
    DistanceMatrix dist;
    int n;
    int max_dim;
    ValueT threshold;
    BinomialTable binom;
    std::vector<int> verts;
    std::vector<SiblingStream> streams;
    HeadHeap heap;

    FiltrationEnumerator(const DistanceMatrix& dist, int max_dim, ValueT threshold)
        : dist(dist), n(dist.n), max_dim(max_dim), threshold(threshold), binom(dist.n, max_dim + 2), verts(max_dim + 1) {
        std::vector<Simplex> roots(n);
        for (int v = 0; v < n; ++v) roots[v] = {0, v, 0.0};
        streams.emplace_back(std::move(roots));
        if (!streams[0].is_empty()) heap.push_head({streams[0].head(), 0});
    }
};

std::unique_ptr<SiblingStream> child_stream(FiltrationEnumerator& fe, const Simplex& s) {
    if (s.dim >= fe.max_dim) return nullptr;
    decode_simplex(fe.verts, s.idx, s.dim, fe.n, fe.binom);
    int m = s.dim + 1;
    int lastv = fe.verts[m - 1];
    std::vector<Simplex> children;
    for (int w = lastv + 1; w < fe.n; ++w) {
        ValueT diam = child_diameter(fe.verts, s.dim, s.diam, w, fe.dist);
        if (diam <= fe.threshold) {
            IndexT child_idx = s.idx + fe.binom.choose(w, m + 1);
            children.push_back({s.dim + 1, child_idx, diam});
        }
    }
    if (children.empty()) return nullptr;
    std::sort(children.begin(), children.end(), simplex_less);
    return std::make_unique<SiblingStream>(std::move(children));
}

bool pop_simplex(FiltrationEnumerator& fe, Simplex& out_s) {
    if (fe.heap.is_empty()) return false;
    HeadEntry head_entry = fe.heap.pop_head();
    out_s = head_entry.simplex;
    int sid = head_entry.id;
    fe.streams[sid].advance();
    if (!fe.streams[sid].is_empty()) fe.heap.push_head({fe.streams[sid].head(), sid});
    auto c_stream = child_stream(fe, out_s);
    if (c_stream) {
        fe.streams.push_back(std::move(*c_stream));
        int new_sid = fe.streams.size() - 1;
        if (!fe.streams[new_sid].is_empty()) fe.heap.push_head({fe.streams[new_sid].head(), new_sid});
    }
    return true;
}

// ============================================================
// Dictionary and Reduction Set Math
// ============================================================
inline void toggle(boost::unordered_flat_set<int>& s, int x) {
    auto it = s.find(x);
    if (it != s.end()) s.erase(it); else s.insert(x);
}
inline void xor_add(boost::unordered_flat_set<int>& dst, const boost::unordered_flat_set<int>& src) {
    for (int x : src) {
	auto [it, inserted] = dst.insert(x);
	if (!inserted) dst.erase(it);
    }
}
inline bool xor_from(boost::unordered_flat_set<int>& dst, const boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, int x) {
    auto it = R.find(x);
    if (it == R.end()) return false;
    xor_add(dst, it->second);
    return true;
}
inline void insert_live(boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, int x) {
    R[x] = {x}; Ri[x] = {x};
}
inline bool is_live(const boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, int x) {
    return R.find(x) != R.end();
}

void remove_maximal(boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, int y, std::vector<int>& removed) {
    removed.clear();
    auto it_users = Ri.find(y);
    if (it_users == Ri.end()) {
        auto it_set_y = R.find(y);
        if (it_set_y != R.end() && it_set_y->second.size() == 1 && it_set_y->second.count(y)) {
            R.erase(it_set_y); removed.push_back(y);
        }
        return;
    }
    std::vector<int> touched(it_users->second.begin(), it_users->second.end());
    Ri.erase(it_users);
    for (int z : touched) {
        auto it_set_z = R.find(z);
        if (it_set_z == R.end()) continue;
        it_set_z->second.erase(y);
        if (it_set_z->second.empty()) {
            R.erase(it_set_z); removed.push_back(z);
        }
    }
}

void replace_pivot(boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, const boost::unordered_flat_set<int>& bar, int j, std::vector<int>& removed) {
    removed.clear();
    auto it_users = Ri.find(j);
    if (it_users == Ri.end()) return;
    std::vector<int> touched(it_users->second.begin(), it_users->second.end());
    std::vector<int> others;
    others.reserve(bar.size());
    for (int y : bar) if (y != j) others.push_back(y);
    for (int z : touched) {
        auto it_set_z = R.find(z);
        if (it_set_z == R.end()) continue;
        xor_add(it_set_z->second, bar);
        for (int y : others) {
            toggle(Ri[y], z);
            if (Ri[y].empty()) Ri.erase(y);
        }
        if (it_set_z->second.empty()) {
            R.erase(it_set_z); removed.push_back(z);
        }
    }
    Ri.erase(j);
}

inline ValueT compute_diameter(std::vector<int>& buf, IndexT idx, int dim, const DistanceMatrix& dist, int n, const BinomialTable& binom) {
    decode_simplex(buf, idx, dim, n, binom);
    int m = dim + 1;
    ValueT diam = 0.0;
    for (int a = 0; a < m; ++a) {
        for (int b = a + 1; b < m; ++b) {
            ValueT d = dist.get(buf[a], buf[b]);
            if (d > diam) diam = d;
        }
    }
    return diam;
}

std::pair<int, ValueT> pick_pivot(const boost::unordered_flat_set<int>& bar, int dim, const DistanceMatrix& dist, int nverts, const BinomialTable& binom, std::vector<int>& buf) {
    bool first_seen = true; int best_idx = 0; ValueT best_diam = 0.0;
    for (int j : bar) {
        ValueT dj = compute_diameter(buf, j, dim, dist, nverts, binom);
        if (first_seen || (dj > best_diam) || ((dj == best_diam) && (j > best_idx))) {
            best_idx = j; best_diam = dj; first_seen = false;
        }
    }
    return {best_idx, best_diam};
}

// ============================================================
// Active Enumeration Logic (Combinadic Simplices Parsing)
// ============================================================
inline IndexT face_index_without(const std::vector<int>& verts, int m, int omit, const BinomialTable& binom) {
    IndexT idx = 0;
    for (int k = 0; k < omit; ++k) idx += binom.choose(verts[k], k + 1);
    for (int k = omit + 1; k < m; ++k) idx += binom.choose(verts[k], k);
    return idx;
}
inline bool contains_vertex(const std::vector<int>& verts, int m, int w) {
    for (int i = 0; i < m; ++i) if (verts[i] == w) return true;
    return false;
}
void remove_face_entry(std::vector<std::pair<int, int>>& lst, int idx) {
    for (size_t k = 0; k < lst.size(); ++k) {
        if (lst[k].first == idx) {
            lst[k] = lst.back(); lst.pop_back(); return;
        }
    }
}
inline int merge_vertex(std::vector<int>& out, const std::vector<int>& xverts, int m, int w) {
    int i = 0, j = 0;
    while (i < m && xverts[i] < w) { out[j] = xverts[i]; i++; j++; }
    out[j] = w; int wpos = j; j++;
    while (i < m) { out[j] = xverts[i]; i++; j++; }
    return wpos;
}
inline ValueT face_diameter(const std::vector<int>& verts, int m, int omit, const DistanceMatrix& dist) {
    ValueT diam = 0.0;
    for (int a = 0; a < m; ++a) {
        if (a == omit) continue;
        for (int b = a + 1; b < m; ++b) {
            if (b == omit) continue;
            ValueT d = dist.get(verts[a], verts[b]);
            if (d > diam) diam = d;
        }
    }
    return diam;
}
bool coface_ready(std::vector<int>& candverts, const std::vector<int>& xverts, int w, int xidx, ValueT xdiam, const DistanceMatrix& dist, const BinomialTable& binom) {
    int m = (int)xverts.size(); int total = m + 1;
    int wpos = merge_vertex(candverts, xverts, m, w);
    for (int omit = 0; omit < total; ++omit) {
        if (omit == wpos) continue;
        IndexT faceidx = face_index_without(candverts, total, omit, binom);
        if (faceidx < xidx) continue;
        ValueT facedia = face_diameter(candverts, total, omit, dist);
        if (facedia >= xdiam) return false;
    }
    return true;
}

// Cache Infrastructure
struct NeighborCache { std::vector<std::vector<int>> order; std::vector<int> limit; };
NeighborCache build_neighbor_cache(const DistanceMatrix& dist, int nverts) {
    std::vector<std::vector<int>> order(nverts); std::vector<int> limit(nverts, 0);
    int nn = nverts - 1;
    std::vector<int> ids(std::max(nn, 0)); std::vector<ValueT> ds(std::max(nn, 0)); std::vector<int> perm(std::max(nn, 0));
    for (int v = 0; v < nverts; ++v) {
        std::vector<int> row(nn); int t = 0;
        for (int w = 0; w < nverts; ++w) {
            if (w == v) continue; ids[t] = w; ds[t] = dist.get(v, w); perm[t] = t; t++;
        }
        std::sort(perm.begin(), perm.end(), [&](int a, int b) { return ds[a] < ds[b]; });
        for (int i = 0; i < nn; ++i) row[i] = ids[perm[i]];
        order[v] = std::move(row);
    }
    return {std::move(order), std::move(limit)};
}
inline int advance_limit(NeighborCache& cache, int v, ValueT r, const DistanceMatrix& dist) {
    const auto& row = cache.order[v]; int limit = cache.limit[v]; int L = (int)row.size();
    while (limit < L && dist.get(v, row[limit]) <= r) limit++;
    cache.limit[v] = limit; return limit;
}
inline std::pair<int, int> choose_pivot_vertex(const std::vector<int>& xverts, int m, ValueT r, const DistanceMatrix& dist, NeighborCache& cache) {
    int best_pos = 0; int best_lim = std::numeric_limits<int>::max();
    for (int i = 0; i < m; ++i) {
        int limit = advance_limit(cache, xverts[i], r, dist);
        if (limit < best_lim) { best_lim = limit; best_pos = i; }
    }
    return {best_pos, best_lim};
}
bool has_leftover_coface(std::vector<bool>& seen_w, const std::vector<int>& xverts, int xidx, ValueT r, const DistanceMatrix& dist, const BinomialTable& binom, NeighborCache& neighbors, std::vector<int>& candverts) {
    int m = (int)xverts.size();
    auto pivot_info = choose_pivot_vertex(xverts, m, r, dist, neighbors);
    int pivot_pos = pivot_info.first, limit = pivot_info.second, pivot = xverts[pivot_pos];
    const auto& row = neighbors.order[pivot];
    for (int t = 0; t < limit; ++t) {
        int w = row[t]; if (contains_vertex(xverts, m, w)) continue; if (seen_w[w]) continue;
        bool ok = true;
        for (int i = 0; i < m; ++i) {
            if (i == pivot_pos) continue;
            if (dist.get(w, xverts[i]) > r) { ok = false; break; }
        }
        if (!ok) continue;
        if (coface_ready(candverts, xverts, w, xidx, r, dist, binom)) return true;
    }
    return false;
}

void register_active(boost::unordered_flat_map<int, std::vector<int>>& active_verts, boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>>& face_to_active, int idx, const std::vector<int>& verts, int n, const BinomialTable& binom) {
    active_verts[idx] = verts; int m = n + 1;
    for (int omit = 0; omit < m; ++omit) {
        IndexT faceidx = face_index_without(verts, m, omit, binom);
        face_to_active[faceidx].push_back({idx, verts[omit]});
    }
}
void unregister_active(boost::unordered_flat_map<int, std::vector<int>>& active_verts, boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>>& face_to_active, const std::vector<int>& doomed, int n, const BinomialTable& binom) {
    int m = n + 1;
    for (int idx : doomed) {
        auto it = active_verts.find(idx); if (it == active_verts.end()) continue;
        const auto& verts = it->second;
        for (int omit = 0; omit < m; ++omit) {
            IndexT faceidx = face_index_without(verts, m, omit, binom);
            auto it_face = face_to_active.find(faceidx);
            if (it_face != face_to_active.end()) {
                remove_face_entry(it_face->second, idx);
                if (it_face->second.empty()) face_to_active.erase(it_face);
            }
        }
        active_verts.erase(it);
    }
}

struct SeedState {
    std::vector<bool> seen; std::vector<int> verts; std::vector<std::vector<int>> seeds;
    SeedState(int nverts) : seen(nverts, false), seeds(nverts) {}
    void reset_seeds() { for (int w : verts) { seen[w] = false; seeds[w].clear(); } verts.clear(); }
};

void collect_candidates(SeedState& state, const std::vector<int>& xverts, ValueT r, const DistanceMatrix& dist, int n, const BinomialTable& binom, const boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>>& face_to_active) {
    state.reset_seeds(); int m = n + 1;
    for (int omit = 0; omit < m; ++omit) {
        IndexT faceidx = face_index_without(xverts, m, omit, binom);
        auto it = face_to_active.find(faceidx); if (it == face_to_active.end()) continue;
        int xextra = xverts[omit];
        for (const auto& pair : it->second) {
            int yidx = pair.first, w = pair.second;
            if (dist.get(xextra, w) <= r) {
                if (!state.seen[w]) { state.seen[w] = true; state.verts.push_back(w); }
                state.seeds[w].push_back(yidx);
            }
        }
    }
}

void process_seeded(boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, std::vector<std::pair<ValueT, ValueT>>& top_pairs, int xidx, const std::vector<int>& other_seed_faces, ValueT top_diam, FiltrationEnumerator& fe, boost::unordered_flat_set<int>& top_bar, std::vector<int>& diam_buf, std::vector<int>& removed, boost::unordered_flat_map<int, std::vector<int>>& active_verts, boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>>& face_to_active, int n) {
    top_bar.clear(); xor_from(top_bar, R, xidx);
    for (int yidx : other_seed_faces) xor_from(top_bar, R, yidx);
    if (top_bar.empty()) return;
    if (top_bar.size() == 1) {
        int j = *top_bar.begin(); remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, fe.binom);
        ValueT jdiam = compute_diameter(diam_buf, j, n, fe.dist, fe.n, fe.binom);
        if (top_diam != jdiam) top_pairs.push_back({jdiam, top_diam});
        return;
    }
    auto pivot_info = pick_pivot(top_bar, n, fe.dist, fe.n, fe.binom, diam_buf);
    int j = pivot_info.first; ValueT jdiam = pivot_info.second;
    replace_pivot(R, Ri, top_bar, j, removed);
    if (!removed.empty()) unregister_active(active_verts, face_to_active, removed, n, fe.binom);
    if (top_diam != jdiam) top_pairs.push_back({jdiam, top_diam});
}

void process_leftover(boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& R, boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, std::vector<std::pair<ValueT, ValueT>>& top_pairs, int xidx, ValueT top_diam, FiltrationEnumerator& fe, std::vector<int>& diam_buf, std::vector<int>& removed, boost::unordered_flat_map<int, std::vector<int>>& active_verts, boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>>& face_to_active, int n) {
    if (!is_live(R, xidx)) return;
    if (is_live(Ri, xidx)) {
        int j = xidx; remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, fe.binom);
        ValueT jdiam = compute_diameter(diam_buf, j, n, fe.dist, fe.n, fe.binom);
        if (top_diam != jdiam) top_pairs.push_back({jdiam, top_diam});
        return;
    }
    boost::unordered_flat_set<int> top_bar = R[xidx]; if (top_bar.empty()) return;
    if (top_bar.size() == 1) {
        int j = *top_bar.begin(); remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, fe.binom);
        ValueT jdiam = compute_diameter(diam_buf, j, n, fe.dist, fe.n, fe.binom);
        if (top_diam != jdiam) top_pairs.push_back({jdiam, top_diam});
        return;
    }
    auto pivot_info = pick_pivot(top_bar, n, fe.dist, fe.n, fe.binom, diam_buf);
    int j = pivot_info.first; ValueT jdiam = pivot_info.second;
    replace_pivot(R, Ri, top_bar, j, removed);
    if (!removed.empty()) unregister_active(active_verts, face_to_active, removed, n, fe.binom);
    if (top_diam != jdiam) top_pairs.push_back({jdiam, top_diam});
}

void append_infinite(std::vector<std::pair<ValueT, ValueT>>& out, const boost::unordered_flat_map<int, boost::unordered_flat_set<int>>& Ri, int dim, const DistanceMatrix& dist, int nverts, const BinomialTable& binom, std::vector<int>& diam_buf) {
    std::vector<int> live_maxima; for (const auto& pair : Ri) live_maxima.push_back(pair.first);
    std::sort(live_maxima.begin(), live_maxima.end());
    ValueT infv = std::numeric_limits<ValueT>::infinity();
    for (int x : live_maxima) {
        ValueT birth = compute_diameter(diam_buf, x, dim, dist, nverts, binom);
        out.push_back({birth, infv});
    }
}

// ============================================================
// Core Execution Block
// ============================================================
std::map<int, std::vector<std::pair<ValueT, ValueT>>> redzed_pipeline(const std::vector<std::vector<ValueT>>& mat_in, int n, ValueT threshold = std::numeric_limits<ValueT>::infinity()) {
    if (n < 1) throw std::invalid_argument("requires n >= 1");
    int nverts = mat_in.size();
    ValueT cone_val = std::numeric_limits<ValueT>::infinity();
    for (int j = 0; j < nverts; ++j) {
        ValueT max_col = 0.0;
        for (int i = 0; i < nverts; ++i) if (mat_in[i][j] > max_col) max_col = mat_in[i][j];
        if (max_col < cone_val) cone_val = max_col;
    }
    if (cone_val < threshold) threshold = cone_val;

    DistanceMatrix dist(nverts, mat_in);
    FiltrationEnumerator fe(dist, n, threshold); BoundaryEnumerator be;
    std::vector<boost::unordered_flat_map<int, boost::unordered_flat_set<int>>> R(n + 1), Ri(n + 1);
    std::vector<std::vector<std::pair<ValueT, ValueT>>> pairs(n + 1);
    boost::unordered_flat_set<int> bar, top_bar;
    std::vector<std::vector<int>> diam_bufs(n + 1);
    for (int d = 0; d <= n; ++d) diam_bufs[d].resize(d + 1);
    std::vector<int> top_buf(n + 1), coface_buf(n + 2), removed;
    boost::unordered_flat_map<int, std::vector<int>> active_verts;
    boost::unordered_flat_map<IndexT, std::vector<std::pair<int, int>>> face_to_active;
    SeedState seeds(fe.n); NeighborCache neighbors = build_neighbor_cache(fe.dist, fe.n);

    Simplex s;
    while (pop_simplex(fe, s)) {
        if (s.dim == 0) { insert_live(R[0], Ri[0], (int)s.idx); continue; }
        bar.clear(); be.set_simplex(s.idx, s.dim, fe.n);
        while (be.has_next_face()) {
            int face = (int)be.next_face_index(fe.binom);
            xor_from(bar, R[s.dim - 1], face);
        }
        if (bar.empty()) {
            int xidx = (int)s.idx; insert_live(R[s.dim], Ri[s.dim], xidx);
            if (s.dim == n) {
                decode_simplex(top_buf, s.idx, n, fe.n, fe.binom);
                collect_candidates(seeds, top_buf, s.diam, fe.dist, n, fe.binom, face_to_active);
                for (int w : seeds.verts) {
                    if (coface_ready(coface_buf, top_buf, w, xidx, s.diam, fe.dist, fe.binom)) {
                        process_seeded(R[n], Ri[n], pairs[n], xidx, seeds.seeds[w], s.diam, fe, top_bar, diam_bufs[n], removed, active_verts, face_to_active, n);
                    }
                }
                if (is_live(R[n], xidx)) {
                    if (has_leftover_coface(seeds.seen, top_buf, xidx, s.diam, fe.dist, fe.binom, neighbors, coface_buf)) {
                        process_leftover(R[n], Ri[n], pairs[n], xidx, s.diam, fe, diam_bufs[n], removed, active_verts, face_to_active, n);
                    }
                }
                if (is_live(R[n], xidx)) register_active(active_verts, face_to_active, xidx, top_buf, n, fe.binom);
            }
        } else if (bar.size() == 1) {
            int j = *bar.begin(); remove_maximal(R[s.dim - 1], Ri[s.dim - 1], j, removed);
            ValueT jdiam = compute_diameter(diam_bufs[s.dim - 1], j, s.dim - 1, fe.dist, fe.n, fe.binom);
            if (s.diam != jdiam) pairs[s.dim - 1].push_back({jdiam, s.diam});
        } else {
            auto pivot_info = pick_pivot(bar, s.dim - 1, fe.dist, fe.n, fe.binom, diam_bufs[s.dim - 1]);
            int j = pivot_info.first; ValueT jdiam = pivot_info.second;
            replace_pivot(R[s.dim - 1], Ri[s.dim - 1], bar, j, removed);
            if (s.diam != jdiam) pairs[s.dim - 1].push_back({jdiam, s.diam});
        }
    }
    for (int p = 0; p <= n; ++p) append_infinite(pairs[p], Ri[p], p, fe.dist, fe.n, fe.binom, diam_bufs[p]);
    std::map<int, std::vector<std::pair<ValueT, ValueT>>> res;
    for (int d = 0; d <= n; ++d) res[d] = pairs[d];
    return res;
}

// ============================================================
// Robust Input Parsers
// ============================================================
std::vector<std::vector<ValueT>> load_point_cloud(const std::string& filename) {
    std::ifstream file(filename); std::string line; std::vector<std::vector<ValueT>> points;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::replace(line.begin(), line.end(), ',', ' '); std::stringstream ss(line);
        ValueT coord; std::vector<ValueT> point;
        while (ss >> coord) point.push_back(coord);
        if (!point.empty()) points.push_back(point);
    }
    size_t n = points.size(); std::vector<std::vector<ValueT>> dist_mat(n, std::vector<ValueT>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            ValueT sum_sq = 0.0;
            for (size_t d = 0; d < points[i].size(); ++d) { ValueT diff = points[i][d] - points[j][d]; sum_sq += diff * diff; }
            dist_mat[i][j] = std::sqrt(sum_sq);
        }
    }
    return dist_mat;
}

std::vector<std::vector<ValueT>> load_lower_distance_matrix(const std::string& filename) {
    std::ifstream file(filename); std::string line; std::vector<ValueT> values;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::replace(line.begin(), line.end(), ',', ' '); std::stringstream ss(line);
        ValueT val; while (ss >> val) values.push_back(val);
    }
    size_t M = values.size(); size_t discriminant = 1 + 8 * M;
    size_t sqrt_disc = std::round(std::sqrt(discriminant));
    if (sqrt_disc * sqrt_disc != discriminant) throw std::runtime_error("Invalid strictly lower distance shape.");
    size_t n = (1 + sqrt_disc) / 2;
    std::vector<std::vector<ValueT>> dist_mat(n, std::vector<ValueT>(n, 0.0)); size_t idx = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) { dist_mat[i][j] = values[idx]; dist_mat[j][i] = values[idx]; idx++; }
    }
    return dist_mat;
}

std::vector<std::vector<ValueT>> load_full_matrix(const std::string& filename) {
    std::ifstream file(filename); std::string line; std::vector<std::vector<ValueT>> mat;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::replace(line.begin(), line.end(), ',', ' '); std::stringstream ss(line);
        ValueT val; std::vector<ValueT> row;
        while (ss >> val) row.push_back(val);
        if (!row.empty()) mat.push_back(row);
    }
    return mat;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <format: pc|txt|full> <max_dim> <threshold> <input_file>\n";
        return 1;
    }
    std::string format = argv[1]; int max_dim = std::stoi(argv[2]);
    ValueT threshold = (std::string(argv[3]) == "inf") ? std::numeric_limits<ValueT>::infinity() : std::stod(argv[3]);
    std::string input_file = argv[4];

    std::vector<std::vector<ValueT>> dist_matrix;
    if (format == "pc") dist_matrix = load_point_cloud(input_file);
    else if (format == "txt") dist_matrix = load_lower_distance_matrix(input_file);
    else if (format == "full") dist_matrix = load_full_matrix(input_file);
    else { std::cerr << "Unknown format profile.\n"; return 1; }

    std::cout << "Dataset loaded. Matrix size: " << dist_matrix.size() << "x" << dist_matrix.size() << "\n";

    auto start = std::chrono::high_resolution_clock::now();
    auto results = redzed_pipeline(dist_matrix, max_dim, threshold);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Execution Time (C++ Redzed): " << elapsed.count() << " ms\n";

    std::ofstream f("output_cpp.txt");
    f << std::setprecision(16);
    for (int d = 0; d <= max_dim; ++d) {
        f << "persistence intervals in dimension " << d << ":\n";
        auto sorted_pairs = results[d];
        std::sort(sorted_pairs.begin(), sorted_pairs.end(), [](const std::pair<ValueT, ValueT>& a, const std::pair<ValueT, ValueT>& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
        for (const auto& pair : sorted_pairs) {
            f << " [" << pair.first << "," << pair.second << ")\n";
        }
    }
    std::cout << "Results written to output_cpp.txt\n";
    return 0;
}
