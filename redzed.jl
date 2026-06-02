#=
main function is redzed(dist, n; threshold)
    dist = distance matrix 
    n = max homology dimension 
    threshold = radius threshold 

    output is a dictionary of the form: dimension => pairs of that dimension
=#


const ValueT = Float64
const IndexT = Int

# ============================================================
# Distance matrix wrapper
# ============================================================

struct DistanceMatrix
    n::Int # number of rows/cols stored for fast access
    mat::Matrix{ValueT}
end

# Convert abstract matrix to dense distance matrix
function DistanceMatrix(mat::AbstractMatrix)
    n = size(mat, 1)
    @assert size(mat, 2) == n "distance matrix must be square"
    return DistanceMatrix(n, ValueT.(mat))
end

# ============================================================
# Binomial table / combinadics
# ============================================================

# Binomial table stores precomputed values to avoid re computing 
struct BinomialTable
    data::Matrix{IndexT}
    nmax::Int
    kmax::Int
end

# Compute the binomail table to support all (n choose k) up to (nmax choose kmax)
function BinomialTable(nmax::Int, kmax::Int)
    data = zeros(IndexT, kmax + 1, nmax + 1)

    @inbounds for n in 0:nmax
        data[1, n + 1] = 1
        for k in 1:min(n, kmax)
            if k == n
                data[k + 1, n + 1] = 1
            else
                data[k + 1, n + 1] = data[k, n] + data[k + 1, n] # Pascal's triangle: (n k) = (n-1 k-1) + (n-1 k)
            end
        end
    end

    return BinomialTable(data, nmax, kmax)
end

# returns (n choose k) from the binomial table lookup
@inline function choose(binom::BinomialTable, n::Int, k::Int)::IndexT
    if k > n
        return 0 # zero by convention
    end
    return @inbounds binom.data[k + 1, n + 1]
end

# for decodng simplices, gets largest n such that (n choose k) < idx
# top_n = maximal possible n 
function get_max_vertex(binom::BinomialTable, idx::IndexT, k::Int, top_n::Int)::Int
    low = k - 1
    high = top_n

    while low < high
        mid = (low + high + 1) >>> 1 # binary search 
        if choose(binom, mid, k) <= idx # if (n k) still smaller, upper half
            low = mid
        else
            high = mid - 1 # else lower half
        end
    end

    return low
end

# reconstruct vertex list from simplex idx and dim
# out = reusuable buffer to avoid extra allocs
function decode_simplex!(
    out::Vector{Int},
    idx::IndexT,
    dim::Int,
    n::Int,
    binom::BinomialTable,
)
    x = idx
    top_n = n - 1

    @inbounds for k in (dim + 1):-1:2
        v = get_max_vertex(binom, x, k, top_n) # decode k-th simplex 
        out[k] = v # store in buffer
        x -= choose(binom, v, k) # strip off (v k) from idx 
        top_n = v - 1 # next vertex will be lower idx
    end

    out[1] = x
    return out
end

# ============================================================
# Boundary enumerator
# ============================================================

# struct for boundary enumeration
mutable struct BoundaryEnumerator
    below::IndexT # below partial sum
    above::IndexT # above partial sum
    j::Int # upper search bound
    k::Int # binomial level to peel off
    idx::IndexT
    dim::Int
end

BoundaryEnumerator() = BoundaryEnumerator(0, 0, 0, 0, 0, 0)

# initialize be given new simplex
function set_simplex!(be::BoundaryEnumerator, idx::IndexT, dim::Int, n::Int)
    be.below = idx
    be.above = 0
    be.j = n - 1
    be.k = dim
    be.idx = idx
    be.dim = dim
    return be
end

# used for iterating 
@inline has_next_face(be::BoundaryEnumerator) = be.k >= 0

# decode then next face index and update enumerator 
function next_face_index!(be::BoundaryEnumerator, binom::BinomialTable)::IndexT
    j = get_max_vertex(binom, be.below, be.k + 1, be.j)
    c1 = choose(binom, j, be.k + 1)
    face_index = be.above - c1 + be.below

    be.j = j
    be.below -= c1
    be.above += choose(binom, j, be.k)
    be.k -= 1

    return face_index
end

# ============================================================
# Filtration enumeration
# ============================================================
# Overall filrtation logic:
# Heap with entries head of sibling streams
# Starts with verts in one stream with head as entry in heap 
# On a step, pop the root of heap, next simplex in filtration 
# Stream that the simplex came from is advanced, next head pushed into heap 
# Then generate children from popped simplex, add head into heap 

# Simplex stored as dimension, idx, diam
struct Simplex
    dim::Int
    idx::IndexT
    diam::ValueT
end

# Tests if a < b according to diam, dim, idx
@inline function simplex_less(
    a::Simplex,
    b::Simplex,
)::Bool
    return (a.diam < b.diam) ||
           ((a.diam == b.diam) && (a.dim < b.dim)) ||
           ((a.diam == b.diam) && (a.dim == b.dim) && (a.idx < b.idx))
end

# Given verts and parent diameter, determines the diameter of a simplex formed by adding w
@inline function child_diameter!(
    verts::Vector{Int},
    parent_dim::Int,
    parent_diam::ValueT,
    w::Int,
    dist::DistanceMatrix,
)::ValueT
    diam = parent_diam
    m = parent_dim + 1
    mat = dist.mat
    wp1 = w + 1

    @inbounds for i in 1:m
        d = mat[wp1, verts[i] + 1]
        if d > diam
            diam = d
        end
    end

    return diam
end

# List to store siblings (formed from same simplex)
mutable struct SiblingStream
    items::Vector{Simplex}
    pos::Int # current position
end

# Empty test
@inline function is_empty(st::SiblingStream)::Bool
    return st.pos > length(st.items)
end

# Get current item
@inline function head(st::SiblingStream)::Simplex
    return @inbounds st.items[st.pos]
end

# Increase position
@inline function advance!(st::SiblingStream)
    st.pos += 1
    return st
end

# Heap entry
# Current head simplex of stream + id of stream it came from
struct HeadEntry
    simplex::Simplex
    id::Int
end

# returns head of a < b
@inline function entry_less(a::HeadEntry, b::HeadEntry)::Bool
    return simplex_less(a.simplex, b.simplex)
end

# Heap structure
# Sorts each stream by head simplex
mutable struct HeadHeap
    data::Vector{HeadEntry}
end

HeadHeap() = HeadHeap(HeadEntry[])
Base.isempty(h::HeadHeap) = isempty(h.data)

# Adds x to heap and increases while less than parent
function push_head!(h::HeadHeap, x::HeadEntry)
    push!(h.data, x)
    i = length(h.data)

    @inbounds while i > 1
        p = i >>> 1
        if entry_less(h.data[i], h.data[p])
            h.data[i], h.data[p] = h.data[p], h.data[i]
            i = p
        else
            break
        end
    end

    return h
end

# pop head and re sort
function pop_head!(h::HeadHeap)::HeadEntry
    n = length(h.data)

    @inbounds h.data[1], h.data[n] = h.data[n], h.data[1]
    x = pop!(h.data)
    n -= 1

    i = 1
    @inbounds while true
        l = i << 1
        r = l + 1

        if l > n
            break
        end

        best = l
        if r <= n && entry_less(h.data[r], h.data[l])
            best = r
        end

        if entry_less(h.data[best], h.data[i])
            h.data[i], h.data[best] = h.data[best], h.data[i]
            i = best
        else
            break
        end
    end

    return x
end

# filtration enumerator struct
mutable struct FiltrationEnumerator
    dist::DistanceMatrix
    n::Int
    max_dim::Int
    threshold::ValueT
    binom::BinomialTable
    verts::Vector{Int} # reusable buffer
    streams::Vector{SiblingStream} # sibling streams created
    heap::HeadHeap
end

# Initialize the heap by adding all verts into a simbling stream
function FiltrationEnumerator(
    dist::DistanceMatrix;
    max_dim::Int = 1,
    threshold::Real = Inf,
)
    n = dist.n
    binom = BinomialTable(n, max_dim + 2)
    verts = Vector{Int}(undef, max_dim + 1)

    streams = SiblingStream[]
    heap = HeadHeap()

    roots = Vector{Simplex}(undef, n)
    @inbounds for v in 0:(n - 1)
        roots[v + 1] = Simplex(0, v, 0.0)
    end

    root_stream = SiblingStream(roots, 1)
    push!(streams, root_stream)

    if !is_empty(root_stream)
        push_head!(heap, HeadEntry(head(root_stream), 1))
    end

    return FiltrationEnumerator(
        dist,
        n,
        max_dim,
        ValueT(threshold),
        binom,
        verts,
        streams,
        heap,
    )
end

# Builds children of a simplex
function child_stream!(
    fe::FiltrationEnumerator,
    s::Simplex,
)
    if s.dim >= fe.max_dim # dont build n+1 streams since they're handled differently
        return nothing
    end

    decode_simplex!(fe.verts, s.idx, s.dim, fe.n, fe.binom)

    m = s.dim + 1
    lastv = fe.verts[m]

    children = Simplex[]
    sizehint!(children, max(0, fe.n - lastv - 1))

    @inbounds for w in (lastv + 1):(fe.n - 1) # w > lastv to ensure each simplex only created once
        diam = child_diameter!(fe.verts, s.dim, s.diam, w, fe.dist)
        if diam <= fe.threshold
            child_idx = s.idx + choose(fe.binom, w, m + 1) # idx of child simplex
            push!(children, Simplex(s.dim + 1, child_idx, diam))
        end
    end

    isempty(children) && return nothing

    sort!(children; lt = simplex_less)
    return SiblingStream(children, 1)
end

# Add stream to the enumerator 
# Add head to heap if nonempty
@inline function add_stream!(
    fe::FiltrationEnumerator,
    st::SiblingStream,
)
    push!(fe.streams, st)
    sid = length(fe.streams)

    if !is_empty(st)
        push_head!(fe.heap, HeadEntry(head(st), sid))
    end

    return nothing 
end

# Advances the stream at position sid up one and adds new head to heap
@inline function advance_head!(
    fe::FiltrationEnumerator,
    sid::Int,
)
    st = fe.streams[sid]
    advance!(st)

    if !is_empty(st)
        push_head!(fe.heap, HeadEntry(head(st), sid))
    end

    return nothing
end

# Returns the next simplex and updates the enumerator
function pop_simplex!(fe::FiltrationEnumerator)
    isempty(fe.heap) && return nothing

    head = pop_head!(fe.heap)
    s = head.simplex

    advance_head!(fe, head.id)

    child_stream = child_stream!(fe, s)
    if child_stream !== nothing
        add_stream!(fe, child_stream)
    end

    return s
end

# Convenience wrapper
function next_simplex!(fe::FiltrationEnumerator)
    return pop_simplex!(fe)
end

# ============================================================
# Dictionaries
# ============================================================

const MaxSet = Set{Int}

# toggle an memberership in a set
@inline function toggle!(s::Set{Int}, x::Int)
    if x in s
        delete!(s, x)
    else
        push!(s, x)
    end
    return s
end

# symmetric difference
@inline function xor_add!(dst::Set{Int}, src::Set{Int})
    @inbounds for x in src
        if x in dst
            delete!(dst, x)
        else
            push!(dst, x)
        end
    end
    return dst
end

# R[x] xor dst if x a key else nothing
@inline function xor_from!(
    dst::Set{Int},
    R::Dict{Int,MaxSet},
    x::Int,
)
    sx = get(R, x, nothing)
    sx === nothing && return false
    xor_add!(dst, sx)
    return true
end

# insert new generator
@inline function insert_live!(
    R::Dict{Int,MaxSet},
    Ri::Dict{Int,MaxSet},
    x::Int,
)
    sx = Set{Int}()
    push!(sx, x)
    R[x] = sx

    rx = Set{Int}()
    push!(rx, x)
    Ri[x] = rx

    return nothing
end

# checks if active 
@inline function is_live(
    R::Dict{Int,MaxSet},
    x::Int,
)
    return haskey(R, x)
end

# removes y in the case that del_bar x = {y}
function remove_maximal!(
    R::Dict{Int,MaxSet},
    Ri::Dict{Int,MaxSet},
    y::Int,
    removed::Vector{Int},
)
    empty!(removed)

    users = get(Ri, y, nothing)

    if users === nothing
        set_y = get(R, y, nothing)
        if set_y !== nothing && length(set_y) == 1 && (y in set_y)
            delete!(R, y)
            push!(removed, y)
        end
        return removed
    end

    touched = collect(users)
    delete!(Ri, y)

    @inbounds for z in touched
        set_z = get(R, z, nothing)
        set_z === nothing && continue

        delete!(set_z, y)

        if isempty(set_z)
            delete!(R, z)
            push!(removed, z)
        end
    end

    return removed
end

# handles del_bar = {y_1, ..., y_k}
function replace_pivot!(
    R::Dict{Int,MaxSet},
    Ri::Dict{Int,MaxSet},
    bar::Set{Int},
    j::Int,
    removed::Vector{Int},
)
    empty!(removed)

    users = get(Ri, j, nothing)
    users === nothing && return removed

    touched = collect(users)

    others = Int[]
    sizehint!(others, max(length(bar) - 1, 0))
    @inbounds for y in bar
        y == j && continue
        push!(others, y)
    end

    @inbounds for z in touched
        set_z = get(R, z, nothing)
        set_z === nothing && continue

        # xor in the whole bar
        xor_add!(set_z, bar)

        # toggle z in all non-pivot reverse sets
        for y in others
            set_y = get!(Ri, y) do
                Set{Int}()
            end
            toggle!(set_y, z)
            isempty(set_y) && delete!(Ri, y)
        end

        if isempty(set_z)
            delete!(R, z)
            push!(removed, z)
        end
    end

    # remove the pivot reverse-entry after the full toggling process
    delete!(Ri, j)

    return removed
end

# find maximal entry in del_bar
function pick_pivot!(
    bar::Set{Int},
    dim::Int,
    dist::DistanceMatrix,
    nverts::Int,
    binom::BinomialTable,
    buf::Vector{Int},
)::Tuple{Int,ValueT}
    first_seen = true
    best_idx = 0
    best_diam = zero(ValueT)

    @inbounds for j in bar
        dj = diameter!(buf, j, dim, dist, nverts, binom)
        if first_seen || (dj > best_diam) || ((dj == best_diam) && (j > best_idx))
            best_idx = j
            best_diam = dj
            first_seen = false
        end
    end

    return best_idx, best_diam
end

# ============================================================
# Active enumeration helpers
# ============================================================

# compute the diameter of a simplex from idx
@inline function diameter!(
    buf::Vector{Int},
    idx::IndexT,
    dim::Int,
    dist::DistanceMatrix,
    n::Int,
    binom::BinomialTable,
)::ValueT
    decode_simplex!(buf, idx, dim, n, binom)

    m = dim + 1
    diam = zero(ValueT)
    mat = dist.mat

    @inbounds for a in 1:m
        va = buf[a] + 1
        for b in (a + 1):m
            d = mat[va, buf[b] + 1]
            if d > diam
                diam = d
            end
        end
    end

    return diam
end

# given verts and one to omit compute the face idx
@inline function face_index_without!(
    verts::Vector{Int},
    m::Int,
    omit::Int,
    binom::BinomialTable,
)::IndexT
    idx = zero(IndexT)

    @inbounds begin
        for k in 1:(omit - 1)
            idx += choose(binom, verts[k], k)
        end
        for k in (omit + 1):m
            idx += choose(binom, verts[k], k - 1)
        end
    end

    return idx
end

# checks if w is in verts
@inline function contains_vertex(
    verts::Vector{Int},
    m::Int,
    w::Int,
)::Bool
    @inbounds for i in 1:m
        if verts[i] == w
            return true
        end
    end
    return false
end

# remove entry from vect of tuples where first entry is idx
function remove_face_entry!(
    lst::Vector{Tuple{Int,Int}},
    idx::Int,
)
    @inbounds for k in eachindex(lst)
        if lst[k][1] == idx
            lst[k] = lst[end]
            pop!(lst)
            return
        end
    end
    return nothing
end

# merge w into sorted xverts
@inline function merge_vertex!(
    out::Vector{Int},
    xverts::Vector{Int},
    m::Int,
    w::Int,
)::Int
    i = 1
    j = 1

    @inbounds while i <= m && xverts[i] < w
        out[j] = xverts[i]
        i += 1
        j += 1
    end

    out[j] = w
    wpos = j
    j += 1

    @inbounds while i <= m
        out[j] = xverts[i]
        i += 1
        j += 1
    end

    return wpos
end

# compute diam of face by omitting vert in index omit
@inline function face_diameter!(
    verts::Vector{Int},
    m::Int,
    omit::Int,
    dist::DistanceMatrix,
)::ValueT
    diam = zero(ValueT)
    mat = dist.mat

    @inbounds for a in 1:m
        a == omit && continue
        va = verts[a] + 1
        for b in (a + 1):m
            b == omit && continue
            d = mat[va, verts[b] + 1]
            if d > diam
                diam = d
            end
        end
    end

    return diam
end

# check whether xU{w} ready to appear in filtration, i.e. that x is the youngest face
function coface_ready!(
    candverts::Vector{Int},
    xverts::Vector{Int},
    w::Int,
    xidx::Int,
    xdiam::ValueT,
    dist::DistanceMatrix,
    binom::BinomialTable,
)::Bool
    m = length(xverts)
    total = m + 1

    wpos = merge_vertex!(candverts, xverts, m, w)

    @inbounds for omit in 1:total
        omit == wpos && continue

        faceidx = face_index_without!(candverts, total, omit, binom)

        if faceidx < xidx
            continue
        end

        facedia = face_diameter!(candverts, total, omit, dist)
        facedia < xdiam || return false
    end

    return true
end

# register an n-simplex as active. Rather than just using keys store extra info to compute quicker
function register_active!(
    active_verts::Dict{Int,Vector{Int}},
    face_to_active::Dict{IndexT,Vector{Tuple{Int,Int}}},
    idx::Int,
    verts::Vector{Int},
    n::Int,
    binom::BinomialTable,
)
    stored = copy(verts)
    active_verts[idx] = stored # store vertices

    m = n + 1
    @inbounds for omit in 1:m # loop over faces
        faceidx = face_index_without!(stored, m, omit, binom) # face idx
        lst = get!(face_to_active, faceidx) do  
            Tuple{Int,Int}[]
        end
        push!(lst, (idx, stored[omit])) # active simplex contains that face
    end

    return nothing
end

# unregister active face
function unregister_active!(
    active_verts::Dict{Int,Vector{Int}},
    face_to_active::Dict{IndexT,Vector{Tuple{Int,Int}}},
    doomed::Vector{Int},
    n::Int,
    binom::BinomialTable,
)
    m = n + 1

    for idx in doomed
        verts = get(active_verts, idx, nothing)
        verts === nothing && continue

        @inbounds for omit in 1:m
            faceidx = face_index_without!(verts, m, omit, binom)
            lst = face_to_active[faceidx]
            remove_face_entry!(lst, idx)
            isempty(lst) && delete!(face_to_active, faceidx)
        end

        delete!(active_verts, idx)
    end

    return nothing
end

# stores info about top level search 
mutable struct SeedState
    seen::Vector{Bool} 
    verts::Vector{Int} # vertices that have been seen
    seeds::Vector{Vector{Int}}
end

# initializes state
function SeedState(nverts::Int)
    seen = Vector{Bool}(undef, nverts)
    fill!(seen, false)
    seeds = Vector{Vector{Int}}(undef, nverts)
    for i in 1:nverts
        seeds[i] = Int[]
    end
    return SeedState(seen, Int[], seeds)
end

# resets state
function reset_seeds!(state::SeedState)
    for w in state.verts
        state.seen[w + 1] = false
        empty!(state.seeds[w + 1])
    end
    empty!(state.verts)
    return nothing
end

# find candidate vertices that may form n+1 simplex
# for each face, find actives with that face 
# check that radi < r 
# if so, add extra active spx to state.seeds[wi] to record that active vert seeded to this w
function collect_candidates!(
    state::SeedState,
    xverts::Vector{Int},
    r::ValueT,
    dist::DistanceMatrix,
    n::Int,
    binom::BinomialTable,
    face_to_active::Dict{IndexT,Vector{Tuple{Int,Int}}},
)
    reset_seeds!(state)

    m = n + 1
    mat = dist.mat

    @inbounds for omit in 1:m
        faceidx = face_index_without!(xverts, m, omit, binom)
        lst = get(face_to_active, faceidx, nothing)
        lst === nothing && continue

        xextra = xverts[omit] + 1

        for pair in lst
            yidx = pair[1]
            w = pair[2]

            if mat[xextra, w + 1] <= r
                wi = w + 1
                if !state.seen[wi]
                    state.seen[wi] = true
                    push!(state.verts, w)
                end
                push!(state.seeds[wi], yidx)
            end
        end
    end

    return nothing
end

# ============================================================
# Lazy r-close neighbor cache
# ============================================================

# section for if x not killed by looping though actives 

struct NeighborCache
    order::Vector{Vector{Int}} # for each vertex v, other verts sorted by dist to v
    limit::Vector{Int} # how many vertices are w/in current radius of v
end

# initialize cache
function build_neighbor_cache(
    dist::DistanceMatrix,
    nverts::Int,
)::NeighborCache
    order = Vector{Vector{Int}}(undef, nverts)
    nn = nverts - 1
    mat = dist.mat

    ids = Vector{Int}(undef, max(nn, 0))
    ds = Vector{ValueT}(undef, max(nn, 0))
    perm = Vector{Int}(undef, max(nn, 0))

    @inbounds for v in 0:(nverts - 1)
        vp1 = v + 1
        row = Vector{Int}(undef, nn)

        t = 1
        for w in 0:(nverts - 1)
            w == v && continue
            ids[t] = w
            ds[t] = mat[vp1, w + 1]
            perm[t] = t
            t += 1
        end

        sort!(perm, alg = Base.Sort.DEFAULT_STABLE, by = i -> ds[i])

        for i in 1:nn
            row[i] = ids[perm[i]]
        end

        order[vp1] = row
    end

    return NeighborCache(order, zeros(Int, nverts))
end

# advance the limit as r increases
@inline function advance_limit!(
    cache::NeighborCache,
    v::Int,
    r::ValueT,
    dist::DistanceMatrix,
)::Int
    row = cache.order[v + 1]
    limit = cache.limit[v + 1]
    L = length(row)
    mat = dist.mat
    vp1 = v + 1

    @inbounds while limit < L && mat[vp1, row[limit + 1] + 1] <= r
        limit += 1
    end

    cache.limit[v + 1] = limit
    return limit
end

# find vertex with fewest neighbors to minimize checking
@inline function choose_pivot_vertex!(
    xverts::Vector{Int},
    m::Int,
    r::ValueT,
    dist::DistanceMatrix,
    cache::NeighborCache,
)::Tuple{Int,Int}
    best_pos = 1
    best_lim = typemax(Int)

    @inbounds for i in 1:m
        limit = advance_limit!(cache, xverts[i], r, dist)
        if limit < best_lim
            best_lim = limit
            best_pos = i
        end
    end

    return best_pos, best_lim
end

# checks if x has a leftover n+1 coface that kills it 
function has_leftover_coface!(
    seen_w::Vector{Bool},
    xverts::Vector{Int},
    xidx::Int,
    r::ValueT,
    dist::DistanceMatrix,
    binom::BinomialTable,
    neighbors::NeighborCache,
    candverts::Vector{Int},
)
    m = length(xverts)

    pivot_pos, limit = choose_pivot_vertex!(xverts, m, r, dist, neighbors)
    pivot = xverts[pivot_pos]
    row = neighbors.order[pivot + 1]
    mat = dist.mat

    @inbounds for t in 1:limit
        w = row[t]

        contains_vertex(xverts, m, w) && continue
        seen_w[w + 1] && continue

        ok = true
        wp1 = w + 1
        for i in 1:m
            i == pivot_pos && continue
            if mat[wp1, xverts[i] + 1] > r
                ok = false
                break
            end
        end

        ok || continue

        if coface_ready!(candverts, xverts, w, xidx, r, dist, binom)
            return true
        end
    end

    return false
end

# ============================================================
# Active enumeration
# ============================================================

# process a top level coface
function process_seeded!(
    R::Dict{Int,MaxSet},
    Ri::Dict{Int,MaxSet},
    top_pairs::Vector{Tuple{ValueT,ValueT}},
    xidx::Int,
    other_seed_faces::Vector{Int},
    top_diam::ValueT,
    fe::FiltrationEnumerator,
    top_bar::Set{Int},
    diam_buf::Vector{Int},
    removed::Vector{Int},
    active_verts::Dict{Int,Vector{Int}},
    face_to_active::Dict{IndexT,Vector{Tuple{Int,Int}}},
    n::Int,
)
    empty!(top_bar)

    xor_from!(top_bar, R, xidx)

    @inbounds for yidx in other_seed_faces
        xor_from!(top_bar, R, yidx)
    end

    isempty(top_bar) && return nothing

    if length(top_bar) == 1
        j = first(top_bar)
        remove_maximal!(R, Ri, j, removed)
        unregister_active!(
            active_verts,
            face_to_active,
            removed,
            n,
            fe.binom,
        )

        jdiam = diameter!(
            diam_buf,
            j,
            n,
            fe.dist,
            fe.n,
            fe.binom,
        )

        if top_diam != jdiam
            push!(top_pairs, (jdiam, top_diam))
        end

        return nothing
    end

    j, jdiam = pick_pivot!(
        top_bar,
        n,
        fe.dist,
        fe.n,
        fe.binom,
        diam_buf,
    )

    replace_pivot!(R, Ri, top_bar, j, removed)

    if !isempty(removed)
        unregister_active!(
            active_verts,
            face_to_active,
            removed,
            n,
            fe.binom,
        )
    end

    if top_diam != jdiam
        push!(top_pairs, (jdiam, top_diam))
    end

    return nothing
end

# process leftover if x still alive
function process_leftover!(
    R::Dict{Int,MaxSet},
    Ri::Dict{Int,MaxSet},
    top_pairs::Vector{Tuple{ValueT,ValueT}},
    xidx::Int,
    top_diam::ValueT,
    fe::FiltrationEnumerator,
    diam_buf::Vector{Int},
    removed::Vector{Int},
    active_verts::Dict{Int,Vector{Int}},
    face_to_active::Dict{IndexT,Vector{Tuple{Int,Int}}},
    n::Int,
)
    haskey(R, xidx) || return nothing

    if haskey(Ri, xidx)
        j = xidx

        remove_maximal!(R, Ri, j, removed)
        unregister_active!(
            active_verts,
            face_to_active,
            removed,
            n,
            fe.binom,
        )

        jdiam = diameter!(
            diam_buf,
            j,
            n,
            fe.dist,
            fe.n,
            fe.binom,
        )

        if top_diam != jdiam
            push!(top_pairs, (jdiam, top_diam))
        end

        return nothing
    end

    top_bar = Set(R[xidx])

    isempty(top_bar) && return nothing

    if length(top_bar) == 1
        j = first(top_bar)

        remove_maximal!(R, Ri, j, removed)
        unregister_active!(
            active_verts,
            face_to_active,
            removed,
            n,
            fe.binom,
        )

        jdiam = diameter!(
            diam_buf,
            j,
            n,
            fe.dist,
            fe.n,
            fe.binom,
        )

        if top_diam != jdiam
            push!(top_pairs, (jdiam, top_diam))
        end

        return nothing
    end

    j, jdiam = pick_pivot!(
        top_bar,
        n,
        fe.dist,
        fe.n,
        fe.binom,
        diam_buf,
    )

    replace_pivot!(R, Ri, top_bar, j, removed)

    if !isempty(removed)
        unregister_active!(
            active_verts,
            face_to_active,
            removed,
            n,
            fe.binom,
        )
    end

    if top_diam != jdiam
        push!(top_pairs, (jdiam, top_diam))
    end

    return nothing
end

# ============================================================
# Infinite intervals
# ============================================================

# find infinite intervals
function append_infinite!(
    out::Vector{Tuple{ValueT,ValueT}},
    Ri::Dict{Int,MaxSet},
    dim::Int,
    dist::DistanceMatrix,
    nverts::Int,
    binom::BinomialTable,
    diam_buf::Vector{Int},
)
    live_maxima = collect(keys(Ri))
    sort!(live_maxima)

    infv = ValueT(Inf)

    @inbounds for x in live_maxima
        birth = diameter!(
            diam_buf,
            x,
            dim,
            dist,
            nverts,
            binom,
        )
        push!(out, (birth, infv))
    end

    return out
end

# ============================================================
# Main algorithm
# ============================================================

function redzed(
    dist::DistanceMatrix,
    n::Int;
    threshold::Real = Inf
)
    n >= 1 || throw(ArgumentError("requires n >= 1"))

    # value where it becomes a cone
    cone_val=minimum(maximum(dist.mat, dims=1))
    if cone_val<threshold
        threshold=cone_val
    end

    fe = FiltrationEnumerator(dist; max_dim = n, threshold = threshold)
    be = BoundaryEnumerator()

    R = [Dict{Int,MaxSet}() for _ in 0:n]
    Ri = [Dict{Int,MaxSet}() for _ in 0:n]
    pairs = [Vector{Tuple{ValueT,ValueT}}() for _ in 0:n]

    bar = Set{Int}()
    top_bar = Set{Int}()

    diam_bufs = [Vector{Int}(undef, d + 1) for d in 0:n]
    top_buf = Vector{Int}(undef, n + 1)
    coface_buf = Vector{Int}(undef, n + 2)

    active_verts = Dict{Int,Vector{Int}}()
    face_to_active = Dict{IndexT,Vector{Tuple{Int,Int}}}()

    seeds = SeedState(fe.n)
    neighbors = build_neighbor_cache(fe.dist, fe.n)
    removed = Int[]

    # main loop
    while true
        s = next_simplex!(fe)
        s === nothing && break # check that ther is a next simplex

        if s.dim == 0 # special case for vertex
            insert_live!(R[1], Ri[1], Int(s.idx))
            continue
        end

        empty!(bar) # clear del_bar
        set_simplex!(be, s.idx, s.dim, fe.n)

        while has_next_face(be) # enumerate faces
            face = Int(next_face_index!(be, fe.binom))
            xor_from!(bar, R[s.dim], face)
        end

        if isempty(bar) # case 1
            xidx = Int(s.idx)
            insert_live!(R[s.dim + 1], Ri[s.dim + 1], xidx)

            if s.dim == n # if n dimensional search for top level simplices
                decode_simplex!(top_buf, s.idx, n, fe.n, fe.binom)

                collect_candidates!(
                    seeds,
                    top_buf,
                    s.diam,
                    fe.dist,
                    n,
                    fe.binom,
                    face_to_active,
                )

                for w in seeds.verts
                    if coface_ready!(
                        coface_buf,
                        top_buf,
                        w,
                        xidx,
                        s.diam,
                        fe.dist,
                        fe.binom,
                    )
                        process_seeded!(
                            R[n + 1],
                            Ri[n + 1],
                            pairs[n + 1],
                            xidx,
                            seeds.seeds[w + 1],
                            s.diam,
                            fe,
                            top_bar,
                            diam_bufs[n + 1],
                            removed,
                            active_verts,
                            face_to_active,
                            n,
                        )
                    end
                end

                if is_live(R[n + 1], xidx)
                    if has_leftover_coface!(
                        seeds.seen,
                        top_buf,
                        xidx,
                        s.diam,
                        fe.dist,
                        fe.binom,
                        neighbors,
                        coface_buf,
                    )
                        process_leftover!(
                            R[n + 1],
                            Ri[n + 1],
                            pairs[n + 1],
                            xidx,
                            s.diam,
                            fe,
                            diam_bufs[n + 1],
                            removed,
                            active_verts,
                            face_to_active,
                            n,
                        )
                    end
                end

                if is_live(R[n + 1], xidx)
                    register_active!(
                        active_verts,
                        face_to_active,
                        xidx,
                        top_buf,
                        n,
                        fe.binom,
                    )
                end
            end

        elseif length(bar) == 1 # case 2a (faster to break these cases up)
            j = first(bar)

            remove_maximal!(R[s.dim], Ri[s.dim], j, removed)

            jdiam = diameter!(
                diam_bufs[s.dim],
                j,
                s.dim - 1,
                fe.dist,
                fe.n,
                fe.binom,
            )

            if s.diam != jdiam
                push!(pairs[s.dim], (jdiam, s.diam))
            end

        else # case 2b
            j, jdiam = pick_pivot!(
                bar,
                s.dim - 1,
                fe.dist,
                fe.n,
                fe.binom,
                diam_bufs[s.dim],
            )

            replace_pivot!(R[s.dim], Ri[s.dim], bar, j, removed)

            if s.diam != jdiam
                push!(pairs[s.dim], (jdiam, s.diam))
            end
        end
    end

    # infinite intervals
    @inbounds for p in 0:n
        append_infinite!(
            pairs[p + 1],
            Ri[p + 1],
            p,
            fe.dist,
            fe.n,
            fe.binom,
            diam_bufs[p + 1],
        )
    end

    return Dict(d => pairs[d + 1] for d in 0:n)
end

# for abstract matrix
function redzed(
    dist::AbstractMatrix,
    n::Int;
    threshold::Real = Inf
)
    return redzed(DistanceMatrix(dist), n; threshold) 
end
