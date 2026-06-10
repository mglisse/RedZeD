# Import the provided redzed implementation
include("redzed.jl")

using LinearAlgebra

# Parse a point cloud file (space or comma-separated coordinates, one point per line)
function load_point_cloud(filename)
    points = Vector{Vector{Float64}}()
    for line in eachline(filename)
        line = strip(line)
        if isempty(line) || startswith(line, "#") continue end
        row = parse.(Float64, split(replace(line, "," => " ")))
        push!(points, row)
    end
    
    n = length(points)
    dist_mat = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        dist_mat[i, j] = sqrt(sum((points[i] .- points[j]).^2))
    end
    return dist_mat
end

# Parse a STRICTLY lower-triangular distance matrix file
function load_lower_distance_matrix(filename)
    values = Float64[]
    for line in eachline(filename)
        line = strip(line)
        if isempty(line) || startswith(line, "#") continue end
        row_vals = parse.(Float64, split(replace(line, "," => " ")))
        append!(values, row_vals)
    end
    
    M = length(values)
    discriminant = 1 + 8 * M
    sqrt_disc = round(Int, sqrt(discriminant))
    if sqrt_disc * sqrt_disc != discriminant
        error("The number of entries ($M) does not form a valid strictly lower triangular matrix.")
    end
    n = (1 + sqrt_disc) ÷ 2
    
    dist_mat = zeros(Float64, n, n)
    idx = 1
    for i in 1:n
        for j in 1:(i-1)
            dist_mat[i, j] = values[idx]
            dist_mat[j, i] = values[idx]
            idx += 1
        end
    end
    return dist_mat
end

function main()
    if length(ARGS) < 4
        println("Usage: julia run_redzed.jl <format: pc|txt> <max_dim> <threshold> <input_file>")
        println("  pc:  Point Cloud format")
        println("  txt: Strictly Lower triangular distance matrix format")
        exit(1)
    end

    format = ARGS[1]
    max_dim = parse(Int, ARGS[2])
    threshold = parse(Float64, ARGS[3])
    input_file = ARGS[4]

    # 1. Load Data
    if format == "pc"
        dist_matrix = load_point_cloud(input_file)
    elseif format == "txt"
        dist_matrix = load_lower_distance_matrix(input_file)
    else
        error("Unknown format profile. Use 'pc' or 'txt'.")
    end

    println("Dataset loaded. Size: $(size(dist_matrix, 1))x$(size(dist_matrix, 2))")
    println("Running redzed warm-up...")
    # Warm up Julia's JIT compiler to ensure timing accuracy
    _ = redzed(dist_matrix[1:min(5, end), 1:min(5, end)], min(max_dim, 1), threshold=threshold)

    println("Benchmarking core algorithm...")
    # 2. Timing the execution
    start_time = time_ns()
    results = redzed(dist_matrix, max_dim, threshold=threshold)
    end_time = time_ns()
    
    elapsed_ms = (end_time - start_time) / 1_000_000.0
    println("Execution Time (Julia): $(elapsed_ms) ms")

    # 3. Output persistence pairs for validation matching Ripser style
    open("output_julia.txt", "w") do f
        for d in 0:max_dim
            write(f, "persistence intervals in dimension $d:\n")
            if haskey(results, d)
                sorted_pairs = sort(results[d], by = x -> (x[1], x[2]))
                for (birth, death) in sorted_pairs
                    write(f, " [$birth,$death)\n")
                end
            end
        end
    end
    println("Results written to output_julia.txt")
end

main()
