import re
import sys
import math

def parse_ripser_style_file(filepath):
    intervals = {}
    current_dim = -1
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "persistence intervals in dim" in line:
                match = re.search(r'dimension\s+(\d+)', line)
                if match:
                    current_dim = int(match.group(1))
                    intervals[current_dim] = []
            elif line.startswith("[") or line.startswith(" ["):
                # Clean characters and extract pairs
                cleaned = line.replace("[", "").replace(")", "").replace(" ", "")
                parts = cleaned.split(",")
                if len(parts) == 2:
                    b = float(parts[0])
                    d = float(parts[1]) if parts[1].lower() != 'inf' else float('inf')
                    intervals[current_dim].append((b, d))
                    
    # Sort for exact index mapping alignment
    for dim in intervals:
        intervals[dim].sort(key=lambda x: (x[0], x[1]))
    return intervals

def compare_results(file1, file2, tol=1e-5):
    print("Reading", file1)
    data1 = parse_ripser_style_file(file1)
    print("Reading", file2)
    data2 = parse_ripser_style_file(file2)
    
    all_dims = set(data1.keys()).union(set(data2.keys()))
    
    for dim in sorted(all_dims):
        pairs1 = data1.get(dim, [])
        pairs2 = data2.get(dim, [])
        
        if len(pairs1) != len(pairs2):
            print(f"❌ Dimension {dim} mismatch! Count: File1={len(pairs1)}, File2={len(pairs2)}")
            return False
            
        for idx, ((b1, d1), (b2, d2)) in enumerate(zip(pairs1, pairs2)):
            b_match = math.isclose(b1, b2, abs_tol=tol) if not (math.isinf(b1) and math.isinf(b2)) else True
            d_match = math.isclose(d1, d2, abs_tol=tol) if not (math.isinf(d1) and math.isinf(d2)) else True
            
            if not (b_match and d_match):
                print(f"❌ Value Mismatch in Dim {dim} at index {idx}:")
                print(f"   File1: [{b1}, {d1})")
                print(f"   File2: [{b2}, {d2})")
                return False
                
    print("✅ Perfect Match! Persistence pairs match identically within tolerance limits.")
    return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python validate.py <file1.txt> <file2.txt>")
        sys.exit(1)
    print("Comparing", sys.argv[1], sys.argv[2])
    compare_results(sys.argv[1], sys.argv[2])
