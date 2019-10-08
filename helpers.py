import math

def divide_into_equal_passes(begin, end, max_increment, include_first = True):
    distance = float(end - begin)
    num_steps = int(math.ceil(abs(distance) / max_increment))
    actual_increment = distance / num_steps
    passes = [begin + i * actual_increment for i in range(num_steps + 1)]
    passes[-1] = end # Remove rounding errors
    first_idx = 0 if include_first else 1
    return passes[first_idx:]

def center_equally_spaced_points_in_range(begin, end, increment, min_edge_size = 0.0):
    min_val = min(begin, end)
    max_val = max(begin, end)
    distance = float(max_val - min_val - 2 * min_edge_size)
    num_steps = abs(int(math.floor(distance / increment)))
    covered_distance = increment * num_steps
    remainder = distance - covered_distance
    actual_edge_size = min_edge_size + remainder / 2.0
    edges = [min_val + actual_edge_size + i * increment for i in range(num_steps+1)]
    direction = -1 if begin > end else 1
    return edges[::direction]

def divide_with_clearout_depths(max_depth, depth_increment, clearout_depth):
    # First, divide depth evenly
    basic_depths = divide_into_equal_passes(0.0, max_depth, depth_increment, include_first = False)
    # Second, add clearout passes
    depths = []
    cleared_depth = clearout_depth
    for depth_idx, depth in enumerate(basic_depths):
        depths.append(depth)
        if depth >= cleared_depth and depth_idx != len(basic_depths)-1:
            cleared_depth += clearout_depth
            depths.append(0.5*depth)
            depths.append(depth)
    print(depths)
    return depths