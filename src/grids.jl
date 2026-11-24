# Grid generator.

"""
    grid_refined(L; h_min, h_max=h_min, growth=1.1) -> Vector{Float64}

Node coordinates of the interval [0, L], refined toward both ends. The cell size is `h_min`
at x = 0 and x = L and grows geometrically by the factor `growth` toward the centre, up to
`h_max`. The grid is symmetric.
"""
function grid_refined(L::Real; h_min::Real, h_max::Real=h_min, growth::Real=1.1)
    L > 0 || throw(ArgumentError("The length must be positive, got $L"))
    0 < h_min <= h_max || throw(ArgumentError("The cell sizes must satisfy 0 < h_min ≤ h_max"))
    growth >= 1 || throw(ArgumentError("The growth factor must be ≥ 1, got $growth"))
    h_min < L / 2 || throw(ArgumentError("h_min must be smaller than L/2"))
    half = L / 2
    widths = Float64[]
    h, total = Float64(h_min), 0.0
    while total + h < half
        push!(widths, h)
        total += h
        h = min(h * growth, h_max)
    end
    # The remaining length is either merged with the last cell or added as a cell
    rest = half - total
    if rest < 0.5 * widths[end]
        widths[end] += rest
    else
        push!(widths, rest)
    end
    edges = cumsum(widths)
    edges[end] = half
    left = vcat(0.0, edges)
    right = L .- reverse(vcat(0.0, edges[1:(end - 1)]))
    return vcat(left, right)
end
