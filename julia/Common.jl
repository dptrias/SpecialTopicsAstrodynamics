# Common.jl
using StaticArrays, LinearAlgebra

"""
    Structure to represent a solar sail with lightness number β and normal vector n.
"""
const SolarSail = Tuple{Float64,SVector{2,Float64}}

""" 
    load_solar_system_parameters() -> Dict{String, Any}
"""
function load_solar_system_parameters()::Dict{String,Any}
    return Dict(
        "μ_SE" => 3.0542e-6,
        "μ_SM" => 3.227154996101724e-7,
        "μ_S" => 132712e15,          # m^3/s^2
        "S_sun" => 1361.0,             # W/m^2 at 1 AU
        "c" => 299792458.0,        # m/s
        "LU_SE" => 149597870.7e3,     # m
        "LU_SM" => 208321282.0e3,     # m
    )
end


"""
    classical_lagrange_points(μ, point) -> Float64

Return the x-coordinate (in normalised units) of the L1 (point=1) or L2 (point=2)
Lagrange point for a CR3BP system with mass parameter `μ`.

Uses Newton-Raphson iteration on the collinear equilibrium condition.
"""
function classical_lagrange_points(μ::Float64, point::Int)::Float64

    f(x) = x - (1 - μ) / (x + μ)^2 - μ / (x - (1 - μ))^2

    df(x) = 1 + 2 * (1 - μ) / (x + μ)^3 + 2 * μ / (x - (1 - μ))^3

    if point == 1
        x = 0.5
    elseif point == 2
        x = 1.5
    else
        error("Invalid point: choose 1 for L1 or 2 for L2.")
    end

    tol = 1e-6
    max_iter = 1000

    for _ in 1:max_iter
        fx = f(x)
        dfx = df(x)

        if dfx == 0.0
            error("Derivative is zero. No solution found.")
        end

        x_new = x - fx / dfx

        if abs(x_new - x) < tol
            println("L$point point found at x = $(round(x_new; digits=6))")
            return x_new
        end

        x = x_new
    end

    error("Newton-Raphson did not converge within $max_iter iterations.")
end

# Helper functions
"""
    eq_condition_solar_sail(r_AEP, μ) -> Float64, SVector{3,Float64}

Compute the lightness number β required to place an artificial equilibrium
point at `r_AEP` (normalised units) in the CR3BP with mass parameter `μ`.
"""
function eq_condition_solar_sail(r_AEP::SVector{3,Float64}, μ::Float64)::SolarSail
    r₁ = norm(r_AEP + SVector(μ, 0.0, 0.0))
    r₂ = norm(r_AEP - SVector(1 - μ, 0.0, 0.0))
    U_x = (1 - μ) * (r_AEP[1] + μ) / r₁^3 + μ * (r_AEP[1] - (1 - μ)) / r₂^3 - r_AEP[1]
    U_y = (1 - μ) * r_AEP[2] / r₁^3 + μ * r_AEP[2] / r₂^3 - r_AEP[2]
    U_z = (1 - μ) * r_AEP[3] / r₁^3 + μ * r_AEP[3] / r₂^3
    ∇U = SVector(U_x, U_y, U_z)
    n = ∇U / norm(∇U)
    r̂₁ = (r_AEP + SVector(μ, 0.0, 0.0)) / r₁
    return SolarSail(r₁^2 * norm(∇U) / ((1 - μ) * dot(r̂₁, n)^2), n)
end

"""
    solar_radiation_pressure(S, r) -> Float64

Compute the solar radiation pressure (Pa) at distance `r` (m) from the Sun,
given solar irradiance `S` (W/m^2) at 1 AU.
"""
function solar_radiation_pressure(S::Float64, r::Float64)::Float64
    c_light = 299792458.0
    AU = 149597870.7e3   # m
    return S / c_light * (AU / r)^2
end

"""
    linear_stability_analysis(μ, r₀, solar_sail=nothing) -> Vector{ComplexF64}

Perform linear stability analysis of an equilibrium point at `r₀` in the CR3BP with mass parameter `μ`. 
If `solar_sail` is provided, the analysis includes the effect of the solar sail with given lightness number and normal vector.
"""
function linear_stability_analysis(
    μ::Float64,
    r₀::Vector{Float64},
    solar_sail::Union{SolarSail,Nothing}=nothing
)::Vector{ComplexF64}
    N = length(r₀)

    r₁ = r₀
    r₁[1] += μ
    r₂ = r₀
    r₂[1] -= (1 - μ)

    r₁_sq = norm(r₁)^2
    r₂_sq = norm(r₂)^2
    Id = Matrix{Float64}(I, N, N)

    Uᵢⱼ = (1 - μ) * (Id - 3 * r₁ * r₁' / r₁_sq) / r₁_sq^(3 / 2) +
          μ * (Id - 3 * r₂ * r₂' / r₂_sq) / r₂_sq^(3 / 2) - Id
    if N == 3
        Uᵢⱼ[3, 3] += 1
    end

    A = [zeros(N, N) Id;
        -Uᵢⱼ zeros(N, N)]
    A[N+1, N+2] = 2
    A[N+2, N+1] = -2

    if solar_sail !== nothing
        β, n = solar_sail
        n_dot_r = dot(n, r₁)
        ζ = zeros(N, 1)
        for i in 1:N
            ζ[i] = 4 * β * n_dot_r * (1 - μ) / r₁_sq^3 *
                   (n_dot_r * r₁[i] - n[i] * r₁_sq / 2)
        end
        a_sᵢⱼ = ζ * n'
        A[N+1:N+N, 1:N] += a_sᵢⱼ
    end

    eigenvalues = eigvals(A)
    return eigenvalues
end