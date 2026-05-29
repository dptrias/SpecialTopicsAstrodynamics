# Common.jl
using StaticArrays, LinearAlgebra, Roots, DifferentialEquations

"""
    Structure to represent a solar sail with lightness number β and normal vector n.
"""
const SolarSail = Tuple{Float64,Vector{Float64}}

""" 
    load_solar_system_parameters() -> Dict{String, Any}

Return a dictionary of solar system constants and normalisation parameters for
the Sun-Earth (SE) and Sun-Mars (SM) CR3BP systems. Keys include mass parameters
μ_SE and μ_SM, gravitational constants, solar irradiance, speed of light,
characteristic length and time units, and the gravitational parameter of Mars.
"""
function load_solar_system_parameters()::Dict{String,Any}
    return Dict(
        "μ_SE" => 3.0542e-6,
        "μ_SM" => 3.227154996101724e-7,
        "μ_S" => 132712e15,                         # m^3/s^2
        "S_Sun" => 1361.0,                          # W/m^2 at 1 AU
        "c" => 299792458.0,                         # m/s
        "LU_SE" => 149597870.7e3,                   # m
        "LU_SM" => 208321282.0e3,                   # m
        "TU_SE" => 365.256 * 24 * 3600.0 / (2 * π), # s
        "TU_SM" => 8253622.0                        # s    
    )
end

"""
    classical_lagrange_point(μ, point) -> Vector{Float64}

Return the position (in normalised units) of the classical Lagrange point
`point` (1 through 5) for a CR3BP system with mass parameter `μ`.

Collinear points L1, L2, L3 are found via Newton-Raphson iteration on the
collinear equilibrium condition. L4 and L5 are returned analytically.
Returns a 2-element vector [x, y].
"""
function classical_lagrange_point(μ::Float64, point::Int)::Vector{Float64}
    tol = 1e-12
    max_iter = 1000

    if point == 4
        # L4
        return [0.5 - μ, sqrt(3) / 2]

    elseif point == 5
        # L5
        return [0.5 - μ, -sqrt(3) / 2]
    end

    function collinear_newton(f_and_df, x₀)
        x = x₀
        for _ in 1:max_iter
            fx, dfx = f_and_df(x)
            iszero(dfx) && error("Derivative is zero at x = $x.")
            x_new = x - fx / dfx
            abs(x_new - x) < tol && return [x_new, 0.0]
            x = x_new
        end
        error("Newton-Raphson did not converge within $max_iter iterations.")
    end

    if point == 1
        x₀ = 1.0 - (μ / 3)^(1 / 3)
        return collinear_newton(x -> (
                x - (1 - μ) / (x + μ)^2 + μ / (x - (1 - μ))^2,
                1 + 2(1 - μ) / (x + μ)^3 - 2μ / (x - (1 - μ))^3
            ), x₀)

    elseif point == 2
        x₀ = 1.0 + (μ / 3)^(1 / 3)
        return collinear_newton(x -> (
                x - (1 - μ) / (x + μ)^2 - μ / (x - (1 - μ))^2,
                1 + 2(1 - μ) / (x + μ)^3 + 2μ / (x - (1 - μ))^3
            ), x₀)

    elseif point == 3
        return collinear_newton(x -> (
                x + (1 - μ) / (x + μ)^2 + μ / (x - (1 - μ))^2,
                1 - 2(1 - μ) / (x + μ)^3 - 2μ / (x - (1 - μ))^3
            ), -1.0)

    else
        error("Invalid point: choose 1-5.")
    end

end

"""
    classical_lagrange_points(μ) -> Matrix{Float64}

Return a 2×5 matrix whose columns are the positions of all five classical
Lagrange points for a CR3BP system with mass parameter `μ`, in normalised
units. Column order is L1 through L5.
"""
function classical_lagrange_points(μ::Float64)::Matrix{Float64}
    return hcat(
        classical_lagrange_point(μ, 1),
        classical_lagrange_point(μ, 2),
        classical_lagrange_point(μ, 3),
        classical_lagrange_point(μ, 4),
        classical_lagrange_point(μ, 5)
    )
end

# Helper functions
"""
    eq_condition_solar_sail(r_AEP, μ) -> SolarSail

Compute the lightness number β required to place an artificial equilibrium
point at `r_AEP` (normalised units) in the CR3BP with mass parameter `μ`.

Also returns the required sail normal unit vector n. The result is a
`SolarSail` tuple (β, n).
"""
function eq_condition_solar_sail(r_AEP::Vector{Float64}, μ::Float64)::SolarSail
    r₁ = copy(r_AEP)
    r₁[1] += μ
    r₂ = copy(r_AEP)
    r₂[1] -= (1 - μ)

    ∇U = (1 - μ) * r₁ / norm(r₁)^3 + μ * r₂ / norm(r₂)^3 - r_AEP
    if length(r_AEP) == 3
        ∇U[3] += r_AEP[3]
    end
    n = ∇U / norm(∇U)
    r̂₁ = r₁ / norm(r₁)

    return norm(r₁)^2 * norm(∇U) / ((1 - μ) * dot(r̂₁, n)^2), n
end

"""
    collinear_AEP(β, μ, point) -> Tuple{Float64, Vector{Float64}}

Find the x-coordinate of the collinear artificial equilibrium point (AEP)
near classical Lagrange point `point` (1, 2, or 3) for a solar sail with
lightness number `β` in a CR3BP system with mass parameter `μ`.

Returns a tuple (x_AEP, n) where x_AEP is the normalised x-coordinate
and n is the sail normal direction (±x̂).
"""
function collinear_AEP(β::Float64, μ::Float64, point::Int)::Tuple{Float64,Vector{Float64}}
    function residual(x::Float64)::Float64
        r₁ = abs(x + μ)
        r₂ = abs(x - (1 - μ))
        U_x = (1 - μ) * (x + μ) / r₁^3 + μ * (x - (1 - μ)) / r₂^3 - x
        return β * (1 - μ) / r₁^2 - U_x
    end

    x_L = classical_lagrange_point(μ, point)[1]
    x_root = find_zero(residual, x_L - 0.0001)
    U_x = (1 - μ) * (x_root + μ) / abs(x_root + μ)^3 +
          μ * (x_root - (1 - μ)) / abs(x_root - (1 - μ))^3 - x_root
    return x_root, sign(U_x) * [1.0, 0.0, 0.0]
end

"""
    solar_radiation_pressure(S, r) -> Float64

Compute the solar radiation pressure (Pa) at distance `r` (m) from the Sun,
given solar irradiance `S` (W/m²) at 1 AU.
"""
function solar_radiation_pressure(S::Float64, r::Float64)::Float64
    c_light = 299792458.0
    AU = 149597870.7e3   # m
    return S / c_light * (AU / r)^2
end

"""
    linear_stability_analysis(r₀, μ[, solar_sail]) -> Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}

Perform linear stability analysis of an equilibrium point at `r₀` in the
CR3BP with mass parameter `μ`. Returns the eigenvalues and eigenvectors of
the linearised equations of motion.

If `solar_sail` is provided as a `SolarSail` tuple (β, n), the analysis
includes the first-order sail acceleration gradient. The state dimension is
inferred from the length of `r₀` (2D or 3D).
"""
function linear_stability_analysis(
    r₀::Vector{Float64},
    μ::Float64,
    solar_sail::Union{SolarSail,Nothing}=nothing
)::Tuple{Vector{ComplexF64},Matrix{ComplexF64}}
    N = length(r₀)

    r₁ = copy(r₀)
    r₁[1] += μ
    r₂ = copy(r₀)
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
        ζ = -4 * β * n_dot_r * (1 - μ) / r₁_sq^3 *
            (n_dot_r * r₁ - n * r₁_sq / 2)
        a_sᵢⱼ = ζ * n'
        A[N+1:N+N, 1:N] += a_sᵢⱼ
    end

    eigenvalues, eigenvectors = eigen(A)
    return eigenvalues, eigenvectors
end

"""
    cr3bp!(dr, r, p, t)

In-place ODE function for the planar CR3BP equations of motion.

State vector `r` is [x, y, vx, vy] in normalised units. Parameter vector
`p` must satisfy p[1] = μ (mass parameter). Equations include the Coriolis
and centrifugal terms of the rotating frame.
"""
function cr3bp!(dr::Vector{Float64}, r::Vector{Float64}, p::Vector{Float64}, t::Float64)::Nothing
    x, y, vx, vy = r
    mu = p[1]
    r₁ = sqrt((x + mu)^2 + y^2)
    r₂ = sqrt((x - (1 - mu))^2 + y^2)

    dr[1] = vx
    dr[2] = vy
    dr[3] = 2 * vy + x - (1 - mu) * (x + mu) / r₁^3 - mu * (x - (1 - mu)) / r₂^3
    dr[4] = -2 * vx + y - (1 - mu) * y / r₁^3 - mu * y / r₂^3
    return nothing
end

"""
    cr3bp_ss!(dr, r, p, t)

In-place ODE function for the planar CR3BP equations of motion with a solar
sail acceleraton.

State vector `r` is [x, y, vx, vy] in normalised units. Parameter vector
`p` must satisfy p[1] = μ (mass parameter), p[2] = β (sail lightness number),
and p[3] = α (sail cone angle in radians, measured from the Sun-line). The
sail is assumed to be a flat, perfectly reflecting surface.
"""
function cr3bp_ss!(dr::Vector{Float64}, r::Vector{Float64}, p::Vector{Float64}, t::Float64)::Nothing
    x, y, vx, vy = r
    mu = p[1]
    beta = p[2]
    alfa = p[3]
    r₁ = sqrt((x + mu)^2 + y^2)
    r₂ = sqrt((x - (1 - mu))^2 + y^2)

    aₛ = beta * (1 - mu) * cos(alfa)^2 / r₁^2
    dr[1] = vx
    dr[2] = vy
    dr[3] = 2 * vy + x - (1 - mu) * (x + mu) / r₁^3 - mu * (x - (1 - mu)) / r₂^3 + aₛ * ((x + mu) * cos(alfa) - y * sin(alfa)) / r₁
    dr[4] = -2 * vx + y - (1 - mu) * y / r₁^3 - mu * y / r₂^3 + aₛ * ((y) * cos(alfa) + (x + mu) * sin(alfa)) / r₁
    return nothing
end

"""
    compute_classical_manifold(r₀, μ, λ, v, t_f) -> Tuple{Array{Float64}, Array{Float64}}

Propagate the stable or unstable manifold branch associated with eigenvalue
`λ` and eigenvector `v` of the monodromy matrix at equilibrium point `r₀`
in the CR3BP with mass parameter `μ`.

Integration is performed for a normalised time span of `t_f`. For a stable
manifold (Re(λ) < 0) the trajectory is integrated backwards in time. Two
branches are returned, corresponding to positive and negative perturbations
along the eigenvector direction. Each branch is a 4xN array [x; y; vx; vy].
"""
function compute_classical_manifold(
    r₀::Vector{Float64}, μ::Float64,
    λ::ComplexF64, v::Vector{ComplexF64},
    t_f::Float64
)::Tuple{Array{Float64},Array{Float64}}
    stable = real(λ) < 0
    δ₀ = 1e-5 * real(v) / norm(real(v))

    t_f_signed = stable ? -t_f : t_f

    prob_1 = ODEProblem(cr3bp!, [r₀; 0; 0] + δ₀, (0.0, t_f_signed), p=[μ])
    sol_1 = solve(prob_1, Tsit5(), p=[μ], saveat=range(0.0, t_f_signed, length=5000), reltol=1e-12, abstol=1e-12)

    prob_2 = ODEProblem(cr3bp!, [r₀; 0; 0] - δ₀, (0.0, t_f_signed), p=[μ])
    sol_2 = solve(prob_2, Tsit5(), p=[μ], saveat=range(0.0, t_f_signed, length=5000), reltol=1e-12, abstol=1e-12)

    return Array(sol_1), Array(sol_2)
end

"""
    compute_deltav_mars_insertion(sol, r_MPO, SSP) -> Tuple{Float64, Vector{Float64}}

Estimate the delta-V (m/s) required for Mars Parking Orbit (MPO) insertion at
the final state of a manifold trajectory `sol`.

`r_MPO` is the target parking orbit radius in normalised Sun-Mars units.
`SSP` is the solar system parameter dictionary returned by
`load_solar_system_parameters`. The rotating-frame velocity at the trajectory
endpoint is compared against the local circular orbit velocity to compute the
manoeuvre magnitude. Returns (Δv, v_orbit) where v_orbit is the target
circular velocity vector in m/s.
"""
function compute_deltav_mars_insertion(
    sol::AbstractMatrix{Float64},
    r_MPO::Float64,
    SSP::Dict{String,Any}
)::Tuple{Float64,Vector{Float64}}
    x_cross = sol[1, end]
    y_cross = sol[2, end]
    θ_cross = atan(y_cross, x_cross - (1 - SSP["μ_SM"]))

    v_circular = sqrt(SSP["μ_SM"] / r_MPO) # Nondimensional circular velocity at Mars parking orbit

    v_orbit = v_circular * [-sin(θ_cross), cos(θ_cross)]
    # v_orbit_inertial = [v_orbit[1] - y_cross, v_orbit[2] + x_cross]

    v_manifold_inertial = [sol[3, end] - y_cross, sol[4, end] + x_cross - (1 - SSP["μ_SM"])]
    # v_manifold_inertial = [sol[3, end], sol[4, end]]

    # Δv = sqrt((sol[3, end] - v_orbit[1])^2 + (sol[4, end] - v_orbit[2])^2) * (SSP["LU_SM"] / SSP["TU_SM"])
    Δv = sqrt((v_manifold_inertial[1] - v_orbit[1])^2 + (v_manifold_inertial[2] - v_orbit[2])^2) * (SSP["LU_SM"] / SSP["TU_SM"])
    return Δv, v_orbit * (SSP["LU_SM"] / SSP["TU_SM"])
end

"""
    compute_deltav_earth_departure(sol, SSP) -> Tuple{Float64, Vector{Float64}}

Estimate the delta-V (m/s) required for Earth departure at the final state
of a manifold trajectory `sol` expressed in the Sun-Mars CR3BP rotating frame.

`SSP` is the solar system parameter dictionary returned by
`load_solar_system_parameters`. The rotating-frame velocity is converted to
an inertial frame velocity and compared against the Earth circular orbit
velocity at the trajectory endpoint. Returns (Δv, v_orbit) where v_orbit is
the target inertial circular velocity vector in m/s.
"""
function compute_deltav_earth_departure(
    sol::AbstractMatrix{Float64},
    SSP::Dict{String,Any}
)::Tuple{Float64,Vector{Float64}}
    x_cross = sol[1, end]
    y_cross = sol[2, end]
    vx_cross = sol[3, end]
    vy_cross = sol[4, end]

    θ_cross = atan(y_cross, x_cross)

    # Inertial circular orbit velocity 
    v_orbit = [-sin(θ_cross), cos(θ_cross)] * (SSP["LU_SE"] / SSP["TU_SE"])

    # Convert SM-CR3BP rotating frame velocity to inertial
    vx_inertial = (vx_cross - y_cross) * (SSP["LU_SM"] / SSP["TU_SM"])
    vy_inertial = (vy_cross + x_cross) * (SSP["LU_SM"] / SSP["TU_SM"])

    Δv = sqrt((vx_inertial - v_orbit[1])^2 + (vy_inertial - v_orbit[2])^2)
    return Δv, v_orbit
end