# WORK PACKAGE 4.2

using Printf, Julianim, MathTeXEngine, DifferentialEquations

include("Common.jl")

begin
    @printf("=== WORK PACKAGE 4.2: CLASSICAL TRANSFERS ===\n\n")

    # Set plotting theme
    set_publication_theme!()
    update_theme!(
        fonts=Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        ))

    # Output directory
    script_dir = @__DIR__
    figures_dir = joinpath(script_dir, "..", "figures")
    mkpath(figures_dir)

    # Solar System parameters
    SSP = load_solar_system_parameters()
    μ = SSP["μ_SM"]

    r_Earth_SM = SSP["LU_SE"] / SSP["LU_SM"]  # Earth-Sun distance in SM normalised units
    r_MPO_meters = 306666e3 # 182222.0e3 # Mars parking orbit radius in meters
    r_MPO = r_MPO_meters / SSP["LU_SM"]  # Mars parking orbit radius in SM normalised units

    # Compute classical L2 points 
    L₂ = classical_lagrange_point(μ, 2)

    λ, v = linear_stability_analysis(L₂, μ)

    # SM-L2 to Mars parking orbit transfer
    condition(u, t, integrator) = (u[1] - (1 - integrator.p[1]))^2 + u[2]^2 - r_MPO^2
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    prob = ODEProblem(
        cr3bp!,
        [L₂; 0; 0] - 1e-5 * real(v)[:, 4] / norm(real(v[:, 4])),
        (0.0, 5 * 365.256 * 24 * 3600 / SSP["TU_SM"]),
        [μ])
    sol_l2_mpo = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, callback=cb)

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"x \ \mathrm{[LU]}",
        ylabel=L"y \ \mathrm{[LU]}",
        aspect=DataAspect(),
        limits=(0.998, 1.006, -0.004, 0.004)
    )

    l_transfer = lines!(ax,
        sol_l2_mpo[1, :], sol_l2_mpo[2, :];
        color=BLUE)

    sc_l2 = scatter!(ax,
        [L₂[1]], [L₂[2]];
        label=L"L$_2$",
        marker=:cross,
        strokecolor=:black,
        color=:black)

    sc_mars = scatter!(ax,
        [1 - μ], [0];
        color=GREEN, marker=:circle, label="Mars")

    θ = range(0, 2π, length=300)
    l_mpo = lines!(ax,
        r_MPO * cos.(θ) .+ (1 - μ), r_MPO * sin.(θ);
        color=:black, linestyle=:dash, label="Mars Parking Orbit")

    sc_intersection = scatter!(ax,
        [sol_l2_mpo[1, end]], [sol_l2_mpo[2, end]];
        color=BLUE, label="Intersection")

    axislegend(ax; position=:rt, backgroundcolor=:white)

    save(joinpath(figures_dir, "wp42_mars_transfer.png"), fig, px_per_unit=4)

    t_flight = sol_l2_mpo.t[end] * SSP["TU_SM"] / (3600 * 24)  # in days
    @printf("Time of flight from SM-L2 to Mars parking orbit: %.2f days\n", t_flight)

    # Delta-V
    Δv_Mars, _ = compute_deltav_mars_insertion(sol_l2_mpo, r_MPO, SSP)
    @printf("Delta-V required at intersection for circular orbit insertion: %.4f m/s\n", Δv_Mars)

    @printf("\n=== END OF WP4.2 ===\n\n")
end


