# WORK PACKAGE 4.1

using Printf, Julianim, MathTeXEngine

include("Common.jl")

begin
    @printf("=== WORK PACKAGE 4.1: INVARIANT MANIFOLDS ===\n\n")

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

    # Compute classical L2 points
    L₂ = classical_lagrange_point(μ, 2)

    λ, v = linear_stability_analysis(L₂, μ)

    # Compute stable and unstable manifolds
    sol_stable_1, sol_stable_2 = compute_classical_manifold(
        L₂, μ,
        λ[1],
        v[:, 1],
        5 * 365.256 * 24 * 3600 / SSP["TU_SM"])

    sol_unstable_1, sol_unstable_2 = compute_classical_manifold(
        L₂, μ,
        λ[4],
        v[:, 4],
        5 * 365.256 * 24 * 3600 / SSP["TU_SM"])

    # Plotting    
    fig = Figure(size=(1500, 800))
    ax1 = Axis(fig[1, 1];
        xlabel=L"x \ \mathrm{[LU]}",
        ylabel=L"y \ \mathrm{[LU]}",
        aspect=DataAspect(),
        limits=(-0.1, 1.1, -0.6, 0.6)
    )

    l_stable_1 = lines!(ax1,
        sol_stable_1[1, :], sol_stable_1[2, :];
        color=BLUE)
    l_unstable_1 = lines!(ax1,
        sol_unstable_1[1, :], sol_unstable_1[2, :];
        color=RED)

    l_stable_2 = lines!(ax1,
        sol_stable_2[1, :], sol_stable_2[2, :];
        color=BLUE)
    l_unstable_2 = lines!(ax1,
        sol_unstable_2[1, :], sol_unstable_2[2, :];
        color=RED)

    sc_sun = scatter!(ax1, [-μ], [0.0];
        color=ORANGE,
        strokecolor=:black,
        marker=:circle,
        markersize=14)
    sc_mars = scatter!(ax1, [1 - μ], [0.0];
        color=GREEN,
        strokecolor=:black,
        marker=:circle,
        markersize=14)
    leg = axislegend(ax1,
        [l_stable_1, l_unstable_1, sc_sun, sc_mars],
        [L"Stable Manifold $W_s$", L"Unstable Manifold $W_u$", "Sun", "Mars"];
        position=:lt)

    ax2 = Axis(fig[1, 2];
        xlabel=L"x \ \mathrm{[LU]}",
        ylabel=L"y \ \mathrm{[LU]}",
        aspect=DataAspect(),
        limits=(0.9925, 1.0075, -0.0075, 0.0075)
    )

    lagrange_points = classical_lagrange_points(μ)
    scatter!(ax2,
        lagrange_points[1, 1:2],
        lagrange_points[2, 1:2];
        color=:black,
        markersize=10,
        marker=:cross,
        strokecolor=:black,
        strokewidth=1.5,
        markerspace=:pixel)

    labels = [L"L$_1$", L"L$_2$"]
    offsets = [(-25, 6), (-15, 6)]
    for i in 1:2
        text!(ax2, lagrange_points[1, i], lagrange_points[2, i];
            text=labels[i],
            offset=offsets[i])
    end

    l_unstable_1 = lines!(ax2,
        sol_unstable_1[1, :], sol_unstable_1[2, :];
        color=RED,
        linestyle=:dash)
    l_unstable_2 = lines!(ax2,
        sol_unstable_2[1, :], sol_unstable_2[2, :];
        color=RED)
    l_stable_1 = lines!(ax2,
        sol_stable_1[1, :], sol_stable_1[2, :];
        color=BLUE)
    l_stable_2 = lines!(ax2,
        sol_stable_2[1, :], sol_stable_2[2, :];
        color=BLUE,
        linestyle=:dash)

    sc_mars = scatter!(ax2, [1 - μ], [0.0];
        color=GREEN,
        strokecolor=:black,
        marker=:circle,
        markersize=14)

    mars_parking_orbit_radius = 306666e3 / SSP["LU_SM"]
    θ = range(0, 2π, length=300)
    l_mpo = lines!(ax2,
        (1 - μ) .+ mars_parking_orbit_radius .* cos.(θ),
        mars_parking_orbit_radius .* sin.(θ);
        color=:black,
        linestyle=:dash)

    leg = axislegend(ax2,
        [l_stable_1, l_stable_2, l_unstable_1, l_unstable_2, l_mpo, sc_mars],
        [L"$W^+_s$", L"$W^-_s$", L"$W^+_u$", L"$W^-_u$", "MPO", "Mars"];
        position=:lt, nbanks=2, backgroundcolor=:white)

    save(joinpath(figures_dir, "wp41_manifolds.png"), fig, px_per_unit=4)

    @printf("\n=== END OF WP4.1 ===\n\n")
end