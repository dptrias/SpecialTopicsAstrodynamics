# WORK PACKAGE 4.3

using Printf, Julianim, MathTeXEngine, DifferentialEquations, Dates

include("Common.jl")

begin
    @printf("=== WORK PACKAGE 4.3: TRANSFERS WITH SOLAR SAILS ===\n\n")

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
    r_MPO_meters = 306666e3 #    182222.0e3 # Mars parking orbit radius in meters
    r_MPO = r_MPO_meters / SSP["LU_SM"]  # Mars parking orbit radius in SM normalised units

    # Solar sail parameters
    β = 0.025
    α_range = range(-90, 90, length=73)  # Cone angles, 2.5° increments

    # Compute classical L2 point
    L₂ = classical_lagrange_point(μ, 2)

    λ, v = linear_stability_analysis(L₂, μ)
    u₀_leg_1 = [L₂; 0; 0] - 1e-5 * real(v)[:, 1] / norm(real(v[:, 1])) # Initial state for leg 1 (Earth to SM-L2)
    u₀_leg_2 = [L₂; 0; 0] - 1e-5 * real(v)[:, 4] / norm(real(v[:, 4])) # Initial state for leg 2 (SM-L2 to Mars parking orbit)

    Δv_earth = zeros(length(α_range))
    Δv_mars = zeros(length(α_range))
    sols_leg_1 = Vector{Any}(undef, length(α_range))
    tof_leg_1 = zeros(length(α_range))
    sols_leg_2 = Vector{Any}(undef, length(α_range))
    tof_leg_2 = zeros(length(α_range))
    condition_Earth(u, t, integrator) = (u[1] + integrator.p[1])^2 + u[2]^2 - r_Earth_SM^2
    condition_MPO(u, t, integrator) = (u[1] - (1 - integrator.p[1]))^2 + u[2]^2 - r_MPO^2
    affect!(integrator) = terminate!(integrator)

    cb_Earth = ContinuousCallback(condition_Earth, affect!)
    cb_MPO = ContinuousCallback(condition_MPO, affect!)
    for (α_index, α) in enumerate(α_range)
        # Earth to SM-L2 transfer (Leg 1)
        local prob = ODEProblem(
            cr3bp_ss!,
            u₀_leg_1,
            (0.0, -5 * 365.256 * 24 * 3600 / SSP["TU_SM"]),
            [μ, β, deg2rad(α)])
        local sol = solve(prob, Tsit5(),
            reltol=1e-12, abstol=1e-12, callback=cb_Earth)
        if sol.retcode == ReturnCode.Terminated
            sols_leg_1[α_index] = sol
            tof_leg_1[α_index] = -sol.t[end] * SSP["TU_SM"] / (24 * 3600 * 365.256)
            Δv_earth[α_index], _ = compute_deltav_earth_departure(sol, μ, SSP)
        else
            Δv_earth[α_index] = Inf
            tof_leg_1[α_index] = Inf
        end

        # SM-L2 to Mars parking orbit transfer (Leg 2)
        local prob = ODEProblem(
            cr3bp_ss!,
            u₀_leg_2,
            (0.0, 5 * 365.256 * 24 * 3600 / SSP["TU_SM"]),
            [μ, β, deg2rad(α)])
        local sol = solve(prob, Tsit5(),
            reltol=1e-12, abstol=1e-12, callback=cb_MPO)
        if sol.retcode == ReturnCode.Terminated
            sols_leg_2[α_index] = sol
            tof_leg_2[α_index] = sol.t[end] * SSP["TU_SM"] / (24 * 3600 * 365.256)
            Δv_mars[α_index], _ = compute_deltav_mars_insertion(sol, μ, r_MPO, SSP)
        else
            Δv_mars[α_index] = Inf
            tof_leg_2[α_index] = Inf
        end
    end

    idx_leg_1_min = argmin(Δv_earth)
    idx_leg_2_min = argmin(Δv_mars)

    @printf("Optimal Leg 1 (Earth to SM-L2):\n")
    @printf("   Optimal cone angle: %.1f°\n", α_range[idx_leg_1_min])
    @printf("   Minimum Δv: %.4f m/s\n", Δv_earth[idx_leg_1_min])
    @printf("   Time of flight: %.4f years\n\n", tof_leg_1[idx_leg_1_min])
    @printf("Optimal Leg 2 (SM-L2 to Mars parking orbit):\n")
    @printf("   Optimal cone angle: %.1f°\n", α_range[idx_leg_2_min])
    @printf("   Minimum Δv: %.4f m/s\n", Δv_mars[idx_leg_2_min])
    @printf("   Time of flight: %.4f years\n\n", tof_leg_2[idx_leg_2_min])

    valid_1 = isfinite.(Δv_earth)
    valid_2 = isfinite.(Δv_mars)
    fig = Figure(size=(1500, 800))
    ax1 = Axis(fig[1, 1];
        xlabel=L"x \ \mathrm{[LU]}",
        ylabel=L"y \ \mathrm{[LU]}",
        aspect=DataAspect(),
        limits=(-0.8, 1.2, -1, 1)
    )

    l_leg_1 = lines!(ax1,
        sols_leg_1[idx_leg_1_min][1, :], sols_leg_1[idx_leg_1_min][2, :];
        color=BLUE,
        label="Leg 1")

    sc_sun = scatter!(ax1, [-μ], [0.0];
        color=ORANGE,
        strokecolor=:black,
        marker=:circle,
        label="Sun")
    sc_mars = scatter!(ax1, [1 - μ], [0.0];
        color=GREEN,
        strokecolor=:black,
        marker=:circle,
        label="Mars")

    θ = range(0, 2π, length=300)
    l_earth_orbit = lines!(ax1,
        r_Earth_SM * cos.(θ) .- μ, r_Earth_SM * sin.(θ);
        color=:black, linestyle=:dash, label="Earth Orbit")

    axislegend(ax1; position=:lt)

    ax2 = Axis(fig[1, 2];
        xlabel=L"x \ \mathrm{[LU]}",
        ylabel=L"y \ \mathrm{[LU]}",
        aspect=DataAspect(),
        limits=(0.997, 1.007, -0.005, 0.005)
    )

    sc_mars = scatter!(ax2, [1 - μ], [0.0];
        color=GREEN,
        strokecolor=:black,
        marker=:circle, label="Mars")

    l_mpo = lines!(ax2,
        r_MPO * cos.(θ) .+ (1 - μ), r_MPO * sin.(θ);
        color=:black, linestyle=:dash, label="MPO")

    l_leg_1 = lines!(ax2,
        sols_leg_1[idx_leg_1_min][1, :], sols_leg_1[idx_leg_1_min][2, :];
        color=BLUE, label="Leg 1")

    l_leg_2 = lines!(ax2,
        sols_leg_2[idx_leg_2_min][1, :], sols_leg_2[idx_leg_2_min][2, :];
        color=RED, label="Leg 2")

    sc_l2 = scatter!(ax2,
        [L₂[1]], [L₂[2]];
        label=L"L$_2$",
        marker=:cross,
        strokecolor=:black,
        color=:black)

    axislegend(ax2; position=:lt)

    save(joinpath(figures_dir, "wp43_optimal_transfers.png"), fig, px_per_unit=4)

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\alpha \ \mathrm{[^\circ]}",
        ylabel=L"\Delta v \ \mathrm{[km/s]}",
        limits=(minimum(α_range), maximum(α_range), nothing, nothing),
        xticklabelrotation=π / 6
    )

    lines!(ax, collect(α_range)[valid_1], Δv_earth[valid_1]; color=:blue, linewidth=4, label=L"\Delta V_E")
    lines!(ax, collect(α_range)[valid_2], Δv_mars[valid_2]; color=:red, linewidth=4, label=L"\Delta V_M")
    axislegend(ax; position=:lt)

    save(joinpath(figures_dir, "wp43_deltav_vs_alpha.png"), fig, px_per_unit=4)

    # Earth-Mars angle
    θ_Earth = atan(sols_leg_1[idx_leg_1_min][2, end], sols_leg_1[idx_leg_1_min][1, end] + μ)
    @printf("Earth-Mars angle at Earth departure: %.4f°\n", rad2deg(θ_Earth))
    θ₀ = deg2rad(122.35)
    t_launch = (θ_Earth - θ₀) / (1 / SSP["TU_SE"] - 1 / SSP["TU_SM"])
    initial_epoch = DateTime(2030, 1, 1)
    launch_epoch = initial_epoch + Second(round(Int, t_launch))
    @printf("Earliest possible launch date: %s\n", string(launch_epoch))

    @printf("\n=== END OF WP4.3 ===\n\n")
end