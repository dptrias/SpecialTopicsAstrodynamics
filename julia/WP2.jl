# WORK PACKAGE 2

using Printf, Julianim, MathTeXEngine, StaticArrays

include("Common.jl")

begin
    # Set plotting theme
    set_publication_theme!()
    update_theme!(
        fonts=Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular)
        ))

    @printf("=== WORK PACKAGE 2: ARTIFICIAL EQUILIBRIUM POINTS ===\n\n")
    # Output directory
    script_dir = @__DIR__
    figures_dir = joinpath(script_dir, "..", "figures")
    mkpath(figures_dir)

    # Solar System parameters
    SSP = load_solar_system_parameters()
    μ_SE = SSP["μ_SE"]
    μ_SM = SSP["μ_SM"]

    # Compute classical L2 points
    L2_SE = classical_lagrange_points(μ_SE, 2)
    L2_SM = classical_lagrange_points(μ_SM, 2)

    # Desired shift sunward from L2 (normalised units)
    Δx = 0.00050

    β_SE, _ = eq_condition_solar_sail(SVector(L2_SE - Δx, 0.0, 0.0), μ_SE)
    β_SM, _ = eq_condition_solar_sail(SVector(L2_SM - Δx, 0.0, 0.0), μ_SM)

    @printf("Beta for Sun-Earth L2: %.8e\n", β_SE)
    @printf("Beta for Sun-Mars  L2: %.8e\n", β_SM)

    # Physical constants and system scales
    μ_S = 132712e15          # m^3/s^2  (solar gravitational parameter)
    m_s = 100.0              # kg       (spacecraft mass)
    S_sun = 1361.0             # W/m^2    (solar irradiance at 1 AU)
    c = 299792458.0        # m/s
    LU_SE_m = 149597870.7e3    # m
    LU_SM_m = 208321282.0e3    # m

    x_SE = (L2_SE - Δx) * LU_SE_m
    x_SM = (L2_SM - Δx) * LU_SM_m

    A_SE = β_SE * μ_S * m_s / (2 * solar_radiation_pressure(S_sun, x_SE) * x_SE^2)
    A_SM = β_SM * μ_S * m_s / (2 * solar_radiation_pressure(S_sun, x_SM) * x_SM^2)

    @printf("Required sail side length for Sun-Earth L2: %.4f m\n", sqrt(A_SE))
    @printf("Required sail side length for Sun-Mars  L2: %.4f m\n", sqrt(A_SM))

    # Sensitivity analysis over a range of Lagrange point shifts
    Δx_values = LinRange(Δx, 2 * Δx, 21)
    β_SE_vals = [eq_condition_solar_sail(SVector(L2_SE - Δx, 0.0, 0.0), μ_SE)[1] for Δx in Δx_values]
    β_SM_vals = [eq_condition_solar_sail(SVector(L2_SM - Δx, 0.0, 0.0), μ_SM)[1] for Δx in Δx_values]

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\Delta x \ \mathrm{[LU]}",
        ylabel=L"\beta \ \mathrm{[-]}",
        limits=(Δx, 2 * Δx, nothing, nothing),
        xticklabelrotation=π / 6
    )

    lines!(ax, collect(Δx_values), β_SE_vals; color=:blue, linewidth=4, label=L"Sun-Earth L$_2$")
    lines!(ax, collect(Δx_values), β_SM_vals; color=:red, linewidth=4, label=L"Sun-Mars L$_2$")

    axislegend(ax; position=:lt)

    save(joinpath(figures_dir, "beta_sensitivity.png"), fig, px_per_unit=4)

    # Validation test
    beta_test, _ = eq_condition_solar_sail(SVector(0.983867, 0.0, 0.0), μ_SE)
    @printf("Calculated beta for Sun-Earth AEP at point 1: %.8e\n", beta_test)
    # Validation test
    beta_test, _ = eq_condition_solar_sail(SVector(0.983908, -0.00144, 0.0), μ_SE)
    @printf("Calculated beta for Sun-Earth AEP at point 2: %.8e\n", beta_test)

    @printf("\n=== END OF WP2 ===\n\n")
end