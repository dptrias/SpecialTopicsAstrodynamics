# WORK PACKAGE 3

using Printf, Julianim, MathTeXEngine, StaticArrays, DifferentialEquations

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

    @printf("=== WORK PACKAGE 3: STABILITY OF THE AEP ===\n\n")
    # Output directory
    script_dir = @__DIR__
    figures_dir = joinpath(script_dir, "..", "figures")
    mkpath(figures_dir)

    # Solar System parameters
    SSP = load_solar_system_parameters()
    μ_SE = SSP["μ_SE"]
    μ_SM = SSP["μ_SM"]


    @printf("\n=== END OF WP3 ===\n\n")
end