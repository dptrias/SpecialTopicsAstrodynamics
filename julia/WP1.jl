# WORK PACKAGE 1
using Printf

include("Common.jl")

begin
    @printf("=== WORK PACKAGE 1: CLASSICAL LAGRANGE POINTS ===\n\n")

    # Solar System parameters
    SSP = load_solar_system_parameters()

    # Compute L2 points
    L2_SE = classical_lagrange_point(SSP["μ_SE"], 2)
    L2_SM = classical_lagrange_point(SSP["μ_SM"], 2)

    # JPL reference values
    L2_SE_JPL = 1.01009044
    L2_SM_JPL = 1.00476311

    @printf("Sun-Earth L2: %.8f  (JPL: %.8f)\n", L2_SE, L2_SE_JPL)
    @printf("Sun-Mars  L2: %.8f  (JPL: %.8f)\n", L2_SM, L2_SM_JPL)

    @printf("Difference for Sun-Earth L2: %.8e\n", abs(L2_SE - L2_SE_JPL))
    @printf("Difference for Sun-Mars  L2: %.8e\n", abs(L2_SM - L2_SM_JPL))

    @printf("\n=== END OF WP1 ===\n\n")
end