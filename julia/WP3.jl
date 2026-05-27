# WORK PACKAGE 3

using Printf

include("Common.jl")

begin
    @printf("=== WORK PACKAGE 3: STABILITY OF THE AEP ===\n\n")

    # Solar System parameters
    SSP = load_solar_system_parameters()
    μ_SE = SSP["μ_SE"]
    μ_SM = SSP["μ_SM"]

    # Compute classical L2 points
    L2_SE = classical_lagrange_point(μ_SE, 2)
    L2_SM = classical_lagrange_point(μ_SM, 2)

    # Desired shift sunward from L2 (normalised units)
    Δx = 0.00050

    sail_SE = eq_condition_solar_sail(L2_SE - [Δx, 0.0], μ_SE)
    sail_SM = eq_condition_solar_sail(L2_SM - [Δx, 0.0], μ_SM)

    λ_classical_SE, _ = linear_stability_analysis(L2_SE, μ_SE)
    λ_classical_SM, _ = linear_stability_analysis(L2_SM, μ_SM)
    λ_AEP_SE, _ = linear_stability_analysis(L2_SE - [Δx, 0.0], μ_SE, (sail_SE[1], sail_SE[2]))
    λ_AEP_SM, _ = linear_stability_analysis(L2_SM - [Δx, 0.0], μ_SM, (sail_SM[1], sail_SM[2]))

    fmt(z) = @sprintf("%+9.6f%+9.6fim", real(z), imag(z))

    println("Eigenvalues (λ) for classical and artificial equilibrium points:\n")
    println("Sun-Earth Classical: [", join(fmt.(λ_classical_SE), ", "), "]")
    println("Sun-Earth Artificial: [", join(fmt.(λ_AEP_SE), ", "), "]")
    println("Sun-Mars Classical: [", join(fmt.(λ_classical_SM), ", "), "]")
    println("Sun-Mars Artificial: [", join(fmt.(λ_AEP_SM), ", "), "]\n")

    # Validation
    println("VALIDATION:\n")
    # Classical
    μ_val = 3e-6
    x_L2_val = classical_lagrange_point(μ_val, 2)
    λ_L2_val, _ = linear_stability_analysis(x_L2_val, μ_val)
    println("Classical L2: [", join(fmt.(λ_L2_val), ", "), "]")
    println("Value according to Zebehely: λ = ±2.4844225586, s = ±i2.0570784929\n")
    μ_val = 3.0034806e-6
    x_AEP_val, n_AEP_val = collinear_AEP(0.01, μ_val, 2)
    solar_sail_val = eq_condition_solar_sail([x_AEP_val, 0.0], μ_val)
    println("Solar Sail AEP: β = ", solar_sail_val[1], ", n = ", solar_sail_val[2])
    println("Solar Sail AEP: x_SL1 = ", x_AEP_val)
    λ_AEP_val, _ = linear_stability_analysis([x_AEP_val, 0.0], μ_val, solar_sail_val)

    println("Solar Sail AEP: [", join(fmt.(λ_AEP_val), ", "), "]")
    @printf("\n=== END OF WP3 ===\n\n")
end