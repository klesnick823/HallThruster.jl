using HallThruster: HallThruster as het
using Test

function test_landmark_regression()
    case1 = (;
        file = "LANDMARK case 1",
        CFL = 0.5,
        landmark_case = 1,
        thrust = 99.904,
        current = 8.251,
        ion_current = 3.703,
        max_Te = 27.347,
        max_E = 3.931e+4,
        max_nn = 3.888e+19,
        max_ni = 1.186e+18,
        efficiencies = Dict(
            "Mass" => 0.997,
            "Current" => 0.454,
            "Divergence" => 1.0,
            "Voltage" => 0.906,
            "Anode" => 0.417,
        ),
    )

    case2 = (;
        file = "LANDMARK case 2",
        CFL = 0.799,
        landmark_case = 2,
        thrust = 101.136,
        current = 8.085,
        ion_current = 3.673,
        max_Te = 34.06,
        max_E = 4.144e+4,
        max_nn = 5.675e+19,
        max_ni = 4.386e+18,
        efficiencies = Dict(
            "Mass" => 0.948,
            "Current" => 0.456,
            "Divergence" => 1.0,
            "Voltage" => 0.985,
            "Anode" => 0.4,
        ),
    )
    case3 = (;
        file = "LANDMARK case 3",
        CFL = 0.799,
        landmark_case = 3,
        thrust = 99.985,
        current = 7.916,
        ion_current = 3.683,
        max_Te = 35.473,
        max_E = 4.17e+4,
        max_nn = 6.573e+19,
        max_ni = 5.134e+18,
        efficiencies = Dict(
            "Mass" => 0.935,
            "Current" => 0.466,
            "Divergence" => 1.0,
            "Voltage" => 0.987,
            "Anode" => 0.392,
        ),
    )

    return @testset "LANDMARK regression" begin
        check_regression_case(case1)
        check_regression_case(case2)
        check_regression_case(case3)
    end
end
