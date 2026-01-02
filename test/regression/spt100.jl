using HallThruster: HallThruster as het
using Test

function test_spt100_regression()
    return @testset "SPT-100 regression" begin
        baseline = (;
            file = "baseline.json",
            thrust = 86.908,
            current = 4.618,
            ion_current = 3.897,
            max_Te = 24.692,
            max_E = 6.505e+4,
            max_nn = 2.088e+19,
            max_ni = 9.627e+17,
            efficiencies = Dict(
                "Mass" => 0.952,
                "Current" => 0.874,
                "Divergence" => 0.949,
                "Voltage" => 0.66,
                "Anode" => 0.569,
            ),
        )
        check_regression_case(baseline)

        with_plume = (;
            file = "with_plume.json",
            thrust = 106.638,
            current = 5.052,
            ion_current = 4.096,
            max_Te = 27.26,
            max_E = 9.749e+4,
            max_nn = 2.096e+19,
            max_ni = 1.174e+18,
            efficiencies = Dict(
                "Mass" => 0.982,
                "Current" => 0.826,
                "Divergence" => 0.963,
                "Voltage" => 0.717,
                "Anode" => 0.568,
            ),
        )
        check_regression_case(with_plume)
    end
end
