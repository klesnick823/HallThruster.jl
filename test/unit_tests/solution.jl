using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_solution_serialization()

    ncells = 50

    config = het.Config(;
        thruster = het.SPT_100,
        discharge_voltage = 300.0,
        domain = (0.0, 0.8),
        anode_mass_flow_rate = 5.0e-6,
    )

    simparams = het.SimParams(
        grid = het.UnevenGrid(ncells),
        duration = 1.0e-6,
        num_save = 1000,
        verbose = false
    )

    cfg_ser = het.serialize(config)

    sol = het.run_simulation(config, simparams)
    avg = het.time_average(sol) 

    frame = avg.frames[1]
    neutral_state = frame.neutrals[:Xe]

    neu = het.serialize(neutral_state)
    test_roundtrip(het.SpeciesState, neu)

    ions = frame.ions[:Xe]
    @test typeof(ions) == Vector{het.SpeciesState}
    test_roundtrip(typeof(ions), ions)

    ion_dict = frame.ions
    @test typeof(ion_dict) == het.OrderedDict{Symbol, Vector{het.SpeciesState}}

    ion_dict_ser = het.serialize(frame.ions)
    @test typeof(ion_dict_ser) == het.OrderedDict{String, Vector{het.OrderedDict{String, Any}}}

    ion_dict_2 = het.deserialize(typeof(ion_dict), ion_dict_ser)
    for k in keys(ion_dict)
        for i in eachindex(ion_dict[k])
            @test struct_eq(ion_dict[k][i], ion_dict_2[k][i])
        end
    end

    test_roundtrip(het.Frame, frame)
    test_roundtrip(het.Solution, avg)

    return
end

test_solution_serialization()
