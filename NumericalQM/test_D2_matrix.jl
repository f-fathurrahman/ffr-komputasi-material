include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")

function test_D2_3pt()
    D2 = build_D2_matrix_3pt(5, 1.0)
    display(D2)
    println()
end

function test_D2_5pt()
    D2 = build_D2_matrix_5pt(9, 1.0)
    display(D2)
    println()
end

function test_D2_7pt()
    D2 = build_D2_matrix_7pt(9, 1.0)
    display(D2)
    println()
end

test_D2_3pt()
test_D2_5pt()
test_D2_7pt()