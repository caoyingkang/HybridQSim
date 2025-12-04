# Definitions for the Bacon-Shor code.

export BaconShorCode

struct BaconShorCode <: AbstractCode
    num_data_qubits::Int
    num_logi_qubits::Int
    BaconShorCode() = new(9, 1)
end


"""
Qubit index at the specified row and column.
"""
idx(::BaconShorCode, row::Int, col::Int) = (row - 1) * 3 + col


function stabilizer_generators(::BaconShorCode)
    stabs = [
        X ⊗ X ⊗ X ⊗ X ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ X ⊗ X ⊗ X ⊗ X,
        Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I,
        I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z,
        # TODO
    ]
    return stabs
end


function stabilizer_matrix(::BaconShorCode)
    # TODO
end


function logical_Zs(::BaconShorCode)
    # TODO
end


function logical_Xs(::BaconShorCode)
    # TODO
end


function encoding_isometry(::BaconShorCode)
    # TODO
end


function penalty_hamiltonian(code::Code_5_1_3)
    Hpen = (-1.) * sum(
        multisite_operator(Val(9), idx(code, r, c) => Z, idx(code, r, c + 1) => Z)
        for r in 1:3 for c in 1:2
    )
    Hpen += (-1.) * sum(
        multisite_operator(Val(9), idx(code, r, c) => X, idx(code, r + 1, c) => X)
        for r in 1:2 for c in 1:3
    )
    return Hpen
end
