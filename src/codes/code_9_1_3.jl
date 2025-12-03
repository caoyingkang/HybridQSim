# Definitions for the 9-qubit [[9,1,3]] surface code.

export Code_9_1_3

struct Code_9_1_3 <: AbstractCode
    num_data_qubits::Int
    num_logi_qubits::Int
    Code_9_1_3() = new(9, 1)
end


function stabilizer_generators(::Code_9_1_3)
    stabs = [
        Z ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I,
        I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ Z, # boundary stabilizer
        X ⊗ X ⊗ I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I,
        I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I,
        I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I ⊗ X ⊗ X  # plaquette stabilizer
    ]
    return stabs
end


function stabilizer_matrix(::Code_9_1_3)
    Mx = @SMatrix [
        0 0 0 0 0 0 1 1 0;
        0 1 1 0 0 0 0 0 0;
        1 1 0 1 1 0 0 0 0;
        0 0 0 0 1 1 0 1 1
    ]
    Mz = @SMatrix [
        1 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 1;
        0 1 1 0 1 1 0 0 0;
        0 0 0 1 1 0 1 1 0
    ]
    return Mx, Mz
end


function logical_Zs(::Code_9_1_3)
    return [Z ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I]
end


function logical_Xs(::Code_9_1_3)
    return [I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X]
end


function encoding_isometry(::Code_9_1_3) 
    ψ0 = 1 / 4 * ((@ket_from_bits "000000000") + (@ket_from_bits "000001110") + (@ket_from_bits "000011011") +
                  (@ket_from_bits "000010101") + (@ket_from_bits "011000000") + (@ket_from_bits "011001110") +
                  (@ket_from_bits "011011011") + (@ket_from_bits "011010101") + (@ket_from_bits "110100000") +
                  (@ket_from_bits "110101110") + (@ket_from_bits "110111011") + (@ket_from_bits "110110101") +
                  (@ket_from_bits "101100000") + (@ket_from_bits "101101110") + (@ket_from_bits "101111011") +
                  (@ket_from_bits "101110101"))
    ψ1 = (I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X) * ψ0
    Uenc_data = [ψ0.data ψ1.data]
    Uenc_dims = GeneralDimensions((ntuple(_ -> 2, Val(9)), ntuple(_ -> 2, Val(1))))
    Uenc = QuantumObject(Uenc_data, dims=Uenc_dims)
    return Uenc
end


function penalty_hamiltonian(::Code_9_1_3)
    Hpen =( -1. * (Z ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I +
                I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I +
                I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I +
                I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ Z + 
                X ⊗ X ⊗ I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I +
                I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I +
                I ⊗ I ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I +
                I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I ⊗ X ⊗ X ) )
    return Hpen
end
