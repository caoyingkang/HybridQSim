# Definitions for the 5-qubit [[5,1,3]] code.

export Code_5_1_3

struct Code_5_1_3 <: AbstractCode
    num_data_qubits::Int
    num_logi_qubits::Int
    Code_5_1_3() = new(5, 1)
end


function stabilizer_generators(::Code_5_1_3)
    stabs = [
        X ⊗ Z ⊗ Z ⊗ X ⊗ I,
        I ⊗ X ⊗ Z ⊗ Z ⊗ X,
        X ⊗ I ⊗ X ⊗ Z ⊗ Z,
        Z ⊗ X ⊗ I ⊗ X ⊗ Z
    ]
    return stabs
end


function stabilizer_matrix(::Code_5_1_3)
    Mx = @SMatrix [
        1 0 0 1 0;
        0 1 0 0 1;
        1 0 1 0 0;
        0 1 0 1 0
    ]
    Mz = @SMatrix [
        0 1 1 0 0;
        0 0 1 1 0;
        0 0 0 1 1;
        1 0 0 0 1
    ]
    return Mx, Mz
end


function logical_Zs(::Code_5_1_3)
    return [Z ⊗ Z ⊗ Z ⊗ Z ⊗ Z]
end


function logical_Xs(::Code_5_1_3)
    return [X ⊗ X ⊗ X ⊗ X ⊗ X]
end


function encoding_isometry(::Code_5_1_3)
    ψ0 = 1 / 4 * ((@ket_from_bits "00000") + (@ket_from_bits "10010") + (@ket_from_bits "01001") +
                  (@ket_from_bits "10100") + (@ket_from_bits "01010") - (@ket_from_bits "11011") -
                  (@ket_from_bits "00110") - (@ket_from_bits "11000") - (@ket_from_bits "11101") -
                  (@ket_from_bits "00011") - (@ket_from_bits "11110") - (@ket_from_bits "01111") -
                  (@ket_from_bits "10001") - (@ket_from_bits "01100") - (@ket_from_bits "10111") +
                  (@ket_from_bits "00101"))
    ψ1 = (X ⊗ X ⊗ X ⊗ X ⊗ X) * ψ0
    Uenc_data = [ψ0.data ψ1.data]
    Uenc_dims = GeneralDimensions((ntuple(_ -> 2, Val(5)), ntuple(_ -> 2, Val(1))))
    Uenc = QuantumObject(Uenc_data, dims=Uenc_dims)
    return Uenc
end


function penalty_hamiltonian(::Code_5_1_3)
    Hpen = -1. * (X ⊗ Z ⊗ Z ⊗ X ⊗ I +
                  I ⊗ X ⊗ Z ⊗ Z ⊗ X +
                  X ⊗ I ⊗ X ⊗ Z ⊗ Z +
                  Z ⊗ X ⊗ I ⊗ X ⊗ Z +
                  Z ⊗ Z ⊗ X ⊗ I ⊗ X)
    return Hpen
end
