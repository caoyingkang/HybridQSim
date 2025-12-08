# Definitions for the 9-qubit [[9,1,3]] surface code.

export SurfaceCode_9_1_3

struct SurfaceCode_9_1_3 <: AbstractCode
    num_data_qubits::Int
    num_logi_qubits::Int
    SurfaceCode_9_1_3() = new(9, 1)
end


function stabilizer_generators(::SurfaceCode_9_1_3)
    stabs = [
        I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I,
        I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I,
        X ⊗ X ⊗ I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I ⊗ X ⊗ X, # plaquette stabilizer
        Z ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ Z, # boundary stabilizer
        I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I,
        I ⊗ I ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I
    ]
    return stabs
end


function stabilizer_matrix(::SurfaceCode_9_1_3)
    Hx = [
        0 0 0 0 0 0 1 1 0;
        0 1 1 0 0 0 0 0 0;
        1 1 0 1 1 0 0 0 0;
        0 0 0 0 1 1 0 1 1
    ]
    Hz = [
        1 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 1;
        0 1 1 0 1 1 0 0 0;
        0 0 0 1 1 0 1 1 0
    ]
    Mx = SMatrix{8,9,Int}([Hx; zero(Hz)])
    Mz = SMatrix{8,9,Int}([zero(Hx); Hz])
    return Mx, Mz
end


function logical_Zs(::SurfaceCode_9_1_3)
    return [Z ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I]
end


function logical_Xs(::SurfaceCode_9_1_3)
    return [X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I]
end


function encoding_isometry(::SurfaceCode_9_1_3)
    ψ0 = 1 / 4 * ((@ket_from_bits "000000000") + (@ket_from_bits "000000110") + (@ket_from_bits "011000000") +
                  (@ket_from_bits "110110000") + (@ket_from_bits "000011011") + (@ket_from_bits "011000110") +
                  (@ket_from_bits "110110110") + (@ket_from_bits "000011101") + (@ket_from_bits "101110000") +
                  (@ket_from_bits "011011011") + (@ket_from_bits "110101011") + (@ket_from_bits "101110110") +
                  (@ket_from_bits "011011101") + (@ket_from_bits "110101101") + (@ket_from_bits "101101011") +
                  (@ket_from_bits "101101101"))
    ψ1 = (X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I ⊗ X ⊗ I ⊗ I) * ψ0
    Uenc_data = [ψ0.data ψ1.data]
    Uenc_dims = GeneralDimensions((ntuple(_ -> 2, Val(9)), ntuple(_ -> 2, Val(1))))
    Uenc = QuantumObject(Uenc_data, dims=Uenc_dims)
    return Uenc
end


function penalty_hamiltonian(::SurfaceCode_9_1_3)
    Hpen = -1. * (I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I +
                  I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I +
                  X ⊗ X ⊗ I ⊗ X ⊗ X ⊗ I ⊗ I ⊗ I ⊗ I +
                  I ⊗ I ⊗ I ⊗ I ⊗ X ⊗ X ⊗ I ⊗ X ⊗ X +
                  Z ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ I ⊗ I ⊗ I +
                  I ⊗ I ⊗ I ⊗ I ⊗ I ⊗ Z ⊗ I ⊗ I ⊗ Z +
                  I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ I ⊗ I +
                  I ⊗ I ⊗ I ⊗ Z ⊗ Z ⊗ I ⊗ Z ⊗ Z ⊗ I)
    return Hpen
end
