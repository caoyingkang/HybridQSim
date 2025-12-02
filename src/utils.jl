const X = sigmax()
const Y = sigmay()
const Z = sigmaz()
const I = one(sigmax())
const PauliType = typeof(sigmax())

const KET0 = basis(2, 0)
const KET1 = basis(2, 1)


"""
Given a binary tuple like (1,0,1,1), return the computational basis state |1011⟩.
"""
function ket_from_bits(bits::NTuple{N,Int}) where {N}
    @inbounds return foldl(⊗, ntuple(i -> (bits[i] == 0 ? KET0 : KET1), Val(N)))
end


"""
A macro version of the function `ket_from_bits` for convenient usage.
For example, `@ket_from_bits "0101"` expands to `ket_from_bits((0,1,0,1))`.
"""
macro ket_from_bits(bitstring)
    s = String(bitstring)
    bits = Tuple(c == '0' ? 0 : 1 for c in s)
    return :(ket_from_bits($bits))
end


"""
Given an array of stabilizer generators, return the projector onto the code space.
"""
function get_Penc(stabs::Vector{OT}) where {OT<:QuantumObject{Operator}}
    projectors = [(1 + s) / 2 for s in stabs]
    return foldl(*, projectors)
end


const LUTKeyType{LenSynd} = NTuple{LenSynd,Int}
const LUTValueType = Pair{Int,PauliType}
const LUTType{LenSynd} = Dict{LUTKeyType{LenSynd},LUTValueType}


"""
Given the stabilizer matrix in binary symplectic form (Mx, Mz), return the decoding 
lookup table as a dictionary mapping syndromes to correction operations.
Only single-qubit Pauli errors are considered currently.
"""
function get_decoding_lut(
    Mx::SMatrix{NumStabs,NumDataQubits,Int},
    Mz::SMatrix{NumStabs,NumDataQubits,Int}
) where {NumStabs,NumDataQubits}
    lut = LUTType{NumStabs}()
    for i in 1:NumDataQubits
        # X error on qubit i
        sx = Mz[:, i]
        sx_tuple = Tuple(sx)
        @assert !haskey(lut, sx_tuple)
        lut[sx_tuple] = (i => X)

        # Z error on qubit i
        sz = Mx[:, i]
        sz_tuple = Tuple(sz)
        @assert !haskey(lut, sz_tuple)
        lut[sz_tuple] = (i => Z)

        # Y error on qubit i
        sy = (sx + sz) .% 2
        sy_tuple = Tuple(sy)
        @assert !haskey(lut, sy_tuple)
        lut[sy_tuple] = (i => Y)
    end
    return lut
end


"""
Return an array of Lindblad operators for T1 and T2 processes on all of the N qubits.
"""
function get_lindblad_ops(T1::Float64, T2::Float64, ::Val{N}) where {N}
    @assert T2 < 2T1 "T2 must be less than 2 * T1"

    γ1 = 1 / T1 # relaxation rate
    γϕ = 1 / T2 - 1 / (2T1) # pure dephasing rate
    l1 = sqrt(γ1) * sigmam()
    lϕ = sqrt(γϕ / 2) * sigmaz()
    lindblad_ops = [multisite_operator(Val(N), i => l) for i in 1:N for l in (l1, lϕ)]
    return lindblad_ops
end


"""
A CPTP map implementing digital feedback.
"""
struct FeedbackCPTPMap{OT<:QuantumObject}
    kraus_ops::Vector{OT}
end


"""
Constructor for FeedbackCPTPMap from an array of stabilizer generators and a decoding lookup table.
"""
function FeedbackCPTPMap(
    stabs::Vector{OT},
    lut::LUTType{LenSynd},
    ::Val{NumDataQubits}
) where {OT<:QuantumObject{Operator},LenSynd,NumDataQubits}
    Π0 = [(1 + s) / 2 for s in stabs]
    Π1 = [(1 - s) / 2 for s in stabs]
    kraus_ops = OT[]
    for synd in Iterators.product(ntuple(_ -> (0, 1), Val(LenSynd))...)
        P = one(Π0[1])
        for i in 1:LenSynd
            P *= synd[i] == 0 ? Π0[i] : Π1[i]
        end
        if haskey(lut, synd)
            C = multisite_operator(Val(NumDataQubits), lut[synd])
            push!(kraus_ops, C * P)
        else
            push!(kraus_ops, P)
        end
    end
    return FeedbackCPTPMap(kraus_ops)
end


"""
Apply the CPTP map to a density matrix ρ0.
"""
function (f::FeedbackCPTPMap)(ρ0)
    ρ = zero(ρ0)
    for K in f.kraus_ops
        ρ += K * ρ0 * K'
    end
    return ρ
end