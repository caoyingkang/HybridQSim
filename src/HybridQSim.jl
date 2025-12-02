module HybridQSim

export get_decoding_lut, FeedbackCPTPMap, get_lindblad_ops
export stabilizer_generators, stabilizer_matrix, encoding_isometry, penalty_hamiltonian
export Code_5_1_3
export run_simulation

using QuantumToolbox: sigmax, sigmay, sigmaz, sigmam, basis, ⊗, QuantumObject, Operator, Ket,
    GeneralDimensions, multisite_operator, mesolve, fidelity
using StaticArrays: SMatrix, @SMatrix

# Utility functions, types, and constants
include("utils.jl")

# Code definitions
abstract type AbstractCode end
include("codes/code_5_1_3.jl")
export Code_5_1_3

# Main functions for simulation
function run_simulation(;
    H::OT, ψ0::KT, tlist::AbstractVector{<:Real}, c_ops::Vector{OT}, feedback::FeedbackCPTPMap{OT}
) where {OT<:QuantumObject{Operator},KT<:QuantumObject{Ket}}
    sol = mesolve(H, ψ0, tlist, c_ops)
    states = sol.states
    fids_before_feedback = [fidelity(ρ, ψ0) for ρ in states]
    fids_after_feedback = [fidelity(feedback(ρ), ψ0) for ρ in states]
    return fids_before_feedback, fids_after_feedback
end

end # module HybridQSim
