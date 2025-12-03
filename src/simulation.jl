# Main function for simulation

export run_simulation

function run_simulation(
    ; H::OT, ψ0::KT, tlist::AbstractVector{<:Real}, c_ops::Vector{OT}, feedback::FeedbackCPTPMap{OT}
) where {OT<:QuantumObject{Operator},KT<:QuantumObject{Ket}}
    sol = mesolve(H, ψ0, tlist, c_ops)
    states = sol.states
    fids_before_feedback = [fidelity(ρ, ψ0) for ρ in states]
    fids_after_feedback = [fidelity(feedback(ρ), ψ0) for ρ in states]
    return fids_before_feedback, fids_after_feedback
end