using QuantumToolbox
using CairoMakie
using HybridQSim

code = Code_5_1_3()

stabs = stabilizer_generators(code)
Mx, Mz = stabilizer_matrix(code)
Uenc = encoding_isometry(code)
Hpen = penalty_hamiltonian(code)

lut = get_decoding_lut(Mx, Mz)
feedback = FeedbackCPTPMap(stabs, lut, Val(code.num_data_qubits))

H = Hpen # TODO: add coherent error Herr
ψ0_unenc = 1/√(2) * (basis(2, 0) + basis(2, 1))
ψ0_enc = Uenc * ψ0_unenc
tlist = LinRange(0.0, 100.0, 21)
T1, T2 = 100.0, 80.0

fids_before_feedback, fids_after_feedback = run_simulation(
    H=H, ψ0=ψ0_enc, tlist=tlist, feedback=feedback,
    c_ops=get_lindblad_ops(T1, T2, Val(code.num_data_qubits))
)

# For comparison, also run the unencoded simulation
H_unenc = zero(sigmaz())
sol_unenc = mesolve(
    H_unenc, ψ0_unenc, tlist, get_lindblad_ops(T1, T2, Val(code.num_logi_qubits))
)
fids_unenc = [fidelity(ρ, ψ0_unenc) for ρ in sol_unenc.states]


# Plotting
f = Figure()
ax = Axis(f[1,1], xlabel="Analog Evolution Time", ylabel="Fidelity", title="[[5,1,3]] Code")
lines!(ax, tlist, fids_unenc, label="Unencoded Analog Evolve", color=:red, linestyle=:solid)
scatter!(ax, tlist, fids_unenc, color=:red)
lines!(ax, tlist, fids_before_feedback, label="Encoded Analog Evolve", color=:blue, linestyle=:dash)
scatter!(ax, tlist, fids_before_feedback, color=:blue)
lines!(ax, tlist, fids_after_feedback, label="Encoded Analog Evolve + Digital Feedback", color=:blue, linestyle=:solid)
scatter!(ax, tlist, fids_after_feedback, color=:blue)
axislegend(ax, position=:lb)
save(joinpath(@__DIR__, "fig_5_1_3.png"), f)
display(f)
