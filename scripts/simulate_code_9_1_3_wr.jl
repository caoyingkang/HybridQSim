using QuantumToolbox
using CairoMakie
using HybridQSim

code = SurfaceCode_9_1_3_WR()
# ψ0 = 1 / √(2) * (basis(2, 0) + basis(2, 1))
ψ0 = basis(2, 0)
tlist = LinRange(0.0, 100.0, 21)
T1, T2 = 100.0, 80.0
ϵ = 0.1 # coherent error strength
λ = 10.0 # penalty strength
seed = 42

stabs = stabilizer_generators(code)
Mx, Mz = stabilizer_matrix(code)
Uenc = encoding_isometry(code)
Hpen = penalty_hamiltonian(code)

c_ops = get_lindblad_ops(T1, T2, Val(code.num_data_qubits))
Herr = rand_coherent_error(ϵ, Val(code.num_data_qubits); seed=seed)
lut = get_decoding_lut(Mx, Mz)
feedback = FeedbackCPTPMap(stabs, lut, Val(code.num_data_qubits))

Htot = λ * Hpen + Herr
ψ0_enc = Uenc * ψ0

fids_before_feedback, fids_after_feedback = run_simulation(
    H=Htot, ψ0=ψ0_enc, tlist=tlist, feedback=feedback, c_ops=c_ops
)

# For comparison, also run the unencoded simulation
c_ops_unenc = get_lindblad_ops(T1, T2, Val(code.num_logi_qubits))
Herr_unenc = rand_coherent_error(ϵ, Val(code.num_logi_qubits); seed=seed)
Htot_unenc = Herr_unenc
sol_unenc = mesolve(Htot_unenc, ψ0, tlist, c_ops_unenc)
fids_unenc = [fidelity(ρ, ψ0) for ρ in sol_unenc.states]


# Plotting
f = Figure()
ax = Axis(f[1, 1], xlabel="Analog Evolution Time", ylabel="Fidelity", title="[[9,1,3]] Surface Code (WR)")
lines!(ax, tlist, fids_unenc, label="Unencoded", color=:red, linestyle=:solid)
scatter!(ax, tlist, fids_unenc, color=:red)
lines!(ax, tlist, fids_before_feedback, label="Encoded w/o feedback", color=:blue, linestyle=:dash)
scatter!(ax, tlist, fids_before_feedback, color=:blue)
lines!(ax, tlist, fids_after_feedback, label="Encoded w/ feedback", color=:blue, linestyle=:solid)
scatter!(ax, tlist, fids_after_feedback, color=:blue)
axislegend(ax, position=:lb)
figs_dir = joinpath(@__DIR__, "figs")
mkpath(figs_dir)
save(joinpath(figs_dir, "surface_code_9_1_3_wr.png"), f)
display(f)
