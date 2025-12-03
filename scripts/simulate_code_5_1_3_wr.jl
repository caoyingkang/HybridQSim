## simulation task setup
using QuantumToolbox
using CairoMakie
using Random
using Statistics
using HybridQSim

code = Code_5_1_3_WR()
tlist = LinRange(0.0, 100.0, 21)
T1, T2 = 100.0, 80.0
ϵ = 0.1 # coherent error strength
λ = 10.0 # penalty strength
seed = 42
num_samples = 20

rng = Random.default_rng()
Random.seed!(rng, seed)

KET0 = basis(2, 0)
KET1 = basis(2, 1)
KETPLUS = 1 / √(2) * (KET0 + KET1)
KETPLUSI = 1 / √(2) * (KET0 + 1im * KET1)

stabs = stabilizer_generators(code)
Mx, Mz = stabilizer_matrix(code)
Uenc = encoding_isometry(code)
Hpen = penalty_hamiltonian(code)

c_ops = get_lindblad_ops(T1, T2, Val(code.num_data_qubits))
c_ops_unenc = get_lindblad_ops(T1, T2, Val(code.num_logi_qubits))
lut = get_decoding_lut(Mx, Mz)
feedback = FeedbackCPTPMap(stabs, lut, Val(code.num_data_qubits))

ψz = foldl(⊗, ntuple(i -> KET0, Val(code.num_logi_qubits))) # |0⟩⊗|0⟩⊗...⊗|0⟩
ψx = foldl(⊗, ntuple(i -> KETPLUS, Val(code.num_logi_qubits))) # |+⟩⊗|+⟩⊗...⊗|+⟩
ψy = foldl(⊗, ntuple(i -> KETPLUSI, Val(code.num_logi_qubits))) # |+i⟩⊗|+i⟩⊗...⊗|+i⟩
ψz_enc = Uenc * ψz
ψx_enc = Uenc * ψx
ψy_enc = Uenc * ψy

## ------------------------ Encoded Evolution ------------------------
all_fids_before_feedback_z = zeros(length(tlist), num_samples)
all_fids_after_feedback_z = zeros(length(tlist), num_samples)
all_fids_before_feedback_x = zeros(length(tlist), num_samples)
all_fids_after_feedback_x = zeros(length(tlist), num_samples)
all_fids_before_feedback_y = zeros(length(tlist), num_samples)
all_fids_after_feedback_y = zeros(length(tlist), num_samples)

for i in 1:num_samples
    Herr = rand_coherent_error(ϵ, Val(code.num_data_qubits); rng=rng)
    Htot = λ * Hpen + Herr

    fids_before_feedback_z, fids_after_feedback_z = run_simulation(
        H=Htot, ψ0=ψz_enc, tlist=tlist, feedback=feedback, c_ops=c_ops
    )
    all_fids_before_feedback_z[:, i] = fids_before_feedback_z
    all_fids_after_feedback_z[:, i] = fids_after_feedback_z

    fids_before_feedback_x, fids_after_feedback_x = run_simulation(
        H=Htot, ψ0=ψx_enc, tlist=tlist, feedback=feedback, c_ops=c_ops
    )
    all_fids_before_feedback_x[:, i] = fids_before_feedback_x
    all_fids_after_feedback_x[:, i] = fids_after_feedback_x

    fids_before_feedback_y, fids_after_feedback_y = run_simulation(
        H=Htot, ψ0=ψy_enc, tlist=tlist, feedback=feedback, c_ops=c_ops
    )
    all_fids_before_feedback_y[:, i] = fids_before_feedback_y
    all_fids_after_feedback_y[:, i] = fids_after_feedback_y
end


## ------------------------ Unencoded Evolution ------------------------
all_fids_unencoded_z = zeros(length(tlist), num_samples)
all_fids_unencoded_x = zeros(length(tlist), num_samples)
all_fids_unencoded_y = zeros(length(tlist), num_samples)

for i in 1:num_samples
    Herr_unenc = rand_coherent_error(ϵ, Val(code.num_logi_qubits); rng=rng)
    Htot_unenc = Herr_unenc

    sol_z = mesolve(Htot_unenc, ψz, tlist, c_ops_unenc)
    fids_unencoded_z = [fidelity(ρ, ψz) for ρ in sol_z.states]
    all_fids_unencoded_z[:, i] = fids_unencoded_z

    sol_x = mesolve(Htot_unenc, ψx, tlist, c_ops_unenc)
    fids_unencoded_x = [fidelity(ρ, ψx) for ρ in sol_x.states]
    all_fids_unencoded_x[:, i] = fids_unencoded_x

    sol_y = mesolve(Htot_unenc, ψy, tlist, c_ops_unenc)
    fids_unencoded_y = [fidelity(ρ, ψy) for ρ in sol_y.states]
    all_fids_unencoded_y[:, i] = fids_unencoded_y
end


## ------------------------ Plotting ------------------------
fig = Figure(size=(800, 600))

ax = Axis(fig[1, 1:2], ylabel="Fidelity", title="[[5,1,3]] Code (WR), Initial State |0⟩")
μ, σ = mean(all_fids_before_feedback_z, dims=2)[:], std(all_fids_before_feedback_z, dims=2)[:]
l1 = lines!(ax, tlist, μ, label="Encoded w/o feedback", color=:red)
scatter!(ax, tlist, μ, color=:red)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:red, 0.3))
μ, σ = mean(all_fids_after_feedback_z, dims=2)[:], std(all_fids_after_feedback_z, dims=2)[:]
l2 = lines!(ax, tlist, μ, label="Encoded w/ feedback", color=:blue)
scatter!(ax, tlist, μ, color=:blue)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:blue, 0.3))
μ, σ = mean(all_fids_unencoded_z, dims=2)[:], std(all_fids_unencoded_z, dims=2)[:]
l3 = lines!(ax, tlist, μ, label="Unencoded", color=:black)
scatter!(ax, tlist, μ, color=:black)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:black, 0.3))

ax = Axis(fig[2, 1:2], ylabel="Fidelity", title="[[5,1,3]] Code (WR), Initial State |+⟩")
μ, σ = mean(all_fids_before_feedback_x, dims=2)[:], std(all_fids_before_feedback_x, dims=2)[:]
lines!(ax, tlist, μ, label="Encoded w/o feedback", color=:red)
scatter!(ax, tlist, μ, color=:red)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:red, 0.3))
μ, σ = mean(all_fids_after_feedback_x, dims=2)[:], std(all_fids_after_feedback_x, dims=2)[:]
lines!(ax, tlist, μ, label="Encoded w/ feedback", color=:blue)
scatter!(ax, tlist, μ, color=:blue)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:blue, 0.3))
μ, σ = mean(all_fids_unencoded_x, dims=2)[:], std(all_fids_unencoded_x, dims=2)[:]
lines!(ax, tlist, μ, label="Unencoded", color=:black)
scatter!(ax, tlist, μ, color=:black)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:black, 0.3))

ax = Axis(fig[3, 1:2], xlabel="Analog Evolution Time", ylabel="Fidelity", title="[[5,1,3]] Code (WR), Initial State |+i⟩")
μ, σ = mean(all_fids_before_feedback_y, dims=2)[:], std(all_fids_before_feedback_y, dims=2)[:]
lines!(ax, tlist, μ, label="Encoded w/o feedback", color=:red)
scatter!(ax, tlist, μ, color=:red)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:red, 0.3))
μ, σ = mean(all_fids_after_feedback_y, dims=2)[:], std(all_fids_after_feedback_y, dims=2)[:]
lines!(ax, tlist, μ, label="Encoded w/ feedback", color=:blue)
scatter!(ax, tlist, μ, color=:blue)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:blue, 0.3))
μ, σ = mean(all_fids_unencoded_y, dims=2)[:], std(all_fids_unencoded_y, dims=2)[:]
lines!(ax, tlist, μ, label="Unencoded", color=:black)
scatter!(ax, tlist, μ, color=:black)
band!(ax, tlist, μ .- σ, μ .+ σ, color=(:black, 0.3))

Legend(fig[1, 3], [l1, l2, l3], ["Encoded w/o feedback", "Encoded w/ feedback", "Unencoded"])

display(fig)
figs_dir = joinpath(@__DIR__, "figs")
mkpath(figs_dir)
save(joinpath(figs_dir, "code_5_1_3_wr.png"), fig)
