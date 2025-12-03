# Definitions for abstract code interface.

export stabilizer_generators, stabilizer_matrix, encoding_isometry, penalty_hamiltonian

abstract type AbstractCode end

"""
Return an array of independent stabilizer generators.
"""
function stabilizer_generators(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end


"""
Return the stabilizer matrix in binary symplectic form, separated into X part and Z part.
"""
function stabilizer_matrix(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end


"""
Return the array of logical Z operators, one for each logical qubit.
"""
function logical_Zs(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end


"""
Return the array of logical X operators, one for each logical qubit.
"""
function logical_Xs(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end


"""
Return the encoding isometry.
"""
function encoding_isometry(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end


"""
Return the penalty Hamiltonian, not including the penalty coefficient.
"""
function penalty_hamiltonian(::AbstractCode)
    throw(ErrorException("Not implemented!"))
end