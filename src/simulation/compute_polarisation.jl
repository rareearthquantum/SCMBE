
function compute_polarisation_column!(polarisation_column, atomic_state, rabi_frequency_column, unrotate_sigma_grid, t_vec, rotate_sigma, cache)
    fill!(polarisation_column, 0)

    atomic_state[:, 1, 1] .= 0.0
    atomic_state[:, 2, 1] .= -1.0
    @views rk4_custom!(atomic_state, t_vec, (rabi_frequency_column, rotate_sigma), cache)

    @inbounds for j in axes(atomic_state, 3)
        value = zero(eltype(polarisation_column))
        for d in axes(atomic_state, 1)
            value += atomic_state[d, 1, j] * unrotate_sigma_grid[d, j]
        end
        polarisation_column[j] = value / size(atomic_state, 1)
    end

    return nothing
end


function make_rotate_sigma(d_vec, t_vec)
    rotate_sigma = Array{ComplexF64}(undef, length(d_vec), 2*length(t_vec))
    halfdt = step(t_vec)/2
    for j in eachindex(t_vec)
        @. rotate_sigma[:, 2j-1] = cis(d_vec*t_vec[j])
        @. rotate_sigma[:, 2j] = cis(d_vec*(t_vec[j]+halfdt))
    end
    return rotate_sigma
end