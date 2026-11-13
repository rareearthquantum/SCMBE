
function run_0d_simulation(d_vec, t_vec, input_field)

    rabi_frequency = zeros(ComplexF64, length(t_vec))
    rabi_frequency .= input_field

    atomic_state = zeros(ComplexF64, length(d_vec), 2, length(t_vec))
    polarisation = zeros(ComplexF64, length(t_vec))

    unrotate_sigma_grid = [cis(-detuning * t) for detuning in d_vec, t in t_vec]
    rotate_sigma = make_rotate_sigma(d_vec, t_vec)

    atom_cache = @views RK4Cache(atomic_state[:, :, 1])

    #-------------------------------------
    @views compute_polarisation_column!(polarisation, atomic_state, rabi_frequency, unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    return polarisation, atomic_state
end