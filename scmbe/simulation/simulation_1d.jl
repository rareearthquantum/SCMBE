function run_1d_simulation(d_vec, t_vec, z_vec, input_field, alpha)

    dz = step(z_vec)

    rabi_frequency = zeros(ComplexF64, length(t_vec), length(z_vec))
    rabi_frequency[:, 1] .= input_field

    atomic_state = zeros(ComplexF64, length(d_vec), 2, length(t_vec))

    polarisation = zeros(ComplexF64, length(t_vec))

    unrotate_sigma_grid = [cis(-detuning * t) for detuning in d_vec, t in t_vec]
    field_cache = @views AB2Cache(rabi_frequency[:, 1])

    rotate_sigma = make_rotate_sigma(d_vec, t_vec)

    atom_cache = @views RK4Cache(atomic_state[:, :, 1])

    #-------------------------------------
    @views compute_polarisation_column!(polarisation[:], atomic_state, rabi_frequency[:, 1], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    @inbounds for j in 2:length(z_vec)

        @views ab2_step!(field_2d!, rabi_frequency[:, j], rabi_frequency[:, j-1], z_vec, dz, (alpha, polarisation, 1.0), field_cache)

        @views compute_polarisation_column!(polarisation, atomic_state, rabi_frequency[:, j], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    end

    return rabi_frequency
end