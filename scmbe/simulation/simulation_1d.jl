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
    #field at z[1] complete so we do polarisation too to make it even
    @views compute_polarisation_column!(polarisation[:], atomic_state, rabi_frequency[:, 1], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    #first step of ab2 has to be an euler step, so might as well do it outside
    @views field_euler_step!(rabi_frequency[:, 2], rabi_frequency[:, 1], dz, alpha, polarisation, 1.0, field_cache)

    #evening things out again
    @views compute_polarisation_column!(polarisation, atomic_state, rabi_frequency[:, 2], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    @inbounds for j in 3:length(z_vec)

        @views field_ab2_step!(rabi_frequency[:, j], rabi_frequency[:, j-1], dz, alpha, polarisation, 1.0, field_cache)

        @views compute_polarisation_column!(polarisation, atomic_state, rabi_frequency[:, j], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)

    end

    return rabi_frequency
end