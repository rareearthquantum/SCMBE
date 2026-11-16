
function run_cylindrical_simulation(d_vec, t_vec, z_vec, q, input_field, alpha, beta)

    dz = step(z_vec)

    rabi_frequency = zeros(ComplexF64, length(t_vec), length(z_vec), length(q.r))
    rabi_frequency[:, 1, :] .= input_field
    rabi_frequency_spectrum = similar(rabi_frequency[:, 1, :])

    atomic_state = zeros(ComplexF64, length(d_vec), 2, length(t_vec))

    polarisation = zeros(ComplexF64, length(t_vec), length(q.r))
    polarisation_spectrum = similar(polarisation)

    unrotate_sigma_grid = cis.(-d_vec * transpose(t_vec))
    rotate_sigma = make_rotate_sigma(d_vec, t_vec)

    diffrac_rotate = cis.(z_vec * beta.*transpose(q.k).^2)
    inverse_diffrac_rotate = conj.(diffrac_rotate)

    field_caches = [@views AB2Cache(rabi_frequency[:, 1, l]) for l in eachindex(q.r)]
    atom_cache = @views RK4Cache(atomic_state[:, :, 1])

    #-------------------------------------------------
    #field at z[1] complete so we do polarisation too to make it even
    for l in eachindex(q.r)
        @views compute_polarisation_column!(polarisation[:, l], atomic_state, rabi_frequency[:, 1, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)
    end

    #first step of ab2 has to be an euler step, so might as well do it outside
    compute_radial_field_spectrum!(field_euler_step!, q, polarisation_spectrum, polarisation, rabi_frequency_spectrum, rabi_frequency, 2,
        alpha, diffrac_rotate, inverse_diffrac_rotate, field_caches, dz)

    #evening things out again
    for l in eachindex(q.r)
        @views compute_polarisation_column!(polarisation[:, l], atomic_state, rabi_frequency[:, 2, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)
    end

    #has to start at 3
    @inbounds for j in 3:length(z_vec)

        compute_radial_field_spectrum!(field_ab2_step!, q, polarisation_spectrum, polarisation, rabi_frequency_spectrum, rabi_frequency, j,
            alpha, diffrac_rotate, inverse_diffrac_rotate, field_caches, dz)

        for l in eachindex(q.r)
            @views compute_polarisation_column!(polarisation[:, l], atomic_state, rabi_frequency[:, j, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)
        end

    end

    return rabi_frequency
end
