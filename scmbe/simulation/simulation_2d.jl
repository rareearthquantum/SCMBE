using FFTW

function run_2d_simulation(d_vec, t_vec, z_vec, y_vec, input_field, alpha, beta)

    dz = step(z_vec)

    rabi_frequency = zeros(ComplexF64, length(t_vec), length(z_vec), length(y_vec))
    rabi_frequency[:,1,:] .= input_field
    rabi_frequency_spectrum = similar(rabi_frequency[:, 1, :])

    atomic_state = zeros(ComplexF64, length(d_vec), 2, length(t_vec))

    polarisation = zeros(ComplexF64, length(t_vec), length(y_vec))
    polarisation_spectrum = similar(polarisation)

    unrotate_sigma_grid = [cis(-detuning * t) for detuning in d_vec, t in t_vec]
    rotate_sigma = make_rotate_sigma(d_vec,t_vec)

    ky_grid = 2pi.*fftfreq(length(y_vec), 1/step(y_vec))
    diffrac_rotate = @. cis(beta*ky_grid^2*dz)
    inverse_diffrac_rotate = conj.(diffrac_rotate)

    field_caches = [@views AB2Cache(rabi_frequency[:, 1, l]) for l in eachindex(y_vec)]
    atom_cache = @views RK4Cache(atomic_state[:, :, 1])

    #-------------------------------------------------
    for l in eachindex(y_vec)
        @views compute_polarisation_column!(polarisation[:,l], atomic_state, rabi_frequency[:, 1, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)
    end

    @inbounds for j in 2:length(z_vec)

        polarisation_spectrum .= fft(polarisation, 2)
        rabi_frequency_spectrum .= fft(rabi_frequency[:, j-1, :], 2)
        for l in eachindex(ky_grid)
            @views ab2_step!(
                field_2d!,
                rabi_frequency[:, j, l],
                rabi_frequency_spectrum[:, l],
                z_vec[j],
                dz,
                (alpha, polarisation_spectrum[:, l], diffrac_rotate[l]),
                field_caches[l]
            )
        end
        for i in eachindex(t_vec)
            rabi_frequency[i, j, :] .*= inverse_diffrac_rotate
        end
        @views ifft!(rabi_frequency[:, j, :], 2)


        for l in eachindex(y_vec)
            @views compute_polarisation_column!(polarisation[:,l], atomic_state, rabi_frequency[:, j, l], unrotate_sigma_grid, t_vec, rotate_sigma, atom_cache)
        end

    end

    return rabi_frequency
end
