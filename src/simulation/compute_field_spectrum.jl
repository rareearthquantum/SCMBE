using FFTW

function compute_field_spectrum!(f!, polarisation_spectrum, polarisation, rabi_frequency_spectrum, rabi_frequency, z_index,
    alpha, diffrac_rotate, inverse_diffrac_rotate, field_caches, dz)
    polarisation_spectrum .= fft(polarisation, 2)
    rabi_frequency_spectrum .= fft(rabi_frequency[:, z_index-1, :], 2)
    for l in axes(rabi_frequency, 3)
        @views f!(
            rabi_frequency[:, z_index, l],
            rabi_frequency_spectrum[:, l],
            dz,
            alpha, polarisation_spectrum[:, l], diffrac_rotate[l],
            field_caches[l]
        )
    end
    for i in axes(rabi_frequency, 1)
        rabi_frequency[i, z_index, :] .*= inverse_diffrac_rotate
    end
    @views ifft!(rabi_frequency[:, z_index, :], 2)

    return nothing
end