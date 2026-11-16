using Hankel

function compute_radial_field_spectrum!(f!,q, polarisation_spectrum, polarisation, rabi_frequency_spectrum, rabi_frequency, z_index,
    alpha, diffrac_rotate, inverse_diffrac_rotate, field_caches, dz)

    for i in axes(polarisation_spectrum, 1)
        polarisation_spectrum[i, :] .= q * polarisation[i,:]
        rabi_frequency_spectrum .= q * rabi_frequency[i, z_index-1, :]
    end
    
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
        rabi_frequency[i, z_index, :] .\= q
    end

    return nothing
end