mutable struct AB2Cache{T}
    fprev::T
    fcurr::T
end

AB2Cache(u) = AB2Cache(similar(u), similar(u))

function field_ab2_step!(out_field, in_field, dz, alpha, polarisation_spectrum, rotfactor, cache::AB2Cache)
    @. cache.fcurr = im * alpha * polarisation_spectrum * rotfactor
    @. out_field = in_field + (dz/2) * (3 * cache.fcurr - cache.fprev)
    copyto!(cache.fprev, cache.fcurr)

    return nothing
end

function field_euler_step!(out_field, in_field, dz, alpha, polarisation_spectrum, rotfactor, cache::AB2Cache)
    @. cache.fprev = im * alpha * polarisation_spectrum * rotfactor
    @. out_field = in_field + dz * cache.fprev

    return nothing
end