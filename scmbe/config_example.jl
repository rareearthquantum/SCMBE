
#this just outlines all of the required variables

struct Config

    d_vec::LinRange{Float64}
    t_vec::LinRange{Float64}
    z_vec::LinRange{Float64}
    y_vec::LinRange{Float64}

    input_rabi_field::Array{Complex64,2}

    alpha::Float64
    beta::Float64

end
