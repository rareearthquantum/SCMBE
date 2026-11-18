
function get_alpha_beta(omega_L, d, n)

    const c = 3e8
    const hbar = 6.13e-16 #?
    const mu_0 = 3424 #add

    k = omega_L/c

    alpha = d^2*n*mu_0*omega_L*c/hbar
    beta = 1/(2*k)


    return alpha, beta
end


function convert_to_dimensionless!(T,L,Y,t_vec,z_vec,y_vec,d_vec,input_field,alpha,beta)
    #if you want to not use all arguments, just chuck in some dummy variable like 1 or 0 in those ones
    
    t_vec ./= T
    z_vec ./= L
    y_vec ./= Y
    d_vec .*= T

    input_field .*= T

    alpha .*= (L*T)
    beta .*= (L/Y^2)

    return nothing
end

function convert_back_to_dimensioned!(T,L,Y,t_vec,z_vec,y_vec,d_vec,input_field,alpha,beta)
    #if you want to not use all arguments, just chuck in some dummy variable like 1 or 0 in those ones
    
    t_vec .*= T
    z_vec .*= L
    y_vec .*= Y
    d_vec ./= T

    input_field ./= T

    alpha ./= (L*T)
    beta ./= (L/Y^2)

    return nothing
end