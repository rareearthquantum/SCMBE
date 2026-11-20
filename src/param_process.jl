
function get_alpha_beta(omega_L, d, n)

    c = 299792458.0
    hbar = 1.054571817e-34
    mu_0 = 1.25663706127e-6

    k = omega_L/c

    alpha = d^2*n*mu_0*omega_L*c/hbar
    beta = 1/(2*k)


    return alpha, beta
end


function convert_to_dimensionless(T,L,Y,t_vec,z_vec,y_vec,d_vec,input_field,alpha,beta,ky_vec,q)
    #if you want to not use all arguments, just chuck in some dummy variable like 1 or 0 in those ones
    
    new_t_vec = t_vec ./ T
    new_z_vec = z_vec ./ L
    new_y_vec = y_vec ./ Y
    new_d_vec = d_vec .* T

    new_input_field = input_field .* T

    new_alpha = alpha * (L*T)
    new_beta = beta * (L/Y^2)

    new_ky_vec = ky_vec .* Y

    new_q = QDHT(q.R / Y, length(q.r))

    return (new_t_vec,new_z_vec,new_y_vec,new_d_vec,new_input_field,new_alpha,new_beta,new_ky_vec,new_q)
end