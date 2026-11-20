using SCMBE, Plots, Hankel
include("./helpers.jl")


#-User Inputs---------------------------------------
#input physical params
omega_L = 3.11e15
dipole_moment = 3.7e-32
number_density = 2e20
sample_length = 10e-3
first_pulse_width = 1e-6
second_pulse_width = first_pulse_width/4
transverse_pulse_width = 100e-6
characteristic_time = first_pulse_width
characteristic_transverse_length = transverse_pulse_width
detuning_half_width = 2e6

#make grids
d_vec = 2pi .* LinRange(-detuning_half_width, detuning_half_width, 256)
t_vec = LinRange(0, 25e-6, 1024)
z_vec = LinRange(0, sample_length, 64)
q = QDHT(500e-6, 64)

t01 = 5e-6
area1 = pi/2

t02 = 12e-6
area2 = pi

#-Making things from user inputs----------------------------------------------------
#make input field
W(w) = 1/(2*w^2) #width factor within the exponential
N(w) = 1/(w*sqrt(2*pi)) #normalises gaussian
transverse_profile = @. exp(-W(transverse_pulse_width)*q.r^2)
first_pulse = (@. area1 * N(first_pulse_width) * exp(-W(first_pulse_width)*(t_vec-t01)^2)) * transpose(transverse_profile)
second_pulse = (@. area2 * N(second_pulse_width) * exp(-W(second_pulse_width)*(t_vec-t02)^2)) * transpose(transverse_profile)
input_field = first_pulse + second_pulse

#params
alpha, beta = get_alpha_beta(omega_L, dipole_moment, number_density)

#-Make dimensionless and simulate----------------------------------------------
#converting to dimensionless params
new_t, new_z, new_y, new_d, new_field, new_alpha, new_beta, new_ky_vec, new_q =
    convert_to_dimensionless(characteristic_time, sample_length, characteristic_transverse_length, t_vec, z_vec, 1:2, d_vec, input_field, alpha, beta, 1:2, q)

#run sim
output_field = @time run_cylindrical_simulation(new_d, new_t, new_z, new_q, new_field, new_alpha, new_beta)

#re-dimensioning the output field
output_field ./= characteristic_time

#-Processing---------------------------------------------------------------------
intensity_disk_integrated = 2pi .* dropdims(integrateR(abs2.(output_field), q; dim=3); dims=3)

#make and display plots
fig1 = heatmap(t_vec, z_vec, transpose(intensity_disk_integrated); c=cgrad(:viridis, scale=:exp)) #watch the scale
fig2 = plot();
[plot!(t_vec, intensity_disk_integrated[:, j]; legend=false, c=cgrad(:rainbow, rev=true), line_z=j, cb=false) for j in floor.(Int, LinRange(1, length(z_vec), 10))];
fig3 = plot(z_vec, echo_efficiency(t_vec, z_vec, t01, t02, first_pulse_width, intensity_disk_integrated))
fig = plot(fig1, fig2, fig3; layout=(2, 2), size=(1600, 1200))

#saving
plot_dir = joinpath(@__DIR__, "plots")
mkpath(plot_dir)
savefig(fig, joinpath(plot_dir, "dimensioned_example.pdf"))

#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
display(fig)
println("Press Enter to exit...")
readline()

