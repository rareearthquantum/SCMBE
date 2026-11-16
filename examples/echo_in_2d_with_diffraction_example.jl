using SCMBE, Plots

#make grids
t_vec = LinRange(0, 1, 1024)
z_vec = LinRange(0, 1, 128)
y_vec = LinRange(-1, 1, 64)
#careful when making detunings
d_vec = LinRange(-100, 100, 256)


#make input field
W(w) = 1/(2*w^2) #width factor within the exponential
N(w) = 1/(w*sqrt(2*pi)) #normalises gaussian

w_y = 1.0
transverse_profile_unit_max_amplitude = @. exp(-W(w_y)*y_vec^2)

t01 = t_vec[1*end÷10]
area1 = pi/2
w1 = t_vec[end]/50
first_pulse = @. area1 * N(w1) * exp(-W(w1)*(t_vec-t01)^2)
first_pulse = first_pulse * transpose(transverse_profile_unit_max_amplitude)

t02 = t_vec[4*end÷10]
area2 = pi
w2 = w1/4
second_pulse = @. area2 * N(w2) * exp(-W(w2)*(t_vec-t02)^2)
second_pulse = second_pulse * transpose(transverse_profile_unit_max_amplitude)

input_field = first_pulse + second_pulse

#params
alpha = 500
beta = 0

#run sim
output_field = @time run_2d_simulation(d_vec, t_vec, z_vec, y_vec, input_field, alpha, beta)
output_intensity = abs2.(output_field)
intensity_summed_over_y = dropdims(sum(output_intensity; dims=3); dims=3)

#make and display plots
fig1 = heatmap(t_vec, z_vec, transpose(intensity_summed_over_y); c=:viridis)
fig2 = plot(); [plot!(t_vec, intensity_summed_over_y[:, j]; legend=false, c=cgrad(:rainbow, rev=true), line_z=j) for j in floor.(Int, LinRange(1, length(z_vec), 10))];
fig = plot(fig1, fig2; layout=(1, 2), size=(1600, 600))
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()