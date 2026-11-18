using SCMBE, Plots

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0, 1, 128)
z_vec = LinRange(0, 1, 64)

#make input field
t0 = t_vec[end÷2]
A = pi/16
w = z_vec[end]/5
input_field = @. A*(1/(w*sqrt(2pi)))*exp(-(1/(2*w^2))*(t_vec-t0)^2)

#param
alpha = 100

#run sim
output_field = run_1d_simulation(d_vec, t_vec, z_vec, input_field, alpha)
output_intensity = abs2.(output_field) #cant plot complex numbers directly, so plot intensity

#display plot
fig1 = heatmap(t_vec, z_vec, transpose(output_intensity); c=:viridis)
fig2 = plot(z_vec, vec(sum(output_intensity, dims=1)))
fig = plot(fig1, fig2; layout=(1, 2), size=(1600, 600))
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()