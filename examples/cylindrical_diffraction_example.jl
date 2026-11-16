using SCMBE, Plots, Hankel

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0, 1, 128)
z_vec = LinRange(0, 1, 64)
q = QDHT(2.0, 128)

#make input field
t0 = t_vec[end÷2]
A = 0.1
w_t = 0.1
w_y = 0.1
input_field = A .* exp.(-(1/w_t) .* (t_vec .- t0) .^ 2) * exp.(-(1/w_y) .* (transpose(q.r)) .^ 2)

#params
alpha = 0
beta = 0.05

#run sim
output_field = @time run_cylindrical_simulation(d_vec, t_vec, z_vec, q, input_field, alpha, beta)
output_intensity = abs2.(output_field[end÷2, :, :])

#make and display plot
fig = heatmap(q.r, z_vec, output_intensity; c=:viridis)
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()