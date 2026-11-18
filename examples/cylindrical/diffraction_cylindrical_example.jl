using SCMBE, Plots, Hankel

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0, 1, 128)
z_vec = LinRange(0, 1, 128)
q = QDHT(1.0, 256)

#make input field
W(w) = 1/(2*w^2) #width factor within the exponential
N(w) = 1/(w*sqrt(2*pi)) #normalises gaussian

t0 = t_vec[end÷2]
A = pi/16
w_t = 0.1
w_y = 0.2
input_field = A .* N(w_t) .* exp.(-W(w_t) .* (t_vec .- t0) .^ 2) * exp.(-W(w_y) .* (transpose(q.r)) .^ 2)

#params
alpha = 0
beta = 0.04

#run sim
output_field = @time run_cylindrical_simulation(d_vec, t_vec, z_vec, q, input_field, alpha, beta)
output_intensity = abs2.(output_field[end÷2, :, :])

#make and display plot
fig = heatmap(q.r, z_vec, output_intensity; c=:viridis)
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()