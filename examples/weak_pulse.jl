using Plots

include("../scmbe/main_include.jl")

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0,1,128)
z_vec = LinRange(0,1,64)

#make input field
t0 = t_vec[end÷2]
A = 0.1
w = 0.1
input_field = @. A*exp(-(1/w)*(t_vec-t0)^2)

#param
alpha = 50

#run sim
output_field = run_1d_simulation(d_vec,t_vec,z_vec,input_field,alpha)
#cant plot the complex numb, so plot intensity
output_intensity = abs2.(output_field)

#display plot
fig = heatmap(t_vec,z_vec,transpose(output_intensity))
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()