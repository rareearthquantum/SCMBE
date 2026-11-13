using Plots

include("../scmbe/main_include.jl")

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0,1,128)
z_vec = LinRange(0,1,64)
y_vec = LinRange(-1,1,64)

#make input field
t0 = t_vec[end÷2]
A = 0.1
w_t = 0.1
w_y = 0.1
input_field = A.*exp.(-(1/w_t).*(t_vec.-t0).^2) * exp.(-(1/w_y).*(transpose(y_vec)).^2)

#params
alpha = 0
beta = 0.05

#run sim
output_field = run_2d_simulation(d_vec,t_vec,z_vec,y_vec,input_field,alpha,beta)
#cant plot the complex numb, so plot intensity
output_intensity = abs2.(output_field[end÷2,:,:])

#make and display plot
fig = heatmap(y_vec,z_vec,output_intensity)
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()