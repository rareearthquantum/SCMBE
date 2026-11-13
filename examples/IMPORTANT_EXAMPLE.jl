# THIS EXAMPLE IS IMPORTANT

#-PACKAGE IMPORT------------------------------------------------------------------------------------------------------------------------------------------------

#=
To have access to the core functions, import or use SCMBE (this project)
Use or import anything else you require
=#
using SCMBE, Plots

#-SOLVER SETUP------------------------------------------------------------------------------------------------------------------------------------------------

#=
You have to define the grid yourself and then pass it in to the solver
This depends on which solver you are using, check carefully
Do it however you wish, with magic numbers, with named variables, your choice
=#
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0,1,128)
z_vec = LinRange(0,1,64)
y_vec = LinRange(-1,1,64)

# Define the input field to pass in, again check solver inputs
t0 = t_vec[end÷2]
A = 0.1
w_t = 0.1
w_y = 0.1
input_field = A.*exp.(-(1/w_t).*(t_vec.-t0).^2) * exp.(-(1/w_y).*(transpose(y_vec)).^2)

# Set the alpha and beta params, again depends on scenario
alpha = 1
beta = 0.05

#-SOLVE------------------------------------------------------------------------------------------------------------------------------------------------

#=
Use one of the "run simulation" functions to get the output field
Note that in the 0d case, the polarisation and atomic_states will be outputted instead
=#
output_field = run_2d_simulation(d_vec,t_vec,z_vec,y_vec,input_field,alpha,beta)

#-PROCESS OUTPUT------------------------------------------------------------------------------------------------------------------------------------------------

# Do whatever you wish with the data, plot it, save it, process it, etc.

# The intensity is a useful quantity, completely real and positive and corresponds closely to energy (this comment needs refinement/revision)
output_intensity = abs2.(output_field[end÷2,:,:])

# Plotting for demonstration
fig = heatmap(y_vec,z_vec,output_intensity)
display(fig)
println("Press Enter to exit...")
readline()