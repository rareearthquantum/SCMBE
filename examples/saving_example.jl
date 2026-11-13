using Plots, JLD2

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

#then we can save the actual data itself in a jld2 file
jldsave("../data/examples/jld2/saving_example.jld2"; rabi_frequency=output_field)
output_field = nothing #remove all previous data just for demonstrative purposes

#then load it
rabi_frequency = load("../data/examples/jld2/saving_example.jld2", "rabi_frequency")

#make plot with the loaded data
fig = heatmap(t_vec,z_vec,transpose(abs2.(rabi_frequency)))

#saving fig, make sure the path exists by making it and then save the fig
mkpath("../data/examples/")
savefig(fig, "../data/examples/plots/saving_example.png") #regular png
savefig(fig, "../data/examples/plots/saving_example.pdf") #can also do pdf for putting in latex since its vector graphics