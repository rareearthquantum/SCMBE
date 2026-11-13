using SCMBE, Plots, JLD2

# Simulation code is identical to weak_pulse_example.jl

d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0, 1, 128)
z_vec = LinRange(0, 1, 64)

t0 = t_vec[end÷2]
A = 0.1
w = 0.1
input_field = @. A*exp(-(1/w)*(t_vec-t0)^2)

alpha = 50

output_field = run_1d_simulation(d_vec, t_vec, z_vec, input_field, alpha)

#=
Make sure when saving that you know exactly where you are saving to and from

Use @__DIR__ to get the directory that the file running the command is in, for this example it would be "DriveLetter:/blah/blahblah/SCMBE/examples/"
Use pwd() to get the path julia was launched from, jldsave and savefig use that internally
Use pkgdir(SCMBE) for the root directory of the package (careful)

@__DIR__ is the most useful at the moment
=#
plot_dir = joinpath(@__DIR__, "plots")
data_dir = joinpath(@__DIR__, "data")


# Then we can save the actual data itself in a jld2 file
mkpath(data_dir) #make sure path exists by creating it
jldsave(joinpath(data_dir, "saving_example.jld2"); rabi_frequency=output_field)
output_field = nothing #remove all previous data just for demonstrative purposes

# Then load it
rabi_frequency = load(joinpath(data_dir, "saving_example.jld2"), "rabi_frequency")

# Make plot with the loaded data
ENV["GKSwstype"] = "100" #turns off the plot popup, useful for cluster
fig = heatmap(t_vec, z_vec, transpose(abs2.(rabi_frequency)))

# Saving fig, make sure the path exists by making it and then save the fig
mkpath(plot_dir)
savefig(fig, joinpath(plot_dir, "saving_example.png")) #regular png
savefig(fig, joinpath(plot_dir, "saving_example.pdf")) #can also do pdf for putting in latex since its vector graphics


#just cleanup, normally you don't need this since you want to keep your files
println("Please view your saved files before they are removed.\nPress Enter to exit and remove your saved files...")
readline()
rm(data_dir, force=true, recursive=true)
rm(plot_dir, force=true, recursive=true)