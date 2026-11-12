using Plots

include("../scmbe/main_include.jl")

#make grids
d_vec = LinRange(-50, 50, 100)
t_vec = LinRange(0, 1, 512)

#make input field
W(w) = 1/(2*w^2) #width factor within the exponential
N(w) = 1/(w*sqrt(2*pi)) #normalises gaussian

t01 = t_vec[2*end÷10]
area1 = pi/2
w1 = 0.1
first_pulse = @. area1*N(w1)*exp(-W(w1)*(t_vec-t01)^2)

t02 = t_vec[5*end÷10]
area2 = pi
w2 = 0.01
second_pulse = @. area2*N(w2)*exp(-W(w2)*(t_vec-t02)^2)

input_field = first_pulse + second_pulse


#run sim
polarisation, atomic_state = run_0d_simulation(d_vec, t_vec, input_field) #outputs atomic state because the 0d run is to just look at if it works or not
polarisation_intensity = abs2.(polarisation) #not really the intensity but we need to plot real numbers

#display plot
fig = plot(t_vec, polarisation_intensity)
display(fig)
#makes the plot stay up, without this it appears for a couple milliseconds and then vanishes
println("Press Enter to exit...")
readline()