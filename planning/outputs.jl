#well the simulations have to output something that can be processed later
#these functions should ideally be used within simulation.jl as some variables might need to be saved as they go

#IMPORTANT, so what will be discarded in the simulations? nothing? and what about the overhead of opening and closing stuff?

#ANSWER, starting point will be just to store the field in either a jld2 or hs5 (or whatever its called)