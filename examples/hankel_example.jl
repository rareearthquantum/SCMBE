using Hankel



R = 100.0
Nr = 128
q = QDHT(R,Nr)

r_vec = q.r
pulse = @. exp(-r_vec^2/2)

@show integral = integrateR(abs2.(pulse),q)

hankeled_pulse = q*pulse
unhankeled_pulse = q\hankeled_pulse

@show pulse ≈ unhankeled_pulse