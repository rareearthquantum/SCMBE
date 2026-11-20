function echo_efficiency(t_vec, z_vec, t01, t02, w1, field)
    tau = t02 - t01
    echo_time = t01 + 2tau
    input_peaks = zeros(length(z_vec))
    echo_peaks = zeros(length(z_vec))
    input_index = searchsortedfirst(t_vec, t01)
    echo_index = searchsortedfirst(t_vec, echo_time)
    width_index = searchsortedfirst(t_vec, w1)
    inputrange = clamp.((input_index-3width_index):(input_index+3width_index), 1, length(t_vec))
    for i in 1:length(z_vec)
        input_peaks[i] = sum(field[inputrange, i])
        echo_peaks[i] = sum(field[(echo_index-2width_index):end, i])
    end
    efficiency = echo_peaks ./ input_peaks[1]
    return efficiency
end

