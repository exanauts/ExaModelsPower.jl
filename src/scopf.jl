using GOC3Benchmark, JSON

function scopf_model(
    filename, uc_filename;
    backend = nothing,
    T = Float64,
    kwargs...
        )


    uc_data = JSON.parsefile(uc_filename)
    data = GOC3Benchmark.get_data_from_file(filename)
    sc_data, lengths = parse_sc_data(data, uc_data)
    return convert_data(sc_data, backend)
end

