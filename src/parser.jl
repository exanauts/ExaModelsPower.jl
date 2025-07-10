convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(convert_array(d, backend) for d in data)

# FOR CONVENIENCE, the MATPOWER file spec
# https://matpower.app/manual/matpower/DataFileFormat.html

struct BusData{T <: Real}
    bus_i :: Int
    type :: Int
    pd :: T
    qd :: T
    gs :: T
    bs :: T
    area :: T
    Vm :: T
    Va :: T
    baseKV :: T
    zone :: T
    vmax :: T
    vmin :: T
end
struct BranchData{T <: Real}
    fbus :: Int
    tbus :: Int
    r :: T
    x :: T
    b :: T
    ratea ::Int
    rateb :: Int
    ratec :: Int
    ratio :: T
    angle :: T
    status :: Int
    angmin :: T
    angmax :: T
    ratea_sq :: Int
    c1 :: T
    c2 :: T
    c3 :: T
    c4 :: T
    c5 :: T
    c6 :: T
    c7 :: T
    c8 :: T
end
struct StorageData{T <: Real}
    storage_bus :: T
    ps :: Int
    qs :: T
    energy :: T
    energy_rating :: T
    charge_rating :: T
    discharge_rating :: T
    charge_efficiency :: T
    discharge_efficiency :: T
    thermal_rating :: T
    qmin :: T
    qmax :: T
    r :: T
    x :: T
    p_loss :: T
    q_loss :: T
    status :: Int
end
struct GenData{T <: Real}
    bus :: Int
    pg :: T
    qg :: T
    qmax :: T
    qmin :: T
    vg :: T
    mBase :: T
    status :: Int
    pmax :: T
    pmin :: T
    startup :: Union{T, Nothing}
    shutdown :: Union{T, Nothing}
    n :: Union{Int, Nothing}
    c :: Union{NTuple{3, T}, Nothing}
    i :: Int
end

struct Data{T <: Real}
    version :: String
    baseMVA :: Float64
    bus :: Vector{BusData{T}}
    gen :: Vector{GenData{T}}
    branch :: Vector{BranchData{T}}
    storage :: Vector{StorageData{T}}
    ratea :: Vector{Int}
    arc :: Vector{Tuple{Int, Int, Int}}
    ref_buses :: Vector{Int}
end

MATPOWER_VAR_TYPES = [
    (name = "version", type = String);
    (name = "baseMVA", type = Float64);
    (name = "bus", type = BusData);
    (name = "gen", type = GenData);
    (name = "gencost", type = GenData);
    (name = "branch", type = BranchData);
    (name = "storage", type = StorageData);
]

function parse_matpower_file(::Type{T}, fname :: String) :: Data{T} where T <: Real
    lines = split(read(open(fname), String), "\n")
    in_array = false
    cur_key = ""
    version = ""
    baseMVA = 0
    bus :: Vector{BusData{T}} = []
    gen :: Vector{GenData{T}} = []
    branch :: Vector{BranchData{T}} = []
    storage :: Vector{StorageData{T}} = []
    line_ind = 1
    line = lines[line_ind]
    data_patt = r"[^\s=;\[\]]+|[=;\[\]]"
    col_patt = r"[^\s]+"
    type = Missing
    row_num = 1
    comment = ""
    col_inds = Dict()

    while true
        words = [m.match for m in eachmatch(data_patt, line)]

        if in_array && length(line) != 0 && line[1] != '%'
            squares = findall(s -> s == "]", words)
            first_sq = length(squares) == 0 ? typemax(Int64) : squares[1]
            first_semi = length(squares) == 0 ? typemax(Int64) : squares[1]
            items :: Vector{T} = map(s -> parse(T, s), words[1:min(min(first_semi - 1, first_sq - 1), length(words) - 1)])

            if length(items) == 0 && length(words) >= 2 && words[length(words)-1] == "]" && words[length(words)] == ";"
                in_array = false
            elseif length(words) != 0 && last(words) != ";"
                error("Invalid matpower file. Line $(line_ind) array doesn't end with ; or ];")
            elseif length(items) != 0
                if cur_key == "bus"
                    push!(bus, BusData(
                        round(Int, items[col_inds["bus_i"]]),
                        round(Int, items[col_inds["type"]]),
                        items[col_inds["Pd"]],
                        items[col_inds["Qd"]],
                        items[col_inds["Gs"]],
                        items[col_inds["Bs"]],
                        items[col_inds["area"]],
                        items[col_inds["Vm"]],
                        items[col_inds["Va"]],
                        items[col_inds["baseKV"]],
                        items[col_inds["zone"]],
                        items[col_inds["Vmax"]],
                        items[col_inds["Vmin"]],
                    ))
                elseif cur_key == "gen"
                    push!(gen, GenData(
                        round(Int, items[col_inds["bus"]]),
                        items[col_inds["Pg"]],
                        items[col_inds["Qg"]],
                        items[col_inds["Qmax"]],
                        items[col_inds["Qmin"]],
                        items[col_inds["Vg"]],
                        items[col_inds["mBase"]],
                        round(Int, items[col_inds["status"]]),
                        items[col_inds["Pmax"]],
                        items[col_inds["Pmin"]],
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        row_num,
                    ))
                elseif cur_key == "gencost"
                    first_cost_col = col_inds["n"] + 1
                    print(col_inds)
                    gen[row_num] = GenData(
                        round(Int, gen[row_num].bus),
                        gen[row_num].pg,
                        gen[row_num].qg,
                        gen[row_num].qmax,
                        gen[row_num].qmin,
                        gen[row_num].vg,
                        gen[row_num].mBase,
                        gen[row_num].status,
                        gen[row_num].pmax,
                        gen[row_num].pmin,
                        items[col_inds["startup"]],
                        items[col_inds["shutdown"]],
                        round(Int, items[col_inds["n"]]),
                        Tuple(items[first_cost_col:first_cost_col+2]),
                        row_num,
                    )
                elseif cur_key == "branch"
                    ratea_i = round(Int, items[col_inds["rateA"]])
                    r = items[col_inds["r"]]
                    xraw = items[col_inds["x"]]
                    braw = items[col_inds["b"]]
                    x = r + im * xraw
                    xi = inv(x)
                    y = ifelse(isfinite(xi), xi, zero(xi))
                    g = real(y)
                    b = imag(y)
                    # tap / shift take default vals (1, 0 respectively) in pglib
                    ## tr = branch.tap .* cos.(branch.shift)
                    ## ti = branch.tap .* sin.(branch.shift)
                    ## ttm = tr^2 + ti^2
                    ## g_fr = branch["g_fr"]
                    ## b_fr = branch["b_fr"]
                    ## g_to = branch["g_to"]
                    ## b_to = branch["b_to"]
                    ## branch.c1 = (-g * tr - b * ti) / ttm
                    ## branch.c2 = (-b * tr + g * ti) / ttm
                    ## branch.c3 = (-g * tr + b * ti) / ttm
                    ## branch.c4 = (-b * tr - g * ti) / ttm
                    ## branch.c5 = (g + g_fr) / ttm
                    ## branch.c6 = (b + b_fr) / ttm
                    ## branch.c7 = (g + g_to)
                    ## branch.c8 = (b + b_to)
                    # g_fr, g_to, b_fr, b_to are all not in pglib
                    b_fr :: T = braw / 2.0
                    b_to :: T = b_fr
                    push!(branch, BranchData(
                        round(Int, items[col_inds["fbus"]]),
                        round(Int, items[col_inds["tbus"]]),
                        r,
                        xraw,
                        braw,
                        ratea_i,
                        round(Int, items[col_inds["rateB"]]),
                        round(Int, items[col_inds["rateC"]]),
                        items[col_inds["ratio"]],
                        items[col_inds["angle"]],
                        round(Int, items[col_inds["status"]]),
                        items[col_inds["angmin"]],
                        items[col_inds["angmax"]],
                        ratea_i ^ 2,
                        -g,
                        -b,
                        -g,
                        -b,
                        g,
                        b + b_fr,
                        g,
                        b + b_to
                    ))
                elseif cur_key == "storage"
                    push!(storage, StorageData(
                        items[col_inds["storage_bus"]],
                        items[col_inds["ps"]],
                        items[col_inds["qs"]],
                        items[col_inds["energy"]],
                        items[col_inds["energy_rating"]],
                        items[col_inds["charge_rating"]],
                        items[col_inds["discharge_rating"]],
                        items[col_inds["charge_efficiency"]],
                        items[col_inds["discharge_efficiency"]],
                        items[col_inds["thermal_rating"]],
                        items[col_inds["qmin"]],
                        items[col_inds["qmax"]],
                        items[col_inds["r"]],
                        items[col_inds["x"]],
                        items[col_inds["p_loss"]],
                        items[col_inds["q_loss"]],
                        items[col_inds["status"]],
                    ))
                    set_fields!(storage[row_num], items, 1)
                end
                row_num += 1
            end
        elseif length(line) != 0 && line[1] != '%' && words[1] != "function"
            col_inds = Dict()
            columns = [m.match for m in eachmatch(col_patt, comment)][2:end]
            comment = ""
            for (i, column) in enumerate(columns)
                merge!(col_inds, Dict(column => i))
            end
            print(col_inds)
            cur_key = ""
            type = Any

            for cur in MATPOWER_VAR_TYPES
                full_name = "mpc.$(cur.name)"
                idxs = findall(s -> s == full_name, words) 
                if idxs != []
                    cur_key = cur.name
                    type = cur.type
                    break
                end
            end

            if cur_key == ""
                error("Error parsing data. Invalid variable assignment on line $(line_ind).")
            end
            if cur_key == "version"
                raw_data = words[length(words)-1]
                version = String(raw_data[2:length(raw_data)-1])
            elseif cur_key == "baseMVA"
                baseMVA = parse(Float64, words[length(words)-1])
            else
                in_array = true
                row_num = 1
                line = String(join(words[findall(s -> s == "[", words)[1]+1:end], " "))
                continue
            end
        elseif length(line) != 0 && line[1] == '%'
            comment = line
        end

        if line_ind < length(lines)
            line = lines[line_ind += 1]
        else
            break
        end
    end

    arc_from = [(i, b.fbus, b.tbus) for (i, b) in enumerate(branch)]
    arc_to = [(i, b.tbus, b.fbus) for (i, b) in enumerate(branch)]
    arc = [arc_from; arc_to]
    ref_buses = filter(i -> (bus[i]).type != 3, 1:length(bus))
    ratea = [branch[l].ratea for (l, _, _) in arc]

    return Data(version, baseMVA, bus, gen, branch, storage, ratea, arc, ref_buses)
end

function standardize_cost_terms!(data :: Data{T}, order) where T <: Real
    gen_order = 1
    for (_, gen) in enumerate(data.gen)
        max_ind = 1
        for i in 1:length(gen.c)
            max_ind = i
            if gen.c[i] != 0
                break
            end
            gen_order = max(gen_order, length(gen.c) - max_ind + 1)
        end
    end
    gen_order = max(gen_order, order + 1)
    for (i, gen) in enumerate(data.gen)
        if length(gen.c) == gen_order
            continue
        end
        std_cost = [0.0 for _ in 1:gen_order]
        cur_cost = reverse(gen.c)
        for i in 1:min(gen_order, length(cur_cost))
            std_cost[i] = cur_cost[i]
        end
        c = reverse(std_cost)
        n = length(gen_order)
        data.gen[i] = GenData(
            data.gen[i].bus,
            data.gen[i].pg,
            data.gen[i].qg,
            data.gen[i].qmax,
            data.gen[i].qmin,
            data.gen[i].vg,
            data.gen[i].mBase,
            data.gen[i].status,
            data.gen[i].pmax,
            data.gen[i].pmin,
            data.gen[i].startup,
            data.gen[i].shutdown,
            n,
            c,
            data.gen[i].i
        )
    end
end

function calc_thermal_limits!(data :: Data{T}) where T <: Real
    for branch in filter(branch -> branch.ratea <= 0, data.branch)
        xi = inv(branch.r + im * branch.x)
        y_mag = abs.(ifelse(isfinite(xi), xi, zero(xi)))

        fr_vmax = data.bus[branch.fbus].Vmax
        to_vmax = data.bus[branch.tbus].Vmax
        m_vmax = max(fr_vmax, to_vmax)

        theta_max = max(abs(branch.angmin), abs(branch.angmax))
        c_max = sqrt(fr_vmax^2 + to_vmax^2 - 2*fr_vmax*to_vmax*cos(theta_max))

        branch.rateA = y_mag * m_vmax * c_max
    end
end

function parse_ac_power_data(::Type{T}, filename) :: Data{T} where T <: Real
    _, f = splitdir(filename)
    name, ext = splitext(f)
    @info "hi"

    if isfile(joinpath(TMPDIR, name) * ".jld2")
        @info "Loading cached JLD2 file"
        loaded = JLD2.load(joinpath(TMPDIR, name) * ".jld2")
        return loaded["data"]
    else
        ff = if isfile(filename)
            filename
        elseif isfile(joinpath(TMPDIR, name) * ".m")
            joinpath(TMPDIR, name) * ".m"
        else
            @info "Downloading $filename"
            Downloads.download(
                "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/master/$filename",
                joinpath(TMPDIR, name * ".m"),
            )
            joinpath(TMPDIR, name * ".m")
        end
        @info "Loading MATPOWER file"
        return process_ac_power_data(T, ff)
    end
end

function process_ac_power_data(::Type{T}, filename) :: Data{T} where T <: Real
    data = parse_matpower_file(T, filename)
    standardize_cost_terms!(data, 2)
    calc_thermal_limits!(data)


    @info "Saving JLD2 cache file"
    _, f = splitdir(filename)
    name, _ = splitext(f)
    # JLD2.save(joinpath(TMPDIR, name * ".jld2"), "data", data)
    
    return data
end

function struct_to_nt(data :: T) where T <: DataType
    result = ()
    for field in fieldnames(T)
        field_val = getfield(field, data)
        if field_val <: DataType
            field_val = struct_to_nt(field_val)
        end
        result = merge((field = field_val), result)
    end
    return result
end
