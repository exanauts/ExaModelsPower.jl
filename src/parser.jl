convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(convert_array(d, backend) for d in data)

Base.@kwdef mutable struct BusData
    bus_i :: Int = 0
    type :: Int = 0
    pd :: Float64 = 0
    qd :: Float64 = 0
    gs :: Float64 = 0
    bs :: Float64 = 0
    area :: Float64 = 0
    Vm :: Float64 = 0
    Va :: Float64 = 0
    baseKV :: Float64 = 0
    zone :: Float64 = 0
    vmax :: Float64 = 0
    vmin :: Float64 = 0
end
Base.@kwdef mutable struct BranchData
    fbus :: Int = 0
    tbus :: Int = 0
    r :: Float64 = 0
    x :: Float64 = 0
    b :: Float64 = 0
    ratea ::Int = 0
    rateb :: Int = 0
    ratec :: Int = 0
    ratio :: Float64 = 0
    angle :: Float64 = 0
    status :: Int = 0
    angmin :: Float64 = 0
    angmax :: Float64 = 0
    ratea_sq :: Int = 0
    c1 :: Float64 = 0
    c2 :: Float64 = 0
    c3 :: Float64 = 0
    c4 :: Float64 = 0
    c5 :: Float64 = 0
    c6 :: Float64 = 0
    c7 :: Float64 = 0
    c8 :: Float64 = 0
end
Base.@kwdef mutable struct StorageData
    storage_bus :: Float64
    ps :: Int
    qs :: Float64
    energy :: Float64
    energy_rating :: Float64
    charge_rating :: Float64
    discharge_rating :: Float64
    charge_efficiency :: Float64
    discharge_efficiency :: Float64
    thermal_rating :: Float64
    qmin :: Float64
    qmax :: Float64
    r :: Float64
    x :: Float64
    p_loss :: Float64
    q_loss :: Float64
    status :: Int
end
Base.@kwdef mutable struct GenData
    bus :: Int = 0
    pg :: Float64 = 0
    qg :: Float64 = 0
    qmax :: Float64 = 0
    qmin :: Float64 = 0
    vg :: Float64 = 0
    mBase :: Float64 = 0
    status :: Int = 0
    pmax :: Float64 = 0
    pmin :: Float64 = 0
    _2 :: Int = 0
    startup :: Float64 = 0
    shutdown :: Float64 = 0
    n :: Int = 0
    c :: Vector{Float64} = Vector{Float64}()
end

const MatPowerRow = Union{BusData, BranchData, StorageData, GenData}
function set_fields!(x :: MatPowerRow, row :: Vector{Float64}, start :: Int64)
    for (i, field) in enumerate(fieldnames(typeof(x))[start:end])
        if i > length(row)
            return
        end
        if typeof(getfield(x, field)) == Int
            setfield!(x, field, round(Int, row[i]))
        else
            setfield!(x, field, row[i])
        end
    end
end

Base.@kwdef mutable struct Data
    version :: String = ""
    baseMVA :: Float64 = 0
    bus :: Vector{BusData} = Vector{BusData}()
    gen :: Vector{GenData} = Vector{GenData}()
    branch :: Vector{BranchData} = Vector{BranchData}()
    storage :: Vector{StorageData} = Vector{StorageData}()
    ratea :: Vector{Float64} = Vector{Float64}()
    arc :: Vector{Tuple{Int, Int, Int}} = Vector{Tuple{Int, Int, Int}}()
    ref_buses :: Vector{Int} = Vector{Int}()
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

function parse_matpower_file(fname :: String) :: Data
    lines = split(read(open(fname), String), "\n")
    in_array = false
    cur_key = ""
    data :: Data = Data()
    line_ind = 1
    line = lines[line_ind]
    pattern = r"[^\s=;\[\]]+|[=;\[\]]"
    type = Missing
    row_num = 1

    while true
        line = split(line, "%")[1]
        words = [m.match for m in eachmatch(pattern, line)]

        if in_array
            squares = findall(s -> s == "]", words)
            first_sq = length(squares) == 0 ? typemax(Int64) : squares[1]
            first_semi = length(squares) == 0 ? typemax(Int64) : squares[1]
            items :: Vector{Float64} = map(s -> parse(Float64, s), words[1:min(min(first_semi - 1, first_sq - 1), length(words) - 1)])

            if length(items) == 0 && length(words) >= 2 && words[length(words)-1] == "]" && words[length(words)] == ";"
                in_array = false
            elseif length(words) != 0 && last(words) != ";"
                error("Invalid matpower file. Line $(line_ind) array doesn't end with ; or ];")
            elseif length(items) != 0
                if cur_key == "bus"
                    push!(data.bus, BusData())
                    set_fields!(data.bus[row_num], items, 1)
                elseif cur_key == "gen"
                    push!(data.gen, GenData())
                    set_fields!(data.gen[row_num], items, 1)
                elseif cur_key == "gencost"
                    data.gen[row_num].c = items[5:end]
                    set_fields!(data.gen[row_num], items[1:4], 11)
                elseif cur_key == "branch"
                    push!(data.branch, BranchData())
                    set_fields!(data.branch[row_num], items, 1)
                    branch = last(data.branch)
                    branch = data.branch[length(data.branch)]
                    branch.ratea_sq = branch.ratea_sq ^ 2
                    x = branch.r + im * branch.x
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
                    b_fr = branch.b / 2.0
                    b_to = b_fr
                    branch.c1 = -g
                    branch.c2 = -b
                    branch.c3 = -g
                    branch.c4 = -b
                    branch.c5 = g
                    branch.c6 = b + b_fr
                    branch.c7 = g
                    branch.c8 = b + b_to
                elseif cur_key == "storage"
                    push!(data.storage, StorageData())
                    set_fields!(data.storage[row_num], items, 1)
                end
                row_num += 1
            end

        elseif length(line) != 0 && line[1] != '%' && words[1] != "function"
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
                data.version = String(raw_data[2:length(raw_data)-1])
            elseif cur_key == "baseMVA"
                data.baseMVA = parse(Float64, words[length(words)-1])
            else
                in_array = true
                row_num = 1
                line = String(join(words[findall(s -> s == "[", words)[1]+1:end], " "))
                continue
            end
        elseif length(line) != 0 && line[1] == '%'
        end

        if line_ind < length(lines)
            line = lines[line_ind += 1]
        else
            break
        end
    end
    return data
end

function standardize_cost_terms!(data :: Data, order)
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
        gen.c = reverse(std_cost)
        gen.n = length(gen_order)
    end
end

function calc_thermal_limits!(data :: Data)
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

function parse_ac_power_data(filename) :: Data
    d, f = splitdir(filename)
    name, ext = splitext(f)

    if isfile(joinpath(TMPDIR, name) * ".jld2") && false
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
        return process_ac_power_data(ff)
    end
end

function process_ac_power_data(filename) :: Data
    data = parse_matpower_file(filename)
    standardize_cost_terms!(data, 2)
    calc_thermal_limits!(data)

    arc_from = [(i, b.fbus, b.tbus) for (i, b) in enumerate(data.branch)]
    arc_to = [(i, b.tbus, b.fbus) for (i, b) in enumerate(data.branch)]
    data.arc = [arc_from; arc_to]
    data.ref_buses = filter(i -> (data.bus[i]).type != 3, 1:length(data.bus))
    data.ratea = [data.branch[l].ratea for (l, _, _) in data.arc]

    @info "Saving JLD2 cache file"
    _, f = splitdir(filename)
    name, _ = splitext(f)
    JLD2.save(joinpath(TMPDIR, name * ".jld2"), "data", data)
    
    return data
end
