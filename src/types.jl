using CSV, DataFrames, Rematch, MATLAB, StaticArrays, CoordinateTransformations, Dierckx, Unitful


abstract type AbstractExperiment end

abstract type PresentOrgNest <: AbstractExperiment end
abstract type AbsentOrgNest <: AbstractExperiment end

######################### Experiments ##############

struct ClosedNest <: AbsentOrgNest
    path::String
    name::String
end

struct TransferNest <: PresentOrgNest
    path::String
    name::String
end

struct Transfer <: AbsentOrgNest
    path::String
    name::String
end

struct TransferNestBelen <: AbsentOrgNest
    path::String
    name::String
end

struct DawaySandpaper <: PresentOrgNest
    path::String
    name::String
end

struct Daway <: PresentOrgNest
    path::String
    name::String
end

function AbstractExperiment(path::String) 
    _, name = splitdir(path)
    @match name begin 
        "Closed nest Therese" || "Closed nest_Belén"                                                => ClosedNest(path, name)
        "Transfer nest Therese"                                                                     => TransferNest(path, name)
        "Transfer Therese" || "Transfer pellet moved Therese"  || "Transfer pellet removed Therese" => Transfer(path, name)
        "Transfer_Belén"                                                                            => TransferNestBelen(path, name)
        "Daway_Belén" || "Awayd_closed_nest_Belén" || "Dright_Belén" || "Dleft_closednest_Belén" || "Dzerovector_Rot_Elin" || "Daway_Rot_Elin" || "Daway_S_Elin" || "Dforward_S_Elin" || "Dleft_S_Elin" || "Dright_S_Elin" => DawaySandpaper(path, name)
        "Daway Therese" || "Dleft Therese" || "Dtowards Therese" || "Dright Therese"                => Daway(path, name)
    end
end

# filter!(x -> x ∉ bad, runs)

const Coordinate = SVector{2, Float64}

struct Point
    xy::Coordinate
    i::Int
end

struct Final
    name::String
    track::Vector{Coordinate}
    time::StepRangeLen
    angle::Vector{Float64}
    length::Vector{Float64}
    speed::Vector{Float64}
    pellet::Union{Missing, Point}
    turn_point::Point
    nest::Coordinate
    feeder::Coordinate
    org_nest::Coordinate
end

struct Experiment{E <: AbstractExperiment}
    exp_type::E
    runs::Vector{Final}
end

struct Basic{E<:AbstractExperiment}
    exp_type::E
    name::String
    path::String
    Basic{E}(exp_type::E, name::String) where {E<:AbstractExperiment} = new(exp_type, name, joinpath(exp_type.path, name))
end
Basic(exp_type::E, name::String) where {E<:AbstractExperiment} = Basic{E}(exp_type, name)

#=function get_experiment(path::String)
    name, run = splitdir(path)
    _, _experiment = splitdir(name)
    match_experiment(_experiment, run, path)
end=#

# get_experiment_type(d::Final{E, <:FinalPellet}) where {E <: Experiment} = E

function Experiment(path::String)
    exp_type = AbstractExperiment(path)
    runs = Final[]
    for r in readdir(path)
        if r"^[^.]"(r)
            a = Basic(exp_type, r)
            if isdir(a.path) 
                if CHECK
                    msgs = check(a)
                    if !isempty(msgs)
                        print_with_color.(:magenta, "Errors in $(a.exp_type.name)/$(a.name):\n")
                        print_with_color.(:red, msgs, '\n')
                        continue
                    end
                end
                push!(runs, Final(a))
            end
        end
    end
    Experiment(exp_type, runs)
end

abstract type Calibration end

struct Calib <: Calibration
    file::String
end
struct JPG <: Calibration
    file::String
    checker::Float64 # in cm
end
struct MATcalib <: Calibration
    file::String
end

const factor_filename = "factors.csv"

function get_calibration(path, calibfile)
    name, ext = splitext(calibfile)
    return @match ext begin
        ".calib" => Calib(joinpath(path, calibfile))
        ".mat" => MATcalib(joinpath(path, calibfile))
        ".jpg" || ".JPG" where r"calib"(name) => begin
            factors = CSV.read(joinpath(path, factor_filename), header=["Factor", "Level"], types=[String, Union{String, Missing}], transforms = Dict(1 => strip))
            d, u = split(first(factors[factors[:Factor] .== "check", :Level]))
            checker = ustrip(uconvert(Unitful.cm, parse(d)*getfield(Unitful, Symbol(u))))
            JPG(joinpath(path, calibfile), checker)
        end
    end
end
get_calibrated_data(path::String, calib::Calib, mtlbs::Vector{Dict{String, Any}}) = mxcall(:calib_calibrate, 2, joinpath(path, calib.file), mtlbs)
get_calibrated_data(path::String, calib::JPG, mtlbs::Vector{Dict{String, Any}}) = mxcall(:jpg_calibrate, 2, joinpath(path, calib.file), calib.checker, mtlbs)
get_calibrated_data(path::String, calib::MATcalib, mtlbs::Vector{Dict{String, Any}}) = mxcall(:mat_calibrate, 2, joinpath(path, calib.file), mtlbs)

const poi_filename = "columns.csv"

function get_data(a::Basic{<:AbstractExperiment})
    fname = joinpath(a.path, poi_filename)
    pois = CSV.read(fname, header=["resfile", "columnnumber", "poi"], types=[String, Int, String], transforms = Dict(1 => strip, 3 => strip))
    # pellet = "pellet" ∈ pois[:poi] ? BasicPresent(pois) : BasicAbsent(pois)
    # Basic(experiment, pellet)
    fname = joinpath(a.path, "calibrations.csv")
    calibs = CSV.read(fname, header=["resfile", "calibfile"], types=[String, String], transforms = Dict(1 => strip, 2 => strip))
    todo = join(calibs, pois, on = :resfile)
    data = Dict{Symbol,Tuple{Array{Float64,2},Array{Float64,1}}}()
    for ucalib in groupby(todo, :calibfile)
        mtlbs = Dict{String, Any}[]
        for ures in groupby(ucalib, :resfile)
            mtlb = Dict{String, Any}()
            mtlb["resfile"] = joinpath(a.path, first(ures[:resfile]))
            mtlb["pois"] = ures[:columnnumber]
            push!(mtlbs, mtlb)
        end
        calib = get_calibration(a.path, first(ucalib[:calibfile]))
        xys, ts = get_calibrated_data(a.path, calib, mtlbs)
        for (poi, xy, t) in zip(ucalib[:poi], xys, ts)
            _t = t isa Float64 ? [t] : t
            data[Symbol(poi)] = (xy, _t)
        end
    end
    data
end

const Coordinate3D = SVector{3, Float64}

# abstract type DataPellet end

# struct DataAbsent <: DataPellet
# end

struct DataPoint
    xy::Coordinate3D
    t::Float64
end

struct Data{E<:AbstractExperiment}
    exp_type::E
    name::String
    pellet::Union{DataPoint, Missing}
    track::Vector{Coordinate3D}
    time::Vector{Float64}
    data::Dict{Symbol,Tuple{Matrix{Float64},Vector{Float64}}}
end

function retrievetrack(_track::Matrix{Float64}, _time::Vector{Float64})
    track = [Coordinate3D(_track[i,1], _track[i,2], 0.0) for i in 1:size(_track, 1)]
    time = vec(_time)
    track, time
end

function retrievepoint(data::Dict{Symbol,Tuple{Matrix{Float64},Vector{Float64}}}, f::Symbol)
    x, t = data[f]
    Coordinate3D(x[1], x[2], 0.0), t[1]
end

function Data(a::Basic{<:AbstractExperiment})
    data = get_data(a)
    track, time = retrievetrack(data[:track]...)
    pellet = haskey(data, :pellet) ? DataPoint(retrievepoint(data, :pellet)...) : missing
    Data(a.exp_type, a.name, pellet, track, time, data)
end

struct Common
    name::String
    pellet::Union{DataPoint, Missing}
    track::Vector{Coordinate3D}
    time::Vector{Float64}
    feeder::Coordinate3D
    nest::Coordinate3D
    org_nest::Coordinate3D
end

######################### Common methods ###########

function Common(b::Data{ClosedNest})
    feeder, _ = retrievepoint(b.data, :feeder)
    nest, _ = retrievepoint(b.data, :nest)
    Common(b.name, b.pellet, b.track, b.time, feeder, nest, nest)
end

function Common(b::Data{TransferNestBelen})
    southbefore, _ = retrievepoint(b.data, :southbefore)
    northbefore, _ = retrievepoint(b.data, :northbefore)
    feederbefore, _ = retrievepoint(b.data, :feederbefore)
    nestbefore, _ = retrievepoint(b.data, :nestbefore)
    v = northbefore - southbefore
    u = nestbefore - feederbefore 
    azimuth = atan2(v[2], v[1]) - atan2(u[2], u[1])
    nest2feeder = norm(nestbefore - feederbefore)

    south, _ = retrievepoint(b.data, :south)
    north, _ = retrievepoint(b.data, :north)
    v = north - south
    α = atan2(v[2], v[1]) + azimuth
    u = Coordinate3D(cos(α), sin(α), 0)
    feeder, _  = retrievepoint(b.data, :feeder)
    nest = feeder + u*nest2feeder
    Common(b.name, b.pellet, b.track, b.time, feeder, nest, Coordinate3D(NaN, NaN, NaN))
end

function Common(b::Data{Transfer})
    south, _ = retrievepoint(b.data, :south)
    north, _ = retrievepoint(b.data, :north)
    factors = CSV.read(joinpath(b.exp_type.path, b.name, factor_filename), header=["Factor", "Level"], types=[String, Union{String, Missing}], transforms = Dict(1 => strip))
    d, u = split(first(factors[factors[:Factor] .== "nest2feeder", :Level]))
    nest2feeder = ustrip(uconvert(Unitful.cm, parse(d)*getfield(Unitful, Symbol(u))))
    d, u = split(first(factors[factors[:Factor] .== "azimuth", :Level]))
    azimuth = ustrip(uconvert(Unitful.rad, parse(d)*getfield(Unitful, Symbol(u))))
    v = north - south
    α = atan2(v[2], v[1]) + azimuth - π
    u = Coordinate3D(cos(α), sin(α), 0)
    feeder, _ = retrievepoint(b.data, :feeder)
    nest = feeder + u*nest2feeder
    Common(b.name, b.pellet, b.track, b.time, feeder, nest, Coordinate3D(NaN, NaN, NaN))
end

function Common(b::Data{TransferNest})
    btransfer = Data(Transfer(b.exp_type.path, b.exp_type.name), b.name, b.pellet, b.track, b.time, b.data)
    ctransfer = Common(btransfer)
    org_nest, _ = retrievepoint(b.data, :originalnest)
    Common(b.name, ctransfer.pellet, ctransfer.track, ctransfer.time, ctransfer.feeder, ctransfer.nest, Coordinate3D(org_nest))
end

function Common(b::Data{DawaySandpaper})
    org_nest, _ = retrievepoint(b.data, :nest)
    initial = mean(first(retrievepoint(b.data, k)) for k in [:rightdowninitial, :leftdowninitial, :rightupinitial, :leftupinitial])
    final = mean(first(retrievepoint(b.data, k)) for k in [:rightdownfinal, :leftdownfinal, :rightupfinal, :leftupfinal])
    v = final - initial
    nest = org_nest + v
    _feeder, _ = retrievepoint(b.data, :feeder)
    feeder = _feeder + v
    Common(b.name, b.pellet, b.track, b.time, feeder, nest, org_nest)
end

function Common(b::Data{Daway})
    org_nest, _ = retrievepoint(b.data, :nest)
    initialfeeder, _ = retrievepoint(b.data, :initialfeeder)
    feeder, _ = retrievepoint(b.data, :feeder)
    v = feeder - initialfeeder
    nest = org_nest + v
    Common(b.name, b.pellet, b.track, b.time, feeder, nest, org_nest)
end

##################################################

function get_standardization(nest::Coordinate3D, feeder::Coordinate3D)
    trans = Translation(-vec(nest))
    v = feeder - nest
    rot = LinearMap(RotZ(-atan2(v[2], v[1]) - π/2))
    rot∘trans
    # f(x::Coordinate3D) = Coordinate(composed(x)[1:2])
    # return f
end

standardize(f::T, pellet::DataPoint, t1) where {T <: AbstractAffineMap} = DataPoint(f(pellet.xy), pellet.t - t1)
standardize(f::T, pellet::Missing, t1) where {T <: AbstractAffineMap} = pellet
const atol = 1e-11

function standardize(f::T, c::Common) where {T<:AbstractAffineMap}
    track = f.(c.track)
    t1 = c.time[1]
    time = c.time - t1
    pellet = standardize(f, c.pellet, t1)
    feeder = f(c.feeder)
    nest = f(c.nest)
    org_nest = f(c.org_nest)
    @assert isapprox(feeder[1], 0, atol = atol) "feeder x value: $feederX > $atol"
    return Common(c.name, pellet, track, time, feeder, nest, org_nest)
end

const smoothing_factor = 10
function smooth(track::Vector{Coordinate3D}, time::Vector{Float64})
    X = Array{Float64, 2}(2, length(track))
    for i in eachindex(track)
        X[:,i] = track[i][1:2]
    end
    Δt = mean(diff(time))
    ParametricSpline(time, X, s = 0.4smoothing_factor/Δt, k=1)
end
function outlier_fences(data, α = 0.2, q = 6)#3
    q1, q3 = quantile(data,[α, 1 - α])
    iqr = q3-q1
    major = q*iqr#10-14
    q3+major
end
function find_outlier(dang, i1, i)
    old = shift!(i)
    for new in i
        M = outlier_fences(dang[i1[1]:old[1]])
        dang[new] > M && return new
        old = new
    end
    return old
end
function findlocalmaxima(x::Vector{T}) where T <: Real
    inds = Int[]
    for i in 2:length(x) - 1
        if x[i - 1] < x[i] > x[i + 1]
            push!(inds, i)
        end
    end
    return inds
end
function detect_turningpoint(pathlength, angles)
    i1 = findfirst(pathlength .> 20)
    if i1 == 0
        i1 = 1
    end
    i2 = pathlength[end] ≤ 400 ? length(pathlength) : findfirst(pathlength .> 400)
    angs = angles[1:i2]
    dang = abs.(diff(angs))
    i = findlocalmaxima(dang)
    filter!(x -> x > i1, i)
    if isempty(i)
        push!(i, length(dang))
    end
    j = find_outlier(dang, i1, i)
    μ = atan2(sum(sin.(angles[i1:j])), sum(cos.(angles[i1:j])))
    δ = abs.(angles[i1:end] - μ)
    turningpoint_i = findfirst(δ .> π/4) + i1 - 1
    if turningpoint_i == 0
        turningpoint_i = 1
    end
    return turningpoint_i
end

function common(track, time)
    spl = smooth(track, time)
    n = 2length(time)
    times = linspace(time[1], time[end], n)
    tracks = [Coordinate(evaluate(spl, t)...) for t in times]
    angles = Vector{Float64}(n)
    pathlengths = Vector{Float64}(n)
    for i in eachindex(angles)
        dx, dy = derivative(spl, times[i])
        angles[i] = atan2(dy, dx)
        pathlengths[i] = sqrt(dx^2 + dy^2)
    end
    pathlengths .= cumsum(pathlengths*step(times))
    turningpoint_i = detect_turningpoint(pathlengths, angles)
    speeds = pathlengths./times
    speeds[1] = speeds[2]
    # range = 1:turningpoint_i
    # feeder2turning = Coordinate.(tracks[range], times[range], angles[range], pathlengths[range], speeds[range])
    tracks, times, angles, pathlengths, speeds, turningpoint_i
end

# abstract type FinalPellet end

# struct FinalAbsent <: FinalPellet
# end

# struct FinalPresent <: FinalPellet
# xy::Coordinate
# i::Int
# end

function find_pellet_i(pellet::DataPoint, time, track::Vector{Coordinate})
    m = Inf
    ind = 0
    for i in eachindex(time) 
        if time[i] > pellet.t - 2 
            _m = norm(track[i] - pellet.xy[1:2])
            if _m < m
                m = _m
                ind = i
            end
        elseif time[i] > pellet.t + 2
            break
        end
    end
    return Point(Coordinate(pellet.xy[1:2]), ind)
end
find_pellet_i(pellet::Missing, time, track::Vector{Coordinate}) = pellet

larger(i::Int, p::Point) = i ≥ p.i
larger(i::Int, p::Missing) = false
function Final(c::Common)
    f = get_standardization(c.nest, c.feeder)
    d = standardize(f, c)
    track, time, angle, pathlength, speed, turningpoint_i = common(d.track, d.time)
    pellet = find_pellet_i(d.pellet, time, track)
    turningpoint_i = larger(turningpoint_i, pellet) ? pellet.i - 1 : turningpoint_i
    Final(c.name, track, time, angle, pathlength, speed, pellet, Point(track[turningpoint_i], turningpoint_i), Coordinate(d.nest[1:2]), Coordinate(d.feeder[1:2]), Coordinate(d.org_nest[1:2])) 
end

Final(a::Basic) = Final(Common(Data(a)))
