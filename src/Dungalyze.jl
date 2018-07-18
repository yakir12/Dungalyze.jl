module Dungalyze

# package code goes here
# using StaticArrays, Dierckx, CSV, DataFrames, MAT, MATLAB, CoordinateTransformations, Plots, Unitful, Rematch, Distributions, Missings
# include(joinpath(Pkg.dir("Dungalyze"), "src", "calibrations.jl"))
# include(joinpath(Pkg.dir("Dungalyze"), "src", "turning_point.jl"))



include(joinpath(Pkg.dir("Dungalyze"), "src", "types.jl"))
include(joinpath(Pkg.dir("Dungalyze"), "src", "plots.jl"))
include(joinpath(Pkg.dir("Dungalyze"), "src", "checks.jl"))

home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
figfolder = "dungalyze figures"
if !isdir(figfolder)
    mkdir(figfolder)
end

const CHECK = true
for experiment in readdir(home)
    # if !r"elin"i(experiment)
        path = joinpath(home, experiment)
        a = Experiment(path)
        plot(a)
    # end
end


#=df = DataFrame(Experiment = "", Run = "", Feeder_y = 0.0, Turning_point_x = 0.0, Turning_point_y = 0.0, Pellet_x = 0.0, Pellet_y = 0.0, Track_end_x = 0.0, Track_end_y = 0.0, ORG_nest_x = 0.0, ORG_nest_y = 0.0, Gravity_turn2pellet_x = 0.0, Gravity_turn2pellet_y = 0.0, Gravity_pellet2end_x = 0.0, Gravity_pellet2end_y = 0.0, Gravity_x = 0.0, Gravity_y = 0.0, Dist2turn = 0.0, Dis2pellet = 0.0, Dis2end = 0.0, Dis2org_nest = 0.0, Dis2grav_turn2pellet = 0.0, Dis2grav_pellet2end = 0.0, Dis2gravity = 0.0)
allowmissing!(df)
function doall(experiment::String; bad::Vector{String} = String[])
    _, _experiment = splitdir(experiment)
    mkdir(_experiment)
    cd(_experiment)
    ds = Final[]
    runs = get_runs(experiment, bad=bad)
    for r in runs
        d = Final(r)
        push!(df, row(d))
        push!(ds, d)
    end
    plot(ds)
    cd("..")
end
for experiment in readdir(home)
    bad = problems[problems[:Experiment] .== experiment, :Run]
    doall(joinpath(home, experiment), bad = bad)
end
deleterows!(df, 1)
CSV.write("all.csv", df)

experiment = "/home/yakir/google-drive/.shared/Coffee beetle experiments/Transfer pellet moved Therese"


for experiment in readdir(home)
    a = Experiment(joinpath(home, experiment))
    plot(a)
end


for experiment in readdir(home)
    if r"elin"i(experiment)
        f = joinpath(home, experiment)
        # println(f)
        a = Experiment(f)
        plot(a)
    end
end





runs = get_runs(experiment)
experiment = get_experiment(runs[1])
ds = Final{typeof(experiment), <:FinalPellet}[]

for r in runs
    d = Final(r)
    push!(ds, d)
end

plot(ds)


=##=function getpois(path::String)
    b = CSV.read(joinpath(path, poi_filename), header=["File", "Number", "POI"], types=[String, Int, String])
    b[:POI] .= strip.(b[:POI])
    pois = Dict{String, Dict{String, Int}}()
    by(b, :File) do df
        pois[df[1, :File]] = Dict(df[i,:POI] => df[i,:Number] for i in 1:nrow(df))
        end
        return pois
    end=##=
    =##=function extract_raw!(data::Dict{Symbol,Tuple{Matrix{Float64},Vector{Float64}}}, pois::Dict{String, Dict{String, Int}}, path::String, file::String)
        v = pois[file]
        matopen(joinpath(path, file), "r") do i
            X, Y, T, tmp = read(i, "xdata", "ydata", "tdata", "status")
            fr = tmp["FrameRate"]
            for (poi, column) in v
                data[Symbol(poi)] = _extract_raw(X[:, column], Y[:, column], T[:, column], fr)
            end
        end
    end=##=
    function doall(experiment::String; bad::Vector{String} = String[])
        runs = get_runs(experiment, bad=bad)
        ds = Final[]
        for r in runs
            a = sort_experiment(r)
            # d = Final(a)
            b = Calibrated(a)
            # c = Standardized(b)
            d = get_final(b)
            # p1 = plot(d)
            # p2 = plotsmooth(c, d)
            m = match(r"/.shared/(.*)$", r)
            path = joinpath("/home/yakir/tmp", m.captures[1])
            # mkpath(path)
            # savefig(p1, joinpath(path, "summary.pdf"))
            # savefig(p2, joinpath(path, "smooth.pdf"))
            push!(ds, d)
        end
        =##=if length(ds) > 1
            p = plot(ds)
            m = match(r"/.shared/(.*)$", experiment)
            path = joinpath("/home/yakir/tmp", m.captures[1])
            savefig(p, joinpath(path, "summary.pdf"))
        end=##=
    end


    # end # module



    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    a = Dict("Awayd_closed_nest_Belén" => ["64"],
             "Closed nest_Belén" => ["24"],
             "Daway_Belén" => [""],
             "Dright_Belén" => [""],
             "Dleft_closednest_Belén" => [""],
             "Transfer_Belén" => [""])
    for (experiment, bad) in a
        println(experiment)
        doall(joinpath(home, experiment), bad = bad)
    end



    doall.(joinpath.(home, readdir(home)))

    doall(joinpath(home, "Transfer pellet moved Therese"))

    experiments = filter(x -> !(r"^Transfer pellet moved Therese$"(x)), readdir(home))
    doall.(joinpath.(home, experiments))


    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    doall(joinpath(home, "Transfer_Belén"))

    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    doall(joinpath(home, "Transfer nest Therese"))

    experiments = filter(r"Bel[e,é]n", readdir(home))
    doall.(joinpath.(home, experiments))

    experiments = filter(r"therese"i, readdir(home))
    doall.(joinpath.(home, experiments))

    experiment = joinpath(home, "Closed nest_Belén")


    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    experiment = joinpath(home, filter(r"therese"i, readdir(home))[6])
    runs = get_runs(experiment)
    r = runs[1]
    a = sort_experiment(r)

    b = Calibrated(a)
    d = Final(b)



    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    include(joinpath(Pkg.dir("Dungalyze"), "src", "checks.jl"))

    for experiment in readdir(home)
        bad = problems[problems[:Experiment] .== experiment, :Run]
        doall(joinpath(home, experiment), bad = bad)
    end

    home = "/home/yakir/google-drive/.shared/Coffee beetle experiments"
    experiment = joinpath(home, "Closed nest Therese")
    bad = [""]



    experiment_path = joinpath(home, experiment)
    all_runs = get_runs(experiment_path)
    bad_runs = joinpath.(experiment_path, problems[problems[:Experiment] .== experiment, :Run])
    runs = setdiff(all_runs, bad_runs)
    for r in runs
        a = sort_experiment(r)
        b = Calibrated(a)
        d = Final(b)=#


end # module
