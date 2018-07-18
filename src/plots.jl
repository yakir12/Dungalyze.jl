using Plots, Distributions, Missings, DataStructures
import Plots:plot
# pgfplots()
gr()
# import Base:norm
#=get_org_nest(d::Final{<:Experiment, <:FinalPellet}) = (missing, missing), missing
get_org_nest(d::Final{<:PresentOrgNest, <:FinalPellet}) = d.org_nest, norm(d.org_nest)
get_gravitues(d::Final{<:Experiment, <:FinalPresent}) = mean(d.track[d.turningpoint_i:d.pellet.i]), mean(d.track[d.pellet.i:end]), mean(d.track[d.turningpoint_i:end])
get_gravitues(d::Final{<:Experiment, <:FinalAbsent}) = (missing,missing), (missing,missing), mean(d.track[d.turningpoint_i:end])
get_pellet_xy(pellet::FinalPresent) = pellet.xy
get_pellet_xy(pellet::FinalAbsent) = (missing, missing)
norm(x::Tuple{Missing, Missing}) = missing
function row(d::Final{<:Experiment, <:FinalPellet})
    org_nest, norm_org_nest = get_org_nest(d)
    gravity2pellet, gravity2end, gravity = get_gravitues(d)
    pelletxy = get_pellet_xy(d.pellet)
    d.experiment.experiment, d.experiment.run, d.feeder[2], d.track[d.turningpoint_i]..., pelletxy..., d.track[end]..., org_nest..., gravity2pellet..., gravity2end..., gravity..., norm(d.track[d.turningpoint_i]), norm(pelletxy), norm(d.track[end]), norm_org_nest, norm(gravity2pellet), norm(gravity2end), norm(gravity)
end=#

#=function row(d::Final{<:Experiment, FinalAbsent})
    org_nest, norm_org_nest = get_org_nest(d)
    gravity2pellet, gravity2end, gravity = get_gravitues(d)
    d.experiment.experiment, d.experiment.run, d.feeder[2], d.track[d.turningpoint_i]..., missing, missing, d.track[end]..., org_nest..., missing, missing, missing, missing, gravity..., norm(d.track[d.turningpoint_i]), missing, norm(d.track[end]), norm_org_nest, missing, missing, norm(gravity)
end=#


normp(x::Float64, y::Float64) = round(Int, sqrt(x^2 + y^2))
@recipe function f(::Type{Val{:feeder2turning}}, x,y,z)
    seriescolor := :blue
    seriestype := :path
    label := "feeder to turning p."
    x,y
end
@recipe function f(::Type{Val{:turning2pellet}}, x,y,z)
    seriescolor := :green
    seriestype := :path
    label := "turning p. to pellet"
    x,y
end
@recipe function f(::Type{Val{:pellet2end}}, x,y,z)
    seriescolor := :red
    seriestype := :path
    label := "pellet to end"
    x,y
end
@recipe function f(::Type{Val{:nest}}, x,y,z)
    label := "nest"
    markershape := :star5
    markersize := 4
    seriescolor := :yellow
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:feeder}}, x,y,z)
    l = normp(x[1], y[1])
    label := "feeder ($l cm)"
    markershape := :circle
    markersize := 4
    seriescolor := :brown
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:turningpoint}}, x,y,z)
    l = normp(x[1], y[1])
    label := "turning p. ($l cm)"
    markershape := :utriangle
    markersize := 4
    seriescolor := :lightgreen
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:pellet}}, x,y,z)
    l = normp(x[1], y[1])
    label := "pellet ($l cm)"
    markershape := :diamond
    markersize := 4
    seriescolor := :brown
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:trackend}}, x,y,z)
    l = normp(x[1], y[1])
    label := "end ($l cm)"
    markershape := :circle
    markersize := 6
    seriescolor := :black
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:gravity_turn2pellet}}, x,y,z)
    l = normp(x[1], y[1])
    label := "gravity t. to p. ($l cm)"
    markershape := :vline
    markersize := 4
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:gravity_pellet2end}}, x,y,z)
    l = normp(x[1], y[1])
    label := "gravity p. to end ($l cm)"
    markershape := :hline
    markersize := 4
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:gravity}}, x,y,z)
    l = normp(x[1], y[1])
    label := "gravity p. to end ($l cm)"
    # markershape := :cross
    markershape := :+
    markersize := 4
    seriestype := :scatter
    x,y
end
@recipe function f(::Type{Val{:turning2end}}, x,y,z)
    l = normp(x[1], y[1])
    label := " turn. p. to end ($l cm)"
    seriescolor := :green
    seriestype := :path
    x,y
end
@recipe function f(::Type{Val{:orgnest}}, x,y,z)
    l = normp(x[1], y[1])
    label := "org. nest ($l cm)"
    # markershape := :xcross
    markershape := :star8
    seriescolor := :yellow
    markersize := 4
    seriestype := :scatter
    x,y
end

function plot_basic(d::Final, p) 
    xy = d.track[1:d.turn_point.i]
    plot!(p, first.(xy), last.(xy), st = :feeder2turning)
    plot!(p, d.nest[1:1], d.nest[2:2], st = :nest)
    plot!(p, d.feeder[1:1], d.feeder[2:2], st = :feeder)
    xy = d.turn_point.xy
    plot!(p, xy[1:1], xy[2:2], st = :turningpoint)
    xy = d.track[end]
    plot!(p, xy[1:1], xy[2:2], st = :trackend)
    xy = mean(d.track[d.turn_point.i:end])
    plot!(p, xy[1:1], xy[2:2], st = :gravity)
end
function plot_extra(d::Final, p)
    if ismissing(d.pellet)
        xy = d.track[d.turn_point.i:end]
        plot!(p, first.(xy), last.(xy), st = :turning2end)
    else
        xy = d.track[d.turn_point.i:d.pellet.i]
        plot!(p, first.(xy), last.(xy), st = :turning2pellet)
        xy = d.track[d.pellet.i:end]
        plot!(p, first.(xy), last.(xy), st = :pellet2end)
        xy = d.pellet.xy
        plot!(p, xy[1:1], xy[2:2], st = :pellet)
        xy = mean(d.track[d.turn_point.i:d.pellet.i])
        plot!(p, xy[1:1], xy[2:2], st = :gravity_turn2pellet)
        xy = mean(d.track[d.pellet.i:end])
        plot!(p, xy[1:1], xy[2:2], st = :gravity_pellet2end)
    end
end
function plot_org(d::Final, p)
    xy = d.org_nest
    plot!(p, xy[1:1], xy[2:2], st = :orgnest)
end
function plot(::Type{E}, d::Final, exp_name::String) where {E <: AbstractExperiment}
    run_name = d.name
    p = plot(aspect_ratio=1, xlabel="X (cm)", ylabel="Y (cm)", title="$exp_name: $run_name", legend=:outertopright, background_color_legend = RGBA(1,1,1,0.25))#, size = (1000, 1000))
    plot_basic(d, p)
    plot_extra(d, p)
    E <: PresentOrgNest && plot_org(d, p)
    savefig(p, joinpath(figfolder, "$exp_name.$run_name.pdf"))
end


const win = 50
function plot_search(data::Matrix{Float64}, title::String)
    d = fit(DiagNormal, data)
    σ = sqrt.(var(d))
    fwhm = 2sqrt(2log(2))*σ
    p = histogram2d(data[1,:], data[2,:], bins=41)
    plot!([0,0], [-2win,2win], color=:white, label="")
    plot!([-2win,2win], [0,0], color=:white, label="")
    μx, μy = round.(Int, mean(d))
    scatter!(mean(d)[1:1], mean(d)[2:2], label="gravity (x, y): $μx cm, $μy cm", color=:green, markersize=8)
    t = linspace(0, 2pi, 100)'
    xy = mean(d) .+ fwhm/2.*[cos.(t); sin.(t)]
    fwhmx, fwhmy = round.(Int, fwhm)
    plot!(xy[1,:], xy[2,:], label="FWHM (x, y): $fwhmx cm, $fwhmy cm", color=:green, markersize=8)
    plot!(background_color_inside="black", aspect_ratio=1, xlim=(-win,win), ylim=(-win,win), colorbar=false, title=title, legend=false)
    return p
end
function get_search_data(ds::Vector{<:Final}, g::R, f::T) where {R<:Function, T<:Function}
    x = Float64[]
    y = Float64[]
    for d in ds
        for p in g(d)
            _x, _y = p - f(d)
            if -win < _x < win && -win < _y < win
                push!(x, _x)
                push!(y, _y)
            end
        end
    end
    return [x y]'
end

get_points(::Type{<:AbsentOrgNest}, ::Type{Any}) = OrderedDict{String, Any}("nest" => x -> x.nest, "turn. p." => x -> x.track[x.turn_point.i], "gravity of turn2end" => x -> mean(x.track[x.turn_point.i:end]))
function get_points(::Type{<:PresentOrgNest}, ::Type{Any})
    points = get_points(ClosedNest, Any)
    points["org_nest"] = x -> x.org_nest
    points
end
get_points(::Type{E}, Missing) where {E <: AbstractExperiment} = get_points(E, Any)
function get_points(::Type{E}, ::Type{Point}) where {E <: AbstractExperiment}
    points = get_points(E, Any)
    points["pellet"] = x -> x.pellet.xy
    points["gravity of turn2pellet"] = x -> mean(x.track[x.turn_point.i:x.pellet.i])
    points["gravity of pellet2end"] = x -> mean(x.track[x.pellet.i:end])
    points
end

get_segments(::Type{P}) where {P <: Any} = OrderedDict{String, Any}("turn2end" => x -> x.track[x.turn_point.i:end])
function get_segments(::Type{Point}) 
    segments = get_segments(Missing)
    segments["turn. p. to pellet"] = x -> x.track[x.turn_point.i:x.pellet.i]
    segments["pellet to end"] = x -> x.track[x.pellet.i:end]
    segments
end

savefigure(::Type{Any}, p, name) = savefig(p, joinpath(figfolder, "$name search density all.pdf"))
savefigure(::Type{Missing}, p, name) = savefig(p, joinpath(figfolder, "$name search density without pellet.pdf"))
savefigure(::Type{Point}, p, name) = savefig(p, joinpath(figfolder, "$name search density with pellet.pdf"))

function plot_density(::Type{E}, ds::Vector{Final}, name::String, n::Int) where {E<:AbstractExperiment}
    P = eltype(getfield.(ds, :pellet))
    points = get_points(E, P)
    segments = get_segments(P)
    segmentsn = length(segments)
    pointsn = length(points)
    ps = Plots.Plot[]
    for (ks, vs) in segments, (kp, vp) in points
        data = get_search_data(ds, vs, vp)
        push!(ps, plot_search(data, "$ks (centered on $kp) n=$n"))
    end
    p = plot(ps..., size=(700pointsn, 700segmentsn), layout=(segmentsn, pointsn))
    savefigure(P, p, name)
end


function plot(a::Experiment{E}) where {E <: AbstractExperiment}
    for d in a.runs
        plot(E, d, a.exp_type.name)
    end
    P = eltype(getfield.(a.runs, :pellet))
    if P isa Any
        withpellet = Final[]
        withoutpellet = Final[]
        for r in a.runs
            ismissing(r.pellet) ? push!(withoutpellet, r) : push!(withpellet, r)
        end
        for ds in [a.runs, withpellet, withoutpellet]
            n = length(ds)
            n > 1 && plot_density(E, ds, a.exp_type.name, n)
        end
    else
        ds = a.runs
        n = length(ds)
        n > 1 && plot_density(E, ds, a.exp_type.name, n)
    end
end

