using MATLAB

check(a::Basic{ClosedNest}) = checkbasic(a.path, ["nest", "feeder"])
check(a::Basic{TransferNest}) = checkbasic(a.path, ["south", "north", "feeder", "originalnest"])
check(a::Basic{Transfer}) = checkbasic(a.path, ["south", "north", "feeder"])
check(a::Basic{TransferNestBelen}) = checkbasic(a.path, ["southbefore", "northbefore", "feederbefore", "nestbefore", "north", "south"])
check(a::Basic{DawaySandpaper}) = checkbasic(a.path, ["nest", "feeder", "rightdowninitial", "leftdowninitial", "rightupinitial", "leftupinitial", "rightdownfinal", "leftdownfinal", "rightupfinal", "leftupfinal"])
check(a::Basic{Daway}) = checkbasic(a.path, ["nest", "feeder", "initialfeeder"])


check!(x, a::MATcalib) = _check!(x, a.file, ["cameraParams"])
check!(x, a::Calib) = _check!(x, a.file, ["fc", "cc", "alpha_c", "kc", "Tc_ext", "Rc_ext"])
check!(x, a::Calibration) = nothing

function _check!(x, filename::String, fields::Vector{String})
    file = MatFile(filename, "r")
    # println(filename)
    try 
        variables = variable_names(file)
        for field in fields
            x["no $field in the calibration file"] = field ∉ variables
        end
        close(file)
    catch 
        x["unknown error in the calibration file"] = true
    end
end

function checkbasic(path::String, fields::Vector{String})
    push!(fields, "track")
    @assert isdir(path) "Path does not exists"
    x = Dict{String, Bool}()
    x["no calibrations file"] = !isfile(joinpath(path, "calibrations.csv"))
    x["no factors file"] = !isfile(joinpath(path, "factors.csv"))
    x["no columns file"] = !isfile(joinpath(path, "columns.csv"))
    if !x["no calibrations file"]
        res2calib = CSV.read(joinpath(path, "calibrations.csv"), header=["resfile", "calibfile"], types=[String, String], transforms = Dict(1 => strip, 2 => strip))
        res2poi = CSV.read(joinpath(path, "columns.csv"), header=["resfile", "columnnumber", "poi"], types=[String, Int, String], transforms = Dict(3 => strip, 1 => strip))
        todo = join(res2calib, res2poi, on = :resfile)
        for f in fields
            x["no $f in columns file"] = f ∉ todo[:poi]
            x["no column for $f in columns file"] = ismissing(todo[todo[:poi] .== f, :columnnumber])
            x["no resfile for $f in columns file"] = ismissing(todo[todo[:poi] .== f, :resfile])
            x["no calibfile for $f in columns file"] = ismissing(todo[todo[:poi] .== f, :calibfile])
        end
        resfiles = unique(todo[:resfile])
        for file in resfiles
            x["no $file"] = !isfile(joinpath(path, file))
        end
        calibfiles = unique(todo[:calibfile])
        for file in calibfiles
            x["no $file"] = !isfile(joinpath(path, file))
            calib = get_calibration(path, file)
            check!(x, calib)
        end
        for r in eachrow(todo)
            if r[:poi] ∈ [fields; "pellet"]
                file = MatFile(joinpath(path, r[:resfile]), "r")
                X = get_variable(file, "xdata")
                XI = nonzeros(X[:, r[:columnnumber]])
                x["empty column $(r[:columnnumber]) for $(r[:poi]) in $(r[:resfile])"] = isempty(XI)
                close(file)
            end
        end
    end
    if !x["no factors file"]
        factors = CSV.read(joinpath(path, factor_filename), header=["Factor", "Level"], types=[String, Union{String, Missing}], transforms = Dict(1 => strip))
        i = findfirst(factors[:Factor], "check")
        if i ≠ 0
            x["check corrupt in factors file"] = !r"^(?:\d+(\.\d+)?)\s(?:mm|cm|m)$"(strip(factors[i, :Level]))
        end
        i = findfirst(factors[:Factor], "azimuth")
        if i ≠ 0
            x["azimuth corrupt in factors file"] = !r"^(?:\d+(\.\d+)?)\s(?:rad|°)$"(strip(factors[i, :Level]))
        end
        i = findfirst(factors[:Factor], "nest2feeder")
        if i ≠ 0
            x["nest2feeder corrupt in factors file"] = !r"^(?:\d+(\.\d+)?)\s(?:mm|cm|m)$"(strip(factors[i, :Level]))
        end
    end
    [k for (k,v) in x if v]
end


