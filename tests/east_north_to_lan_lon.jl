include("../core/utilities.jl")

using .Utilities
using ArgParse

# set cli arguments

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--easting", "-e"
        "--northing", "-n"
    end
    return parse_args(s)
end


args = parse_commandline()

function main()
    easting = try
        parse(Float64, args["easting"])
    catch
        parse(Int, args["easting"])
    end
    northing = try
        parse(Float64, args["northing"])
    catch
        parse(Int, args["northing"])
    end
    latlon = convert_easting_northing_to_latlon(easting, northing)
    println("lat: $(latlon.lat), lon: $(latlon.lon)")
end

main()
