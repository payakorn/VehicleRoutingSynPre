using Pkg

# # Pkg.add("~/.julia/dev/VRPTW.jl")
# push!(Base.load_path(), "/home/payakorn.s/.julia/dev/VehicleRoutingSynPre")
# Pkg.add(path="/home/payakorn.s/.julia/dev/VehicleRoutingSynPre")
Pkg.activate(".")
Pkg.instantiate()
# Pkg.resolve()

using VehicleRoutingSynPre

VehicleRoutingSynPre.find_opt(["ins25-$i" for i in 1:10]...)