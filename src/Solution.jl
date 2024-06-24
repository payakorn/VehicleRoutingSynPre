struct Problem
    name::String
    num_node::Int64
    num_vehi::Int64
    num_serv::Int64
    mind::Vector{Float64}
    maxd::Vector{Float64}
    a::Array{Int64}
    r::Array{Int64}
    d::Matrix{Float64}
    p::Array{Float64}
    e::Vector{Int64}
    l::Vector{Int64}
    PRE::Vector{Tuple}
    SYN::Vector{Tuple}
end


mutable struct IndexSolution
    vehicle::Int64
    position::Int64
    problem::Problem
end


function load_problem(Name::String)
    num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l = load_data(Name)

    # find Synchronization and Precedence constraints
    serv_r = find_service_request(r)
    # serv_a = find_compat_vehicle_node(a, r)
    SYN = find_SYN(serv_r, mind, maxd)
    PRE = find_PRE(serv_r, mind, maxd)

    return Problem(Name, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
end


function print_sol(sol::IndexSolution)
    println("vehicle: $(sol.vehicle)\nposition: $(sol.position)")
end


function find_starttime(sol::IndexSolution, route::Array)
    node = route[sol.vehicle][sol.position][1]
    serv = route[sol.vehicle][sol.position][2]
    println("node: $node\nserv: $serv\ne[$node]: $(sol.problem.e[node])")
    if sol.position == 1
        if sol.problem.d[1, node] < sol.problem.e[node]
            st = sol.problem.e[node]
        else
            st = sol.problem.d[1, node]
        end
    else
        sol.position -= 1
        st = find_starttime(sol, route)
    end
    return st
end
