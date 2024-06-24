mutable struct Particle
    route::Vector{Vector{Vector{Int64}}}
    starttime::Dict{Int64, Array{Float64}}
    slot::Dict{Int64, Vector{Int64}}
    serv_a::Tuple
    serv_r::Dict{Int64, Vector{Int64}}
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


"""
    load_data(name::String) 

    return number of node (including depot), number of vehicle, number of services, min precedence, max, precedence, matrix 'a', matrix 'r', distance, processing time, earliest, latest duedate

    Example:

    num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l = load_data("ins10-1")
"""
function load_data(name::String)
    f = JLD2.load(joinpath(@__DIR__, "..", "data", "raw_HHCRSP", "$name.jld2"))
    return f["num_node"], f["num_vehi"], f["num_serv"], f["mind"], f["maxd"], f["a"], f["r"], f["d"], f["p"], f["e"], f["l"]
end


"""
    find_group_of_node(num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)

    devide node to compatible vehicle randomly

"""
function find_group_of_node(num_node::Int64, num_vehi::Int64, a::Matrix, r::Matrix)
    node_vehicle = r*a'
    group_node = [Int64[] for i in 1:num_vehi]
    # group_node = Dict(i => Int64[1] for i in 1:num_vehi)
    for i in 2:num_node
        if sum(r[i, :]) == 1
            push!(group_node[rand(findall(x->x==1, node_vehicle[i, :]))], i)
        else
            jobs = findall(x -> x > 0, node_vehicle[i, :])
            for j in jobs
                push!(group_node[j], i)
            end
        end
    end
    return group_node
end


function find_vehicle_by_service(route::Array, node::Int64, service::Int64, num_vehi::Int64)
    for vehi in 1:num_vehi
        for i in route[vehi]
            if issubset((node, service), i)
                return vehi
            end
        end
    end
    return nothing
end


function in_SYN(node::Int64, SYN::Vector)
    # any(in.(node, SYN))
    any(x[1] == node for x in SYN)
end


function in_PRE(node::Int64, PRE::Vector)
    # any(in.(node, SYN))
    any(x[1] == node for x in PRE)
end


function union_route(v::Vector)
    u = v[1]
    for i in 2:length(v)
        union!(u, v[i])
    end
    return u
end


function update_starttime(input_route::Array, slot::Dict{Int64, Vector{Int64}}, node::Int64, service::Int64, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array, st::Dict)
    vehi, loca = find_location_by_node_service(input_route, node, service)
    for lo in loca+1:length(input_route[vehi])
        cnode = input_route[vehi][lo][1]
        cserv = input_route[vehi][lo][2]
        cvehi, cloca = find_location_by_node_service(input_route, cnode, cserv)

        if in_PRE(cnode, PRE)
            st = find_node_starttime(input_route, slot, cnode, cserv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st, vehicle=vehi)
        elseif in_SYN(cnode, SYN)
            st = find_node_starttime(input_route, slot, cnode, cserv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st, vehicle=vehi)
        else
            if st[node][vehi, service] + p[vehi, service, node] + d[node, cnode] < e[cnode]
                st[cnode][cvehi, cserv] = e[cnode]
            else
                st[cnode][cvehi, cserv] = st[node][vehi, service] + p[vehi, service, node] + d[node, cnode]
            end
        end

    end
    return st
end


function find_node_starttime(input_route::Array, slot::Dict{Int64, Vector{Int64}}, node::Int64, service::Int64, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array, st::Dict; vehicle=1, oneback=false)
    if node == 1 && service == 1
        vehi, loca = (vehicle, length(input_route[vehicle]))
    else
        vehi, loca = find_location_by_node_service(input_route, node, service)
    end


    # println("finding starting time of node: $node, service: $service, at vehicle: $vehi in the No. $loca")

    # check if the depot [1, 1]
    if loca != 1
        left_node = input_route[vehi][loca-1][1]
        left_serv = input_route[vehi][loca-1][2]
        # st = find_node_starttime(input_route, slot, left_node, left_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st, vehicle=vehi)
        st_left = st[left_node][vehi, left_serv] + p[vehi, service, left_node] + d[left_node, node]
    else
        st_left = d[1, node]
    end

    # calculate completion time
    if st_left < e[node]
        st[node][vehi, service] = e[node]
    else
        st[node][vehi, service] = st_left
    end

    if in_SYN(node, SYN)
        syn = SYN[findfirst(x->x[1]==node, SYN)]
        other_serv = setdiff(syn, [node, service])[1]
        ovehi, oloca = find_location_by_node_service(input_route, node, other_serv)

        if vehi == ovehi
            return ErrorException("infeasible")
        end
        # update starttime of syn node
        # if st[node][ovehi, other_serv] != st[node][vehi, service]
        #     st = find_node_starttime (input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
        # end

        if st[node][ovehi, other_serv] != 0.0
            if st[node][ovehi, other_serv] < st[node][vehi, service]
                st[node][ovehi, other_serv] = st[node][vehi, service]

                next_node = input_route[ovehi][oloca+1][1]
                next_serv = input_route[ovehi][oloca+1][2]
                next_st = st[next_node][ovehi, next_serv]
                if next_st > e[next_node] || st[node][ovehi, other_serv] + p[ovehi, other_serv, node] + d[node, next_node] > e[next_node]
                    # println("update starttime of syn $node $service to $node $other_serv")
                    st = update_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
                end
            else
                st[node][vehi, service] = st[node][ovehi, other_serv]
            end
        end
    end
    
    if in_PRE(node, PRE)
        pre = PRE[findfirst(x->x[1]==node, PRE)]
        other_serv = setdiff(pre, [node, service])[1]
        if other_serv in slot[node]
            in_position = findfirst(x->x==service, pre) # current position of service
            ovehi, oloca = find_location_by_node_service(input_route, node, other_serv)

            # update starttime of PRE node
            # st = find_node_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)

            if st[node][ovehi, other_serv] != 0.0
                if in_position == 2
                    if (st[node][ovehi, other_serv] - st[node][vehi, service]) < mind[node]
                        st[node][ovehi, other_serv] = mind[node] + st[node][vehi, service]

                        # same car
                        if vehi == ovehi
                            return st
                        end

                        next_node = input_route[ovehi][oloca+1][1]
                        next_serv = input_route[ovehi][oloca+1][2]
                        next_st = st[next_node][ovehi, next_serv]
                        if next_st > e[next_node] || st[node][ovehi, other_serv] + p[ovehi, other_serv, node] + d[node, next_node] > e[next_node]
                            st = update_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
                        end

                        # st = update_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
                    elseif (st[node][ovehi, other_serv] - st[node][vehi, service]) > maxd[node]
                        st[node][vehi, service] = st[node][ovehi, other_serv] - maxd[node]
                    end
                elseif in_position == 3
                    if (st[node][vehi, service] - st[node][ovehi, other_serv]) < mind[node]
                        st[node][vehi, service] = mind[node] + st[node][ovehi, other_serv]
                    elseif (st[node][vehi, service] - st[node][ovehi, other_serv]) > maxd[node]
                        st[node][ovehi, other_serv] = st[node][vehi, service] - maxd[node]

                        # same car
                        if vehi == ovehi
                            return st
                        end

                        next_node = input_route[ovehi][oloca+1][1]
                        next_serv = input_route[ovehi][oloca+1][2]
                        next_st = st[next_node][ovehi, next_serv]
                        if next_st > e[next_node] || st[node][ovehi, other_serv] + p[ovehi, other_serv, node] + d[node, next_node] > e[next_node]
                            st = update_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
                        end

                        # st = update_starttime(input_route, slot, node, other_serv, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st)
                    end
                end
            end
        end
    end
    return st
end


function find_node_starttime(par::Particle, node::Int64, service::Int64, st::Dict; vehicle=1)
    find_node_starttime(par.route, par.slot, node, service, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.PRE, par.SYN, st, vehicle=vehicle)
end


# function find_starttime2(input_route::Array, slot::Dict{Int64, Vector{Int64}}, node::Int64, vehi::Int64, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)
#     notcalculate = Vector{Int64}[]
#     route = input_route[vehi]
#     # for inode in route[node:end]
#     #     @show inode
#     # end
# end


function find_starttime2(par::Particle, node::Int64, vehi::Int64)
    return find_starttime2(par.route, par.slot, node, vehi, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)
end


function find_starttime(route::Array, slot::Dict, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Vector, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Vector, l::Vector, PRE::Vector, SYN::Vector)
    st = initial_starttime(num_node, num_vehi, num_serv)
    # for i in 1:num_vehi
    #     st = find_node_starttime(route, slot, 1, 1, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st, vehicle=i)
    # end
    # return st

    for (i, vehi) in enumerate(route)
        for ns in vehi
            # println("finding starttime of $ns")
            st = find_node_starttime(route, slot, ns[1], ns[2], num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN, st, vehicle=i)
        end
    end
    return st
end


function find_starttime(par::Particle)
    # st = initial_starttime(par.num_node, par.num_vehi, par.num_serv)
    # for i in 1:par.num_vehi
    #     st = find_node_starttime(par, 1, 1, st, vehicle=i)
    # end
    return find_starttime(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)

    # for (i, vehi) in enumerate(par.route)
    #     for ns in vehi
    #         st = find_node_starttime(par, ns[1], ns[2], st, vehicle=i)
    #     end
    # end
    # return st
end


function initial_starttime(num_node::Int64, num_vehi::Int64, num_serv::Int64)
    Dict(i => zeros(Float64, num_vehi, num_serv) for i in 0:num_node)
end


# function find_starttime(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)
#     st = Dict(i => zeros(Float64, num_vehi, num_serv) for i in 0:num_node)

#     for vehi in 1:num_vehi
#         for i in 1:length(route[vehi])
#             if i == 1 # start at node 1 (origin)
#                 left = 1
#                 right = route[vehi][i][1]

#                 if d[left, right] > e[route[vehi][i][1]]
#                     st[right][vehi, route[vehi][i][2]] = d[left, right]
#                 else
#                     st[right][vehi, route[vehi][i][2]] = e[route[vehi][i][1]]
#                 end

#                 # if right != 1
#                 #     for sl in 2:length(slot[right])
#                 #         v = find_vehicle_by_service(route, right, slot[right][sl-1], num_vehi)
#                 #         st[right][v, slot[right][sl]] = maximum(st[right][slot[right][sl-1]])
#                 #     end
#                 # end
#             else
#                 left = route[vehi][i-1][1]
#                 right = route[vehi][i][1]

#                 if st[left][vehi, route[vehi][i-1][2]] + p[vehi, route[vehi][i-1][2], left] + d[left, right] > e[route[vehi][i][1]]
#                     st[right][vehi, route[vehi][i][2]] = st[left][vehi, route[vehi][i-1][2]] + p[vehi, route[vehi][i-1][2], left] + d[left, right]
#                 else
#                     st[right][vehi, route[vehi][i][2]] = e[route[vehi][i][1]]
#                 end
#                 # println("vehi: $vehi st[$right][$vehi, route[$vehi][$i][2]] = $(st[right][vehi, route[vehi][i][2]])" )
#             end

#             if right != 1
#                 for sl in 2:length(slot[right])
#                     v = find_vehicle_by_service(route, right, slot[right][sl-1], num_vehi)
#                     vv = find_vehicle_by_service(route, right, slot[right][sl], num_vehi)
#                     if any([(right, slot[right][sl-1], slot[right][sl]) == i || (right, slot[right][sl], slot[right][sl-1]) == i for i in SYN]) # synchronization
#                         if maximum(st[right][:, slot[right][sl-1]]) > st[right][vv, slot[right][sl]]
#                             st[right][vv, slot[right][sl]] = maximum(st[right][:, slot[right][sl-1]])
#                         elseif maximum(st[right][:, slot[right][sl-1]]) < st[right][vv, slot[right][sl]]
#                             diff_starttime = st[right][vv, slot[right][sl]] - st[right][v, slot[right][sl-1]]
#                             # st[right][v, slot[right][sl]] = st[right][v, slot[right][sl]]
#                             # vv = find_vehicle_by_service(route, right, slot[right][sl-1], num_vehi)
#                             v, loca = find_location_by_node_service(route, right, slot[right][sl-1])
#                             st[right][vv, slot[right][sl-1]] = st[right][v, slot[right][sl]]
#                             for k in loca+1:length(route[vv])
#                                 st[route[vv][k][1]][vv, route[vv][k][2]] += diff_starttime
#                             end
#                         end
#                     elseif any([(right, slot[right][sl-1], slot[right][sl]) == i for i in PRE]) # precedence
#                         # @show right
#                         # @show v
#                         # before_vehi = argmax(st[right][:, slot[right][sl-1]])
#                         # starttime_before_vehi = st[right][:, slot[right][sl-1]][vehi]
#                         # v is vehicle of service sl-1
#                         # vv is vehicle of service sl
#                         if st[right][vv, slot[right][sl]] - st[right][v, slot[right][sl-1]] < mind[right]
#                             if v != vv # if they come from different vcles
#                                 st[right][vv, slot[right][sl]] += mind[right] - (st[right][vv, slot[right][sl]] - st[right][v, slot[right][sl-1]])
#                             else
#                                 if st[right][vv, slot[right][sl]] - (st[right][v, slot[right][sl-1]] + p[v, slot[right][sl-1], right]) < mind[right]
#                                     st[right][vv, slot[right][sl]] += st[right][vv, slot[right][sl]] - (st[right][v, slot[right][sl-1]] + p[v, slot[right][sl-1], right])
#                                 end
#                             end
#                         end
#                     else
#                         st[right][vv, slot[right][sl]] = maximum([maximum(st[right][:, slot[right][sl-1]]) + p[vv, slot[right][sl-1], right], st[right][vv, slot[right][sl]]])
#                     end
#                 end
#             end
#             # println("vehicle: $vehi, depart: $left, arrive: $right, start $(st[right][vehi, route[vehi][i][2]])")
#         end
#     end
#     return st
# end


# function find_starttime(par::Particle)
#     return find_starttime(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)
# end


function create_empty_slot(num_node::Int64)
    return Dict(i => Int64[] for i in 2:num_node)
end


function find_service_request(r::Array)
    # serv_requst = [findall(x->x==1, r[2, :])]
    # for node in 3:size(r, 1)-1
    #     append!(serv_requst, [findall(x->x==1, r[node, :])])
    # end
    # return serv_requst
    return Dict(node => findall(x->x==1, r[node, :]) for node in 2:size(r, 1)-1)
end


function find_service(par::Particle)
    return find_service(par.route, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)
end


function generate_empty_particle(serv_a::Tuple, serv_r::Dict, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)

    # cluter nodes to vehicle
    group_node = find_group_of_node(num_node, num_vehi, a, r)

    # create empty route
    route = [[Int64[1, 1]] for _ in 1:num_vehi]
    # route = [[[1, 1]], [[9, 5], [1, 1]], [[2, 4], [9, 6], [1, 1]]] # for test

    # create empty slot
    slot = create_empty_slot(num_node)
    # slot[9] = [5, 6]

    # just add
    # starttime = find_starttime(route, num_node, slot, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, SYN, PRE)
    starttime = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
    # service = create_empty_slot(num_node)

    return Particle(route, starttime, slot, serv_a, serv_r, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
end


function initial_insert_service(serv_r::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, SYN::Vector{Tuple})
    slot = deepcopy(serv_r)
    for serv in slot
        slot[serv[1]] = [rand(serv[2])]
    end

    # SYN
    for syn in SYN
        slot[syn[1]] = collect(syn[2:end])
    end
    return slot
end


function check_assigned_node(route::Array, node::Int64, service::Int64, num_vehi::Int64)
    for vehi in 1:num_vehi
        if [node, service] in route[vehi]
            return true
        end
    end
    return false
end


function choose_vehicle_from_service(a::Matrix{Int64}, service::Vector{Int64})
    sort_num_service = sortperm([length(findall(x->x==1, a[:, i])) for i in service])
    syn_dict = Dict()
    # avai_vehi = [i for i in 1:size(a, 1)]
    forbid = []
    for s in service[sort_num_service]
        xx = setdiff(findall(a[:, s] .== 1), forbid)
        syn_dict[s] = rand(xx)
        push!(forbid, syn_dict[s])
    end
    return (syn_dict[s] for s in service)
end


function insert_vehicle_to_service(input_route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Vector, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Vector, l::Vector, SYN::Vector, PRE::Vector)

    route = deepcopy(input_route)
    st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)

    for syn in SYN
        for i in 2:length(syn)-1
            if !check_assigned_node(route, syn[1], syn[i], num_vehi) && !check_assigned_node(route, syn[1], syn[i+1], num_vehi) && issubset([syn[i], syn[i+1]], slot[syn[1]])
                vehi1, vehi2 = choose_vehicle_from_service(a, [syn[i], syn[i+1]])
                insert!(route[vehi1], 1, [syn[1], syn[i]])
                insert!(route[vehi2], 1, [syn[1], syn[i+1]])
            end
        end
    end
    
    for node in setdiff(sortperm(e[1:end-1])[2:end], [sl[1] for sl in SYN])
        # println("$node, $(slot[node])")
        for sl in slot[node]
            if !check_assigned_node(route, node, sl, num_vehi)
                vehi = rand(findall(x->x==1, a[:, sl]))
                insert!(route[vehi], length(route[vehi]), [node, sl])
            end
        end
    end
    return route
end


function insert_vehicle_to_service(par::Particle)
    par.route = insert_vehicle_to_service(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.SYN, par.PRE)
    par.starttime = find_starttime(par)
    return par
end


function find_SYN(serv_request::Dict, mind::Vector, maxd::Vector)
    SYN = Tuple[]
    if length(mind) <= 75
        runnum = length(mind)-1
    else
        runnum = length(mind)
    end
    
    for i in 2:runnum
        if mind[i] == 0 && maxd[i] == 0 && length(serv_request[i]) == 2
            push!(SYN, (i, serv_request[i]...))
        end
    end
    return SYN
end

"""
    function find_PRE(serv_request::Dict)

        (node, service1, service2) => service1 must start before service2 at the node 

"""
function find_PRE(serv_request::Dict, mind::Array, maxd::Array)
    PRE = Tuple[]
    for i in 2:length(serv_request)+1
        if length(serv_request[i]) > 1
            if mind[i] != 0 && maxd[i] != 0
                for j in 1:length(serv_request[i])-1
                    push!(PRE, (i, serv_request[i][j], serv_request[i][j+1]))
                end
            end
        end
    end
    return PRE
end


function find_vehicle_service(a::Matrix)
    vehicle_compat = collect(Tuple(findall(x->x==1, a[vehi, :])) for vehi in 1:size(a, 1))
    return vehicle_compat
end


function find_compat_vehicle_node(a::Array, r::Array)
    compat_matrix = a*r'
    return Tuple(Tuple(findall(x->x>0, compat_matrix[i, :])) for i in 1:size(compat_matrix, 1))
end


function generate_particles(name::String)
    # load data
    num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l = load_data(name)

    # find Synchronization and Precedence constraints
    serv_r = find_service_request(r)
    serv_a = find_compat_vehicle_node(a, r)
    SYN = find_SYN(serv_r, mind, maxd)
    PRE = find_PRE(serv_r, mind, maxd)

    par = generate_empty_particle(serv_a, serv_r, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
    # par.slot = initial_insert_service(serv_r, num_node, num_vehi, num_serv, SYN) # also insert synchronization
    par.slot = serv_r
    par = insert_vehicle_to_service(par)
    # par = insert_PRE(par)
end


function find_remain_service(slot::Dict, serv_r::Dict)
    for (i, sl) in slot
        try setdiff!(serv_r[i], slot[i]) catch e; continue end
        if isempty(serv_r[i])
            delete!(serv_r, i)
        end
    end
    return serv_r
end


function find_remain_service(particle::Particle)
    return find_remain_service(particle.slot, particle.serv_r)
end


function find_location_by_node_service(route::Array, node::Int64, service::Int64)
    for (num_v, vehi) in enumerate(route)
        for (num_loca, n) in enumerate(vehi)
            if [node, service] == [n[1], n[2]]
                return (num_v, num_loca)
            end
        end
    end
    return (nothing, nothing)
end


function insert_node(input_route::Array, slot::Dict{Int64, Vector{Int64}}, node, service, num_node::Int64, num_vehi::Int64, serv_a::Tuple, num_serv::Int64, mind::Vector, maxd::Vector, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Vector, l::Vector, serv_r::Dict{Int64, Vector{Int64}}, PRE::Vector{Tuple}, SYN::Vector)
    route = deepcopy(input_route)
    for vehi in route
        for numnode in 1:length(vehi)
            insert!(vehi, numnode, [node, service])
            st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
            if compatibility(route, slot, a, r, serv_a, serv_r, PRE, SYN)
                nothing
            end
        end
    end
end

function check_all_insert_node(particle::Particle)
    return check_all_insert_node(particle.route, particle.r, particle.e)
end



function check_all_insert_node(route::Array, r::Matrix, e::Vector)

    all_route = []
    for vehi in route
        append!(all_route, vehi)
    end

    # all_node = getindex.(findall(BitMatrix(r)), [1 2])
    all_node = findall(BitMatrix(r))
    for x in all_node
        if [x[1], x[2]] in all_route
            continue
        elseif x[1] == 1 || x[2] == 1 || x[1] == length(e) || x[2] == length(e)
            continue
        else
            println("no node $(x[1]) service $(x[2])")
            return false
        end
    end
    return true
end
            

function insert_PRE(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Vector, a::Matrix, r::Matrix, d::Matrix, p::Array, e::Vector, l::Vector, serv_r::Dict{Int64, Vector{Int64}}, PRE::Vector{Tuple}, SYN::Vector)
    serv_r = find_remain_service(slot, serv_r)
    for (node, remain_service) in serv_r
        for serv in remain_service
            # check Precedence
            for pre_service in findall(issubset((node, serv), x) for x in PRE)
                location = findfirst(x->x==serv, PRE[pre_service])
                if location == 2
                    location_of_slot = findfirst(x->x==PRE[pre_service][3], slot[node])
                    (v, loca) = find_location_by_node_service(route, node, slot[node][location_of_slot])
                    for vv in findall(x->x==1, a[:, serv])
                        loca_run = 1
                        while check_assigned_node(route, node, PRE[pre_service][2], num_vehi) == false && loca_run <= length(route[vv])
                            if a[v, serv] == 1
                                insert!(route[v], loca, [node, PRE[pre_service][2]])
                                insert!(slot[node], location_of_slot, PRE[pre_service][2])
                            else
                                insert!(route[vv], loca_run, [node, PRE[pre_service][2]])
                                insert!(slot[node], location_of_slot, PRE[pre_service][2])
                            end
                            st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
                            if !check_PRE(route, st, maxd, PRE)
                                if a[v, serv] == 1
                                    deleteat!(route[vv], loca)
                                    deleteat!(slot[node], location_of_slot)
                                else
                                    deleteat!(route[vv], loca_run)
                                    deleteat!(slot[node], location_of_slot)
                                end
                            end
                            loca_run += 1
                        end
                    end
                else
                    location_of_slot = findfirst(x->x==PRE[pre_service][2], slot[node])
                    (v, loca) = find_location_by_node_service(route, node, slot[node][location_of_slot])
                    for vv in findall(x->x==1, a[:, serv])
                        loca_run = 1
                        while check_assigned_node(route, node, PRE[pre_service][3], num_vehi) == false && loca_run <= length(route[vv])-1
                            if a[v, serv] == 1
                                insert!(route[v], loca+1, [node, PRE[pre_service][3]])
                                insert!(slot[node], location_of_slot+1, PRE[pre_service][3])
                            else
                                # vv = rand(findall(x->x==1, a[:, serv]))
                                insert!(route[vv], loca_run, [node, PRE[pre_service][3]])
                                insert!(slot[node], location_of_slot+1, PRE[pre_service][3])
                            end
                            st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
                            if !check_PRE(route, st, maxd, PRE)
                                if a[v, serv] == 1
                                    deleteat!(route[vv], loca+1)
                                    deleteat!(slot[node], location_of_slot+1)
                                else
                                    deleteat!(route[vv], loca_run)
                                    deleteat!(slot[node], location_of_slot+1)
                                end
                            end
                            loca_run += 1
                        end
                    end
                end
            end
        end
    end
    return route, slot, serv_r
end


function insert_PRE(par::Particle)
    route, slot, serv_r = insert_PRE(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.d, par.p, par.e, par.l, par.serv_r, par.PRE, par.SYN)
    par.route = route
    par.slot = slot
    par.serv_r = serv_r
    par.starttime = find_starttime(par)
    return par
end


function insert_service_to_node_slot(particle::Particle, node::Int64, slot::Int64)
    nothing
end


function compatibility(particle::Particle)
    return compatibility(particle.route, particle.slot, particle.a, particle.r, particle.serv_a, particle.serv_r, particle.PRE, particle.SYN)
end


function check_PRE(route::Array, starttime::Dict, maxd::Vector, PRE::Vector{Tuple})
    for (node, j, k) in PRE
        (v1, lo1) = find_location_by_node_service(route, node, j)
        (v2, lo2) = find_location_by_node_service(route, node, k)

        if isnothing(v1) || isnothing(v2) || isnothing(lo1) || isnothing(lo2)
            continue
        elseif v1 == v2
            if lo1 > lo2 || (starttime[node][v2, k] - starttime[node][v1, j]) > maxd[node]
                return false
            end
        else
            if (starttime[node][v2, k] - starttime[node][v1, j]) > maxd[node]
                return false
            end
        end
    end
    return true
end


function check_PRE(particle::Particle)
    check_PRE(particle.route, particle.starttime, particle.maxd, particle.PRE)
end


function check_SYN(route::Array, SYN::Vector{Tuple})
    for (node, j, k) in SYN
        (v1, _) = find_location_by_node_service(route, node, j)
        (v2, _) = find_location_by_node_service(route, node, k)
        if v1 == v2
            return false
        end
    end
    return true
end


function check_SYN(particle::Particle)
    check_SYN(particle.route, particle.SYN)
end


function check_all_syn_pre(particle::Particle)
    println("ALL: $(check_all_insert_node(particle))")
    println("COM: $(compatibility(particle))")
    println("PRE: $(check_PRE(particle))")
    println("SYN: $(check_SYN(particle))")
end


function compatibility(route::Array, slot::Dict, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, PRE::Vector{Tuple}, SYN::Vector{Tuple})

    for (numvehi, vehi) in enumerate(route)
        for node_service in vehi[1:end-1]
            # println("node $(node_service[1]) service $(node_service[2])")
            if r[node_service[1], node_service[2]] != 1 || a[numvehi, node_service[2]] != 1
                return false
            end
        end
    end
    return true
end


function feasibility(particle::Particle)
    if !compatibility(particle)
        return false
    else
        return true
    end
end


function insert_initial_SYN(particle::Particle)
    nothing
end


function total_distance(route::Array, d::Matrix)
    dis = 0.0
    route1 = deepcopy(route)
    for vehi in route1
        insert!(vehi, 1, [1, 1])
        for node in 1:length(vehi)-1
            dis += d[vehi[node][1], vehi[node+1][1]]
        end
    end
    return dis
end


function total_distance(particle::Particle)
    return total_distance(particle.route, particle.d)
end


function tardiness(route::Array, starttime::Dict, p::Array, l::Vector)
    tardy = 0.0
    max_tardy = 0.0
    for (i, vehi) in enumerate(route)
        for node_service in vehi
            st = starttime[node_service[1]][i, node_service[2]] - l[node_service[1]]
            if st > 0.0
                tardy += st
                if st > max_tardy
                    max_tardy = deepcopy(st)
                end
            end
        end
    end 
    return tardy, max_tardy
end


function tardiness(particle::Particle)
    return tardiness(particle.route, particle.starttime, particle.p, particle.l)
end


function objective_value(route, starttime, p, l, d; lambda=[1/3, 1/3, 1/3])
    tardy, max_tardy = tardiness(route, starttime, p, l)
    return lambda[1]*total_distance(route, d) + lambda[2]*tardy + lambda[3]*max_tardy
end


function objective_value(particles::Vector)
    obj = Float64[]
    for par in particles
        push!(obj, objective_value(par))
    end
    return obj
end


function objective_value(particle::Particle; lambda=[1/3, 1/3, 1/3])
    tardy, max_tardy = tardiness(particle)
    return lambda[1]*total_distance(particle) + lambda[2]*tardy + lambda[3]*max_tardy
end


function example()
    route = [
        [[11, 3], [4, 2], [6, 3], [10, 1], [8, 3], [1, 1]],
        [[9, 6], [1, 1]],
        [[9, 5], [11, 6], [7, 5], [3, 5], [2, 4], [10, 4], [5, 4], [1, 1]]
    ]
    
    slot = Dict(i => Int64[] for i in 2:11)
    slot[2] = [4]
    slot[3] = [5]
    slot[4] = [2]
    slot[5] = [4]
    slot[6] = [3]
    slot[7] = [5]
    slot[8] = [3]
    slot[9] = [5, 6]
    slot[10] = [1, 4]
    slot[11] = [3, 6]

    return route, slot
end


function example2()
    route = [
        [[10, 1], [3, 1], [4, 3], [8, 1], [11, 3], [5, 2], [7, 1], [1, 1]],
        [[2, 5], [9, 4], [11, 5], [1, 1]],
        [[10, 6], [9, 5], [6, 4], [1, 1]]
    ]
    
    slot = Dict(i => Int64[] for i in 2:11)
    slot[2] = [5]
    slot[3] = [1]
    slot[4] = [3]
    slot[5] = [2]
    slot[6] = [4]
    slot[7] = [1]
    slot[8] = [1]
    slot[9] = [4, 5]
    slot[10] = [1, 6]
    slot[11] = [3, 5]

    return route, slot
end


function generate_example()
    route, slot = example()
    par = generate_particles("ins10-1")
    par.route = route
    par.slot = slot
    par.starttime = find_starttime(par)
    return par
end


function complete(route, slot)
    nothing
end


function complete(particle::Particle)
    nothing
end


function List(num_node::Int64, num_vehi::Int64, num_serv::Int64, r::Matrix, SYN::Vector{Tuple}, PRE::Vector{Tuple})
    r[1, :] = zeros(1, num_serv)
    index = getindex.(findall(x->x==1, r[1:num_node, :]), [1 2])
    set_of_all_index = [Tuple(index[i, :]) for i in 1:size(index, 1)]
    # set = Iterators.filter(x -> x[1] != x[2], Iterators.product(set_of_all_index, set_of_all_index))
    set = combinations(set_of_all_index, 2)
    return set
end


function List(par::Particle)
    return List(par.num_node, par.num_vehi, par.num_serv, par.r, par.SYN, par.PRE)
end


function test_forloop(list)
    for (i, j) in enumerate(list)
        if i == 500
            println("reach 500 iterations")
        end
    end
end


function in_same_route(input_route::Array)
    for (vehi, route) in enumerate(input_route)
        all_node = [x[1] for x in route]
        old_len = length(all_node)
        unique!(all_node)
        new_len = length(all_node)
        if old_len > new_len
            # println("Same vehicle")
            return false
        end
    end
    return true
end


function in_same_route(par::Particle)
    return in_same_route(par.route)
end


function swap(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array, list)
    input_route = deepcopy(route)
    for ls in list
        num_v1, num_loca1 = find_location_by_node_service(input_route, ls[1][1], ls[1][2])
        num_v2, num_loca2 = find_location_by_node_service(input_route, ls[2][1], ls[2][2])
        test_route = deepcopy(input_route)
        test_route[num_v1][num_loca1], test_route[num_v2][num_loca2] = test_route[num_v2][num_loca2], test_route[num_v1][num_loca1]
        if compatibility(test_route, slot, a, r, serv_a, serv_r, PRE, SYN) && in_same_route(test_route) && num_v1 != num_v2
            test_st = try find_starttime(test_route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN) catch StackOverflowError; continue end
            st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
            if objective_value(test_route, test_st, p, l, d) < objective_value(route, st, p, l, d) && check_PRE(test_route, test_st, maxd, PRE) && check_SYN(test_route, SYN)
                input_route = deepcopy(test_route)
            end
        end
    end
    return input_route
end


function swap(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)
    list = List(num_node, num_vehi, num_serv, r, SYN, PRE)
    return swap(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, serv_a, serv_r, d, p, e, l, PRE, SYN, list)
end


function swap(par::Particle; list=nothing)
    if isnothing(list)
        route = swap(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.serv_a, par.serv_r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)
    else
        route = swap(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.serv_a, par.serv_r, par.d, par.p, par.e, par.l, par.PRE, par.SYN, list)
    end
    par.route = route
    par.starttime = find_starttime(par)
    return par
end


function move(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array, list)
    input_route = deepcopy(route)
    for ls in list
        num_v1, num_loca1 = find_location_by_node_service(input_route, ls[1][1], ls[1][2])
        num_v2, num_loca2 = find_location_by_node_service(input_route, ls[2][1], ls[2][2])
        test_route = deepcopy(input_route)

        #  move
        moved_item = splice!(test_route[num_v1], num_loca1)
        insert!(test_route[num_v2], num_loca2, moved_item)

        # find starttime
        
        if compatibility(test_route, slot, a, r, serv_a, serv_r, PRE, SYN) && check_SYN(test_route, SYN)

            test_st = try find_starttime(test_route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN) catch StackOverflowError; continue end
            st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)

            if objective_value(test_route, test_st, p, l, d) < objective_value(route, st, p, l, d) && check_PRE(test_route, test_st, maxd, PRE)
                # println("cost reduce in Move")
                input_route = deepcopy(test_route)
            end
        end
    end
    return input_route
end


function random_move(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array, list)
    input_route = deepcopy(route)
    list = collect(list)
    for ls in list[1:Int(round(length(list)/4))]
        num_v1, num_loca1 = find_location_by_node_service(input_route, ls[1][1], ls[1][2])
        num_v2, num_loca2 = find_location_by_node_service(input_route, ls[2][1], ls[2][2])
        test_route = deepcopy(input_route)

        #  move
        moved_item = splice!(test_route[num_v1], num_loca1)
        insert!(test_route[num_v2], num_loca2, moved_item)

        # find starttime
        test_st = find_starttime(test_route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
        st = find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)

        if compatibility(test_route, slot, a, r, serv_a, serv_r, PRE, SYN)
            println("random Move")
            input_route = deepcopy(test_route)
        end
    end
    return input_route
end


function move(route::Array, slot::Dict{Int64, Vector{Int64}}, num_node::Int64, num_vehi::Int64, num_serv::Int64, mind::Vector, maxd::Array, a::Matrix, r::Matrix, serv_a::Tuple, serv_r::Dict, d::Matrix, p::Array, e::Array, l::Array, PRE::Array, SYN::Array)
    list = List(num_node, num_vehi, num_serv, r, SYN, PRE)
    return move(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, serv_a, serv_r, d, p, e, l, PRE, SYN, list)
end


function move(par::Particle; list=nothing)
    if isnothing(list)
        route = move(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.serv_a, par.serv_r, par.d, par.p, par.e, par.l, par.PRE, par.SYN)
    else
        route = move(par.route, par.slot, par.num_node, par.num_vehi, par.num_serv, par.mind, par.maxd, par.a, par.r, par.serv_a, par.serv_r, par.d, par.p, par.e, par.l, par.PRE, par.SYN, list)
    end
    par.route = route
    par.starttime = find_starttime(par)
    return par
end


function path_relinking(route, best_route, a)
    test_route = deepcopy(route)

    rand_vehi = rand(1:length(route))
    remove_route = test_route[rand_vehi]
    remove_best_route = best_route[rand_vehi]
    
    remain_node = setdiff(remove_route, remove_best_route)

    # insert best route
    test_route[rand_vehi] = best_route[rand_vehi]

    # insert
    for ns in remain_node
        choose_vehi = rand(findall(x->x==1, a[:, ns[2]]))
        choose_posi = length(test_route[choose_vehi])
        insert!(test_route[choose_vehi], choose_posi, ns)
    end
    return test_route
end


function path_relinking(par::Particle, par_best::Particle)
    # best_route = rand(par_best.route)
    route = path_relinking(par.route, par_best.route, par.a)
    par.route = route
    par.starttime = find_starttime(par)
    return par
end


function local_search(particle::Particle, best_par::Particle)
    test_par = deepcopy(particle)
    list = List(particle.num_node, particle.num_vehi, particle.num_serv, particle.r, particle.SYN, particle.PRE)
    # first_obj = objective_value(test_par)
    test_par = swap(test_par, list=list)
    test_par = move(test_par, list=list)
    # test_par = path_relinking(test_par, best_par)
    return test_par
end


function local_search_ver2(particle::Particle)
    list = List(particle.num_node, particle.num_vehi, particle.num_serv, particle.r, particle.SYN, particle.PRE)
    test_par = deepcopy(particle)
    first_obj = objective_value(test_par)
    origin_dis = Inf
    while objective_value(test_par) < origin_dis
        origin_dis = objective_value(test_par)
        test_par = swap(test_par, list=list)
        test_par = move(test_par, list=list)
    end
    println("apply local search $(@sprintf("%6.2f", first_obj)) => $(@sprintf("%6.2f", objective_value(test_par)))")
    return test_par
end


function PSO(Name::String; num_par=15, max_iter=150)
    particles = []
    # obj = []

    particles = [generate_particles(Name) for _ in 1:num_par]
    # obj = objective_value.(particles)
    obj = objective_value.(particles)

    # find best particle
    best_index = argmin(obj)
    best_par = particles[best_index]


    # defind parameters
    iter = 1
    old_best = Inf
    new_best = Inf
    not_improve = 1

    # loop
    while iter < max_iter && not_improve < 10

        # save objective_value
        # old_best = new_best

        # local search
        for i in 1:num_par
        # for i in Iterators.filter(x->x!=best_index, 1:num_par)
            println("apply local search to particle $i !!!!")
            particles[i] = local_search(particles[i], best_par)
        end
        # particles = [local_search(particles[i]) for i in 1:num_par]

        # find best particle
        obj = objective_value.(particles)
        best_index = argmin(obj)

        if obj[best_index] < new_best
            new_best = obj[best_index]
            best_par = particles[best_index]
            not_improve = 1
        else
            not_improve += 1
        end

        # if not_improve == 4
        #     particles = [generate_particles(Name) for _ in 1:num_par]
        # end
        
        # if old_best - new_best < 1e4
        #     not_improve += 1
        # end

        println("iter: $iter best[$best_index]: $(@sprintf("%.2f", new_best))")
        iter += 1
    end

    # save solution
    save_particle(best_par, Name)

end


function location_simulation(instance_name::AbstractString; initial=false)
    if initial
        return mkpath("data/simulations/$instance_name/initial")
    else
        return mkpath("data/simulations/$instance_name")
    end
end


function loop_test(num, Name)
    for _ in 1:num
        par = generate_particles(Name)
        local_search(par)
    end
end