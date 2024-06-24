struct Ins
    name::String
    num_node::Int64
    num_vehi::Int64
    num_serv::Int64
    serv_a::Tuple
    serv_r::Dict
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


mutable struct Sol
    route::Array
    slot::Dict
    ins::Ins
end

"""
---
use to load all data from instance including 

number of node, vehicle service, ...

---

## `Example:`

```julia
ins = load_ins("ins10-1")
```
"""
function load_ins(Name::String)
    num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l = load_data(Name)

    # find Synchronization and Precedence constraints
    serv_r = find_service_request(r)
    serv_a = find_compat_vehicle_node(a, r)
    SYN = find_SYN(serv_r, mind, maxd)
    PRE = find_PRE(serv_r, mind, maxd)
    return Ins(Name, num_node, num_vehi, num_serv, serv_a, serv_r, mind, maxd, a, r, d, p, e, l, PRE, SYN)
end

"""
```julia
    sortperm(ins::Ins)
```

---

Example:

```julia
    ins = load_ins("ins10-1")
    sl = sortperm(ins)
```

use to sort the upper time window of instance

"""
function sortperm_l(ins::Ins)
    setdiff(Base.sortperm(ins.l), [1])
end


function sortperm_e(ins::Ins)
    setdiff(Base.sortperm(ins.e), [1])
end

"""
```julia
    all_node_serv(ins::Ins)
```

find all node and service needed to schedule in the from

Vector[[node, service]]
[[1, 1], [10, 2], [5, 3], ...]

---

`Example:`

```
    ins = load_ins("ins10-1")
    all_ns = all_node_serv(ins)
```
"""
function all_node_serv(ins::Ins)
    all_ns = Vector{Int64}[]
    push!(all_ns, [1, 1])
    for i in 2:ins.num_node
        for j in findall(x->x==1, ins.r[i, :])
            push!(all_ns, [i, j])
        end
    end
    # findall(BitMatrix(ins.r))
    return all_ns
end


"""
```julia
    generate_empty_sol(ins::Ins)
```

---

use to generate empty solution of instance; the route = [[[1, 1]], 
                                                        [[1, 1]],
                                                         .
                                                         .
                                                         .
                                                         [[1, 1]]]

---

Example:

```julia
    ins = load_ins("ins10-1")
    sol = generate_empty_sol(ins)
```
"""
function generate_empty_sol(ins::Ins)
    route = [[[1, 1]] for _ in 1:ins.num_vehi]
    slot = Dict{Int64, Vector{Int64}}()
    for syn in ins.SYN
        slot[syn[1]] = Int64[syn[2], syn[3]]
    end
    
    for pre in ins.PRE
        slot[pre[1]] = Int64[pre[2], pre[3]]
    end
    return Sol(route, slot, ins)
end


function generate_particles(ins::Ins)
    sol = generate_empty_sol(ins)
    sol = insert_node_service(sol)
    return sol
end


function inserted_node(sol::Sol)
    an = all_node_serv(sol.ins)
    all_node = Vector{Int64}[]
    for i in 1:sol.ins.num_vehi
        append!(all_node, sol.route[i])
    end

    setdiff!(an, all_node)

    return an
end


function load_example1()
    route, slot = example()
    ins = load_ins("ins10-1")
    sol = generate_empty_sol(ins)
    sol.route = route
    return sol
end


function swap(sol::Sol, list)
    input_route = deepcopy(sol.route)
    for ls in list
        num_v1, num_loca1 = find_location_by_node_service(input_route, ls[1][1], ls[1][2])
        num_v2, num_loca2 = find_location_by_node_service(input_route, ls[2][1], ls[2][2])
        test_route = deepcopy(input_route)
        test_route[num_v1][num_loca1], test_route[num_v2][num_loca2] = test_route[num_v2][num_loca2], test_route[num_v1][num_loca1]
        if compatibility(test_route, sol.slot, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.PRE, sol.ins.SYN) && in_same_route(test_route)
            # test_st = starttime(test_route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
            # st = starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
            # test_st = try starttime(sol, test_route) catch StackOverflowError; continue end
            # st = try starttime(sol, input_route) catch StackOverflowError; continue end 
            test_st = try starttime(sol, test_route) catch StackOverflowError;continue end
            st = starttime(sol, input_route)
            if objective_value(test_route, test_st, sol.ins.p, sol.ins.l, sol.ins.d) < objective_value(input_route, st, sol.ins.p, sol.ins.l, sol.ins.d) && check_PRE(test_route, test_st, sol.ins.maxd, sol.ins.PRE) && check_SYN(test_route, sol.ins.SYN)
                input_route = deepcopy(test_route)
            end
        end
    end
    sol.route = input_route
    return sol
end


function move(sol::Sol, list)
    input_route = deepcopy(sol.route)
    for ls in list
        num_v1, num_loca1 = find_location_by_node_service(input_route, ls[1][1], ls[1][2])
        num_v2, num_loca2 = find_location_by_node_service(input_route, ls[2][1], ls[2][2])

        # if num_v1 == num_v2
        #     continue
        # end

        test_route = deepcopy(input_route)

        #  move
        moved_item = splice!(test_route[num_v1], num_loca1)
        insert!(test_route[num_v2], num_loca2, moved_item)

        
        if compatibility(test_route, sol.slot, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.PRE, sol.ins.SYN) && in_same_route(test_route)
            # find starttime
            test_st = try starttime(sol, test_route) catch StackOverflowError; println("StackOverflowError");continue end
            st = starttime(sol, input_route)

            if objective_value(test_route, test_st, sol.ins.p, sol.ins.l, sol.ins.d) < objective_value(input_route, st, sol.ins.p, sol.ins.l, sol.ins.d) && check_PRE(test_route, test_st, sol.ins.maxd, sol.ins.PRE) && check_SYN(test_route, sol.ins.SYN)
                # println("cost reduce in Move")
                input_route = deepcopy(test_route)
            end
        end
    end
    sol.route = input_route
    return sol
end


function local_search(particle::Sol, best_par::Sol)
    test_par = deepcopy(particle)
    list = List(particle.ins.num_node, particle.ins.num_vehi, particle.ins.num_serv, particle.ins.r, particle.ins.SYN, particle.ins.PRE)
    # first_obj = objective_value(test_par)
    test_par = swap(test_par, list)
    test_par = move(test_par, list)
    test_par = path_relinking(test_par, best_par)
    return test_par
end


function total_distance(sol::Sol)
    return total_distance(sol.route, sol.ins.d)
end


function tardiness(sol::Sol)
    return tardiness(sol.route, starttime(sol), sol.ins.p, sol.ins.l)
end


function objective_value(sol::Sol)
    return objective_value(sol.route, starttime(sol), sol.ins.p, sol.ins.l, sol.ins.d)
end


function find_starttime(sol::Sol)
    # st = initial_starttime(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv)
    find_starttime(sol.route, sol.slot, sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.mind, sol.ins.maxd, sol.ins.a, sol.ins.r, sol.ins.d, sol.ins.p, sol.ins.e, sol.ins.l, sol.ins.PRE, sol.ins.SYN)
end


function insert_node_sevice(sol::Sol, node::Int64, serv::Int64)
    nothing
end


function find_before_after_serv(ins::Ins, node::Int64)
    if any(in.(node, [ins.PRE[i][1] for i in 1:length(ins.PRE)]))
        index_pre = findfirst(in.(node, [ins.PRE[i][1] for i in 1:length(ins.PRE)]))
        return ins.PRE[index_pre][2:3]
    elseif any(in.(node, [ins.SYN[i][1] for i in 1:length(ins.SYN)]))
        index_syn = findfirst(in.(node, [ins.SYN[i][1] for i in 1:length(ins.SYN)]))
        return ins.SYN[index_syn][2:3]
    else
        return findall(x->x==1.0, ins.r[node, :])
    end
    throw("find_before_after_serv => Not found node $node in PRE or SYN")
end
"""
```julia
    insert_node_service(test_sol::Sol)
```

---

use to insert all node service to route in sol

"""
function insert_node_service(test_sol::Sol)
    sol = deepcopy(test_sol)
    remain_node = inserted_node(sol)
    sorted_node = sortperm_e(sol.ins)
    for node in sorted_node
        if any(in.(node, remain_node))
            com_serv = findall(x->x==1, sol.ins.r[node, :])
            if length(com_serv) == 1
                compat_vehicle = rand(findall(x->x==1, sol.ins.a[:, com_serv]))
                insert!(sol.route[compat_vehicle], length(sol.route[compat_vehicle]), [node, com_serv[1]])
                # println("node: $node, has 1 service")
            else
                # println("node: $node, has 2 services")
                serv_all = [i for i in find_before_after_serv(sol.ins, node)]

                com_vehi1 = findall(x->x==1, sol.ins.a[:, serv_all[1]])
                com_vehi2 = findall(x->x==1, sol.ins.a[:, serv_all[2]])
                # com_vehi = []
                if isempty(intersect(com_vehi1, com_vehi2))
                    com_vehi = [rand(com_vehi1), rand(com_vehi2)]
                elseif isempty(setdiff(com_vehi1, com_vehi2))
                    com2 = rand(com_vehi2)
                    com1 = rand(setdiff(com_vehi1, com2))
                    com_vehi = [com1, com2]
                elseif isempty(setdiff(com_vehi2, com_vehi1))
                    com1 = rand(com_vehi1)
                    com2 = rand(setdiff(com_vehi2, com1))
                    com_vehi = [com1, com2]
                elseif !isempty(setdiff(com_vehi1, com_vehi2)) && !isempty(setdiff(com_vehi2, com_vehi1))
                    com1 = rand(setdiff(com_vehi1, com_vehi2))
                    com2 = rand(setdiff(com_vehi2, com_vehi1))
                    com_vehi = [com1, com2]
                end

                # insert 2 services
                insert!(sol.route[com_vehi[1]], length(sol.route[com_vehi[1]]), [node, serv_all[1]])
                insert!(sol.route[com_vehi[2]], length(sol.route[com_vehi[2]]), [node, serv_all[2]])
            end
        end
    end
    return sol
end


function swap(sol::Sol)
    list = List(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.r, sol.ins.SYN, sol.ins.PRE)
    sol.route = swap(sol.route, sol.slot, sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.mind, sol.ins.maxd, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.d, sol.ins.p, sol.ins.e, sol.ins.l, sol.ins.PRE, sol.ins.SYN, list)
    return sol
end


function List(sol::Sol)
    list = List(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.r, sol.ins.SYN, sol.ins.PRE)
end


# function swap(sol::Sol, list)
#     sol.route = swap(sol.route, sol.slot, sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.mind, sol.ins.maxd, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.d, sol.ins.p, sol.ins.e, sol.ins.l, sol.ins.PRE, sol.ins.SYN, list)
#     return sol
# end


function move(sol::Sol)
    list = List(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.r, sol.ins.SYN, sol.ins.PRE)
    sol.route = move(sol.route, sol.slot, sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.mind, sol.ins.maxd, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.d, sol.ins.p, sol.ins.e, sol.ins.l, sol.ins.PRE, sol.ins.SYN, list)
    return sol
end


function random_move(sol::Sol)
    list = List(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.r, sol.ins.SYN, sol.ins.PRE)
    sol.route = random_move(sol.route, sol.slot, sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv, sol.ins.mind, sol.ins.maxd, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.d, sol.ins.p, sol.ins.e, sol.ins.l, sol.ins.PRE, sol.ins.SYN, list)
    return sol
end


function path_relinking_insert(input_route::Vector, remain_node::Vector, a::Matrix{Int64})
    for ns in remain_node
        test_route = deepcopy(input_route)
        choose_vehi = rand(findall(x->x==1, a[:, ns[2]]))
        len_route = length(test_route[choose_vehi])
        choose_posi = randcycle(len_route)

        terminate = false
        posi = 1
        while !terminate && posi <= len_route
            posi = choose_posi[posi]
            insert!(test_route[choose_vehi], posi, ns)
            posi += 1

            # try calculate starttime
            st = try starttime(test_route) catch e; continue end
            inout_route = deepcopy(test_route)

            terminate = true
        end
    end
    return input_route
end


function path_relinking(route::Array, best_route::Array, a::Matrix, sol::Sol)
    test_route = deepcopy(route)

    rand_vehi = rand(1:length(route))
    remove_route = test_route[rand_vehi]
    remove_best_route = best_route[rand_vehi]
    
    remain_node = setdiff(remove_route, remove_best_route)

    # insert best route
    test_route[rand_vehi] = best_route[rand_vehi]

    # insert
    test_route = path_relinking_insert(test_route, remain_node, a)

    # check
    insert_all = true
    for ns in remain_node
        if !check_assigned_node(test_route, ns[1], ns[2], length(test_route))
            insert_all = false
            break
        end
    end
    if insert_all
        # println("can apply path relinking!!!")
        sol.route = test_route
        try starttime(sol) catch e; return route end
        return test_route
    else
        return route
    end
end


function path_relinking(par::Sol, par_best::Sol)
    # best_route = rand(par_best.route)
    route = path_relinking(par.route, par_best.route, par.ins.a, par)
    par.route = route
    return par
end


function initial_starttime(sol::Sol)
    initial_starttime(sol.ins.num_node, sol.ins.num_vehi, sol.ins.num_serv)
end


function find_other_serv_in_syn_pre(node, serv, SET; is_pre=false)
    pre = SET[findfirst(x->x[1]==node, SET)]
    other_serv = setdiff(pre, [node, serv])[1]
    if is_pre
        if findfirst(pre[2:end] .== serv) < findfirst(pre[2:end] .== other_serv)
            return other_serv, false
        else
            return other_serv, true
        end
    else
        return other_serv
    end
end


function check_SYN(sol::Sol)
    check_SYN(sol.route, sol.ins.SYN)
end


function check_PRE(sol::Sol)
    check_PRE(sol.route ,starttime(sol), sol.ins.maxd, sol.ins.PRE)
end


function compatibility(sol::Sol)
    compatibility(sol.route, sol.slot, sol.ins.a, sol.ins.r, sol.ins.serv_a, sol.ins.serv_r, sol.ins.PRE, sol.ins.SYN)
end


function forward_starttime(sol::Sol, st::Dict, start_lo::Int64, stop_lo::Int64, vehi::Int64)
    min_end = minimum([stop_lo, length(sol.route[vehi])])
    for lo in start_lo+1:min_end

        # last node and service
        bnode = sol.route[vehi][lo-1][1]
        bserv = sol.route[vehi][lo-1][2]

        # current node and service
        node = sol.route[vehi][lo][1]
        serv = sol.route[vehi][lo][2]

        arrival_time = st[bnode][vehi, bserv] + sol.ins.p[vehi, bserv, bnode] + sol.ins.d[bnode, node]

        if arrival_time < sol.ins.e[node] # when starttime is not change
            return st
        elseif in_SYN(node, sol.ins.SYN)
            st[node][vehi, serv] = arrival_time
            other_serv = find_other_serv_in_syn_pre(node, serv, sol.ins.SYN)
            ovehi, oloca = find_location_by_node_service(sol.route, node, other_serv)
            st = forward_starttime(sol, st, oloca, lo, ovehi)
        elseif in_PRE(node, sol.ins.PRE)
            st[node][vehi, serv] = arrival_time
            other_serv, before_serv = find_other_serv_in_syn_pre(node, serv, sol.ins.PRE, is_pre=true)
            ovehi, oloca = find_location_by_node_service(sol.route, node, other_serv)
            
            if before_serv
                if st[node][vehi, serv] - st[node][ovehi, other_serv] > sol.ins.maxd[node]
                    st[node][ovehi, other_serv] = st[node][vehi, serv] - sol.ins.maxd[node]
                    st = forward_starttime(sol, st, oloca, lo, ovehi)
                else
                    return st
                end
            else
                if st[node][ovehi, other_serv] - st[node][vehi, serv] < sol.ins.mind[node]
                    st[node][ovehi, other_serv] = sol.ins.mind[node] + st[node][vehi, serv]
                    st = forward_starttime(sol, st, oloca, lo, ovehi)
                else
                    return st
                end
            end
        else
            st[node][vehi, serv] = arrival_time
        end
    end
    return st
end


# function find_starttime(route, slot, num_node, num_vehi, num_serv, mind, maxd, a, r, d, p, e, l, PRE, SYN)
#     st = initial_starttime(num_node, num_vehi, num_serv)

#     # find max column
#     len_route = [length(route[i]) for i in 1:num_vehi]
#     maxcolumn = maximum(len_route)

#     calculated = Int64[]

#     for column in 1:maxcolumn
#         for v in 1:num_vehi
#             if column <= len_route[v]

#                 # find departure and arrival node
#                 if column == 1
#                     bnode = 1
#                     bserv = 1
#                 else
#                     bnode = route[v][column-1][1]
#                     bserv = route[v][column-1][2]
#                 end
#                 node = route[v][column][1]
#                 serv = route[v][column][2]
                
#                 min_st = st[bnode][v, bserv] + d[bnode, node] + p[v, bserv, bnode]
#                 if min_st < e[node]
#                     min_st = e[node]
#                 end


#                 println("calculate column $column vehicle $v node $(route[v][column])")
#                 println("e = $(e[node]), d = $(d[bnode, node]) arrive time = $(st[bnode][v, bserv] + d[bnode, node] + p[v, bserv, bnode])")

#                 if in_SYN(node, SYN) && in(node, calculated)
#                     # find location of syn node
#                     other_serv = find_other_serv_in_syn_pre(node, serv, SYN)
#                     ovehi, oloca = find_location_by_node_service(route, node, other_serv)
#                     if min_st < st[node][ovehi, other_serv]
#                         # println("update min_st")
#                         st[node][v, serv] = st[node][ovehi, other_serv]
#                     elseif min_st > st[node][ovehi, other_serv]
#                         st[node][v, serv] = min_st
#                         st[node][ovehi, other_serv] = min_st
#                         st = forward_starttime(sol, st, oloca, column, ovehi)
#                     else
#                         st[node][v, serv] = min_st
#                     end
#                 elseif in_PRE(node, PRE) && in(node, calculated)

#                     # find location of pre node
#                     other_serv, before_serv = find_other_serv_in_syn_pre(node, serv, PRE, is_pre=true)
#                     ovehi, oloca = find_location_by_node_service(route, node, other_serv)

#                     if before_serv # other service has to process first
#                         if min_st - st[node][ovehi, other_serv] < mind[node]
#                             st[node][v, serv] = st[node][ovehi, other_serv] + mind[node]
#                         elseif min_st - st[node][ovehi, other_serv] > maxd[node]
#                             st[node][v, serv] = min_st
#                             st[node][ovehi, other_serv] = min_st - maxd[node]
#                             st = forward_starttime(sol, st, oloca, column, ovehi)
#                         end
#                     else
#                         if st[node][ovehi, other_serv] - min_st < mind[node]
#                             st[node][v, serv] = min_st
#                             st[node][ovehi, other_serv] = min_st + mind[node]
#                             st = forward_starttime(sol, st, oloca, column, ovehi)
#                         elseif st[node][ovehi, other_serv] - min_st > maxd[node]
#                             st[node][v, serv] = st[node][ovehi, other_serv] - maxd[node]
#                         end
#                     end
#                 else
#                     st[node][v, serv] = min_st
#                 end


#                 push!(calculated, node)
#             end
#             # println(" ")
#         end
#     end
#     return st
# end


function starttime(sol::Sol, route::Array)
    sol.route = route
    return starttime(sol)
end


function starttime(sol::Sol)
    st = initial_starttime(sol)

    # find max column
    len_route = [length(sol.route[i]) for i in 1:sol.ins.num_vehi]
    maxcolumn = maximum(len_route)

    calculated = Int64[]

    for column in 1:maxcolumn
        for v in 1:sol.ins.num_vehi
            if column <= len_route[v]

                # find departure and arrival node
                if column == 1
                    bnode = 1
                    bserv = 1
                else
                    bnode = sol.route[v][column-1][1]
                    bserv = sol.route[v][column-1][2]
                end
                node = sol.route[v][column][1]
                serv = sol.route[v][column][2]
                
                min_st = st[bnode][v, bserv] + sol.ins.d[bnode, node] + sol.ins.p[v, bserv, bnode]
                if min_st < sol.ins.e[node]
                    min_st = sol.ins.e[node]
                end


                # println("calculate column $column vehicle $v node $(sol.route[v][column])")
                # println("e = $(sol.ins.e[node]), d = $(sol.ins.d[bnode, node]) arrive time = $(st[bnode][v, bserv] + sol.ins.d[bnode, node] + sol.ins.p[v, bserv, bnode])")

                if in_SYN(node, sol.ins.SYN) && in(node, calculated)
                    # find location of syn node
                    other_serv = find_other_serv_in_syn_pre(node, serv, sol.ins.SYN)
                    ovehi, oloca = find_location_by_node_service(sol.route, node, other_serv)
                    if min_st < st[node][ovehi, other_serv]
                        # println("update min_st")
                        st[node][v, serv] = st[node][ovehi, other_serv]
                    elseif min_st > st[node][ovehi, other_serv]
                        st[node][v, serv] = min_st
                        st[node][ovehi, other_serv] = min_st
                        st = forward_starttime(sol, st, oloca, column, ovehi)
                    else
                        st[node][v, serv] = min_st
                    end
                elseif in_PRE(node, sol.ins.PRE) && in(node, calculated)

                    # find location of pre node
                    other_serv, before_serv = find_other_serv_in_syn_pre(node, serv, sol.ins.PRE, is_pre=true)
                    ovehi, oloca = find_location_by_node_service(sol.route, node, other_serv)

                    if before_serv # other service has to process first
                        if min_st - st[node][ovehi, other_serv] < sol.ins.mind[node]
                            st[node][v, serv] = st[node][ovehi, other_serv] + sol.ins.mind[node]
                        elseif min_st - st[node][ovehi, other_serv] > sol.ins.maxd[node]
                            st[node][v, serv] = min_st
                            st[node][ovehi, other_serv] = min_st - sol.ins.maxd[node]
                            st = forward_starttime(sol, st, oloca, column, ovehi)
                        else
                            st[node][v, serv] = min_st
                        end
                    else
                        if st[node][ovehi, other_serv] - min_st < sol.ins.mind[node]
                            st[node][v, serv] = min_st
                            st[node][ovehi, other_serv] = min_st + sol.ins.mind[node]
                            st = forward_starttime(sol, st, oloca, column, ovehi)
                        elseif st[node][ovehi, other_serv] - min_st > sol.ins.maxd[node]
                            st[node][v, serv] = st[node][ovehi, other_serv] - sol.ins.maxd[node]
                        else
                            st[node][v, serv] = min_st
                        end 
                    end
                else
                    st[node][v, serv] = min_st
                end


                push!(calculated, node)
            end
            # println(" ")
        end
    end
    return st
end


function PSO(ins::Ins; num_par=15, max_iter=50)
    particles = []
    obj = []

    particles = [generate_particles(ins) for _ in 1:num_par]
    obj = objective_value(particles)

    # find best particle
    best_index = argmin(obj)
    best_par = particles[best_index]


    # defind parameters
    iter = 1
    old_best = Inf
    new_best = Inf
    not_improve = 1

    # find current number of run
    location = "$(location_simulation(ins.name, initial=false))"
    num = length(glob("$(ins.name)*.jld2", location))+1
    io = open(joinpath(dir(), "data", "simulations", ins.name, "$(ins.name)-$num.csv"), "a")
    write(io, "iter,obj,time\n")
    close(io)
    # loop
    while iter < max_iter && not_improve < 5
        t = @elapsed begin
        # save objective_value
        # old_best = new_best

        # local search
        for i in 1:num_par
        # for i in Iterators.filter(x->x!=best_index, 1:num_par)
            println("apply local search to particle $i")
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
        end

        println("iter: $iter best[$best_index]: $(@sprintf("%.2f", new_best)), PRE: $(check_PRE(best_par)), SYN: $(check_SYN(best_par)), Compat: $(compatibility(best_par))")
        io = open(joinpath(dir(), "data", "simulations", ins.name, "$(ins.name)-$num.csv"), "a")
        write(io, "$iter,$(@sprintf("%.2f", new_best)),$(@sprintf("%.2f", t))\n")
        close(io)
        iter += 1
    end

    # save solution
    save_particle(best_par)

    tab = CSV.File(joinpath(dir(), "data", "simulations", ins.name, "$(ins.name)-$num.csv")) |> DataFrame
    sent_email_report("report2014", "iter: $iter, $(@sprintf("%.2f", new_best)), PRE: $(check_PRE(best_par)), SYN: $(check_SYN(best_par)), Compat: $(compatibility(best_par))\n", df=tab)
    ic = open(joinpath(dir(), "data", "simulations", ins.name, "$(ins.name)-$num.md"), "w")
    write(ic, "$(latexify(tab, env=:mdtable))\n")
    close(ic)
end


function find_min_objective(ins_name::AbstractString)
    min_obj = []
    min_iter = []
    min_time = []
    location = "$(location_simulation(ins_name, initial=false))"
    for fo in glob("$(ins_name)*.csv", location)
        df = CSV.File(fo) |> DataFrame
        # println("iter: $(size(df, 1)) $(minimum(df[!, 2])) time: $(sum(df[!, 3]))")
        push!(min_iter, size(df, 1))
        try push!(min_obj, minimum(df[!, 2])) catch e; push!(min_obj, 0.0) end
        try push!(min_time, sum(df[!, 3])) catch e; push!(min_time, 0.0) end
        # push!(min_time, sum(df[!, 3]))
    end
    println("TOTAL: $(length(min_obj))")
    if length(min_iter) == 0
        return (0,0,0,0)
    else
        println("ins: $ins_name")
        println("min obj: $(minimum(min_obj))")
        println("min iter: $(minimum(min_iter))")
        println("min time: $(minimum(min_time))")
        return (minimum(min_iter), minimum(min_obj), minimum(min_time), minimum(min_time)/60)
    end
end


function natural(x, y)
    k(x) = [occursin(r"\d+", s) ? parse(Int, s) : s 
            for s in split(replace(x, r"\d+" => s->" $s "))]
    A = k(x); B= k(y)    
    for (a, b) in zip(A, B)
        if !isequal(a, b)
            return typeof(a) <: typeof(b) ? isless(a, b) :
                   isa(a,Int) ? true : false
        end
    end
    return length(A) < length(B)
end


"""
    find_ins(ins_name::AbstractString)

TBW
"""
function find_ins(ϵ::AbstractString)
    dm = []
    location = joinpath(splitdir(splitdir(Base.find_package("VehicleRoutingSynPre"))[1])[1], "data", "simulations")
    println("pwd: $(pwd())")
    println("location: $location")
    for α in sort(glob("$ϵ*", location), lt=natural)
        for β in glob("*.csv", α)
            df = CSV.File(β) |> DataFrame
            if size(df, 1) > 5
                @show min_iter = size(df, 1)
                @show min_obj = minimum(df[!, 2])
                @show min_time = sum(df[!, 3])
                push!(dm, (splitdir(α)[2], min_obj, min_iter, min_time))
            end
        end
    end
    return dm
end


function create_csv_2014()
    loca = joinpath(splitdir(splitdir(Base.find_package("VehicleRoutingSynPre"))[1])[1], "data", "table")
    df = CSV.File(joinpath(loca, "table.csv")) |> DataFrame
    io = open(joinpath(loca, "our.csv"), "w")
    write(io, "Instance,Obj2014,obj,numiteration,time(second),time(min)\n")
    for (i, ins_name) in enumerate(df[!, 1])
        min_iter, min_obj, min_time, min_time60 = find_min_objective(ins_name)
        obj2014 = try minimum([df[i, 3], df[i, 6], df[i,7], df[i, 9], df[i, 11], df[i, 13]]) catch e; minimum([df[i, 6], df[i, 7], df[i, 9], df[i, 11], df[i, 13]]) end
        mintime = try @sprintf("%.2f", min_time) catch e; 0.0 end
        mintime60 = try @sprintf("%.2f", min_time60) catch e; 0.0 end
        write(io, "$ins_name,$(obj2014),$(min_obj),$(min_iter),$(mintime),$(mintime60)\n")
    end
    dm = find_ins("ins500")
    for dd in dm
        mintime = try @sprintf("%.2f", dd[4]) catch e; 0.0 end
        mintime60 = try @sprintf("%.2f", dd[4]/60) catch e; 0.0 end
        write(io, "$(dd[1]),$(Inf),$(dd[2]),$(dd[3]),$(mintime),$(mintime60)\n")
    end
    close(io)
    return CSV.File(joinpath(loca, "our.csv")) |> DataFrame
end


function sent_email_report(file_name::String, massage; df=nothing)
    create_csv_2014()
    loca = splitdir(splitdir(Base.find_package("VehicleRoutingSynPre"))[1])[1]
    weave(joinpath(loca, "report", "$file_name.Jmd"), doctype="md2html")
    attachments = [joinpath(loca, "report", "$file_name.html")]
    subject = "report"
    sent_email(subject, massage, df=df, attachments=attachments)
end


function sent_email(subject::String, massage; df=nothing, attachments=nothing)
    username = "payakorn.sak@gmail.com"
    opt = SendOptions(
    isSSL = true,
    username = "payakorn.sak@gmail.com",
    passwd = "cdtcdmxydxihuroo")
    # passwd = "daxdEw-kyrgap-2bejge")
    #Provide the message body as RFC5322 within an IO
    if isnothing(df)
        msg = "$massage"
    else
        msg = Markdown.parse(
            """
            $massage

            $(latexify(df, env=:mdtable, latex=false))

            """
        )
    end

    """
    Example:
    ## Julia in a Nutshell

        1. **Fast** - Julia was designed from the beginning for [high performance](https://docs.julialang.org/en/v1/manual/types/).
        1. **Dynamic** - Julia is [dynamically typed](https://docs.julialang.org/en/v1/manual/types/).
        1. **Reproducible** - recreate the same [Julia environment](https://julialang.github.io/Pkg.jl/v1/environments/) every time.
        1. **Composable** - Julia uses [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/) as a paradigm.
        1. **General** - One can build entire [Applications and Microservices](https://www.youtube.com/watch?v=uLhXgt_gKJc) in Julia.
        1. **Open source** - Available under the [MIT license](https://github.com/JuliaLang/julia/blob/master/LICENSE.md), with the [source code](https://github.com/JuliaLang/julia) on GitHub.

        It has *over 5,000* [Julia packages](https://juliahub.com/ui/Packages) and a *variety* of advanced ecosystems. Check out more on [the Julia Programing Language website](https://julialang.org).
    """

    msg = get_mime_msg(msg)
    if isnothing(attachments)
        body = get_body(["<payakornn@gmail.com>"], "You <$username>", subject, msg)
    else
        body = get_body(["<payakornn@gmail.com>"], "You <$username>", subject, msg; attachments)
    end
    # body = IOBuffer(
    # # "Date: Fri, 18 Oct 2013 21:44:29 +0100\r\n" *
    # "From: You <$username>\r\n" *
    # "To: payakornn@gmail.com\r\n" *
    # "Subject: $subject\r\n" *
    # "\r\n" *
    # "$massage\r\n")
    url = "smtps://smtp.gmail.com:465"
    rcpt = ["<payakornn@gmail.com>"]
    from = "<$username>"
    resp = send(url, rcpt, from, body, opt)
end


function load_example1_par()
    par = generate_particles("ins10-1")
    route, slot = example()
    par.route = route
    par.slot = slot
    par.starttime = find_starttime(par)
    return par
end


function dist(x1::Int64, x2::Int64, y1::Int64, y2::Int64)
    x1 = float(x1)
    x2 = float(x2)
    y1 = float(y1)
    y2 = float(y2)
    dist(x1, x2, y1, y2)
end


function dist(x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    sqrt((x1-x2)^2 + (y1-y2)^2)
end


function dist(x::Vector, y::Vector)
    [dist(x[i], x[j], y[i], y[j]) for i in 1:length(x), j in 1:length(y)]
end


function load_ins_text100(num_node::Int64, num_ins::Int64)

    location = joinpath(@__DIR__, "..", "data", "raw_HHCRSP", "HHCRSP")
    s = open(joinpath(location, "InstanzVNS_HCSRP_$(num_node)_$(num_ins).txt"))
    lines = readlines(s)
    close(s)

    variables = ["nbServi",
                "nbVehi",
                "nbCust",
                "nbSynch",
                "e",
                "l",
                "cx",
                "cy",
                "s",
                "delta",
                "p",
                "att"]
    
    line_var = Int64[]
    for (i, var) in enumerate(variables)
        if !isempty(findall(occursin.(var, lines)))
            println("line: $(findfirst(occursin.("$var", lines[i:end]))+i-1) has $(var)")
            push!(line_var, findfirst(occursin.("$var", lines[i:end]))+i-1)
        end
    end

    # parse number of nodes, vehicles, services
    num_serv = parse(Int64, split(replace(lines[1], ";" => ""), "=")[2])
    num_vehi = parse(Int64, split(replace(lines[2], ";" => ""), "=")[2])
    num_cust = parse(Int64, split(replace(lines[3], ";" => ""), "=")[2])
    num_sync = parse(Int64, split(replace(lines[4], ";" => ""), "=")[2])
    num_node = num_cust + num_sync

    e =     parse.(Int64, split(lines[line_var[5]+1:line_var[6]-1][1]))
    l =     parse.(Int64, split(lines[line_var[6]+1:line_var[7]-1][1]))
    cx =    parse.(Int64, split(lines[line_var[7]+1:line_var[8]-1][1]))
    cy =    parse.(Int64, split(lines[line_var[8]+1:line_var[9]-1][1]))

    # r
    r = []
    rr = lines[(line_var[9] + 1):(line_var[10]-1)]
    for (i, j) in enumerate(rr)
        j = String.(split(j))
        push!(r, parse.(Int64, j))
    end
    r = [r[i][j] for i in 1:num_node+1, j in 1:num_serv]


    # delta
    delta = []
    ddelta = lines[(line_var[10] + 1):(line_var[11]-1)]
    for (i, j) in enumerate(ddelta)
        j = String.(split(j))
        push!(delta, parse.(Int64, j))
    end
    delta = [delta[i][j] for i in 1:num_node+1, j in 1:2]
    mind = delta[:, 1]
    maxd = delta[:, 2]

    # p_text = lines[line_var[11]+1:line_var[12]-1]
    p_text = parse.(Int64, split(lines[line_var[11] + num_vehi + 1]))[1]
    p = zeros(Float64, num_vehi, num_serv)
    for i in 2:num_node+1
        p = cat(p, p_text*ones(Float64, num_vehi, num_serv), dims=3)
    end

    # a
    a = []
    att =   lines[line_var[12]+1:end]
    for (i, j) in enumerate(att)
        j = String.(split(j))
        push!(a, parse.(Int64, j))
    end
    a = [a[i][j] for i in 1:num_vehi, j in 1:num_serv]

    d = dist(cx, cy)

    # find Synchronization and Precedence constraints
    serv_r = find_service_request(r)
    serv_a = find_compat_vehicle_node(a, r)
    SYN = find_SYN(serv_r, mind, maxd)
    PRE = find_PRE(serv_r, mind, maxd)
    return Ins("ins$(num_node)-$num_ins", num_node+1, num_vehi, num_serv, serv_a, serv_r, mind, maxd, a, r, d, p, e, l, PRE, SYN)
end


function load_ins_text(num_node::Int64, num_ins::Int64; locc=nothing)

    if num_node >= 100
        load_ins_text100(num_node, num_ins)
    else
        
        location = joinpath(@__DIR__, "..", "data", "raw_HHCRSP", "HHCRSP")
        s = open(joinpath(location, "InstanzCPLEX_HCSRP_$(num_node)_$(num_ins).txt"))
        lines = readlines(s)
        close(s)
        # for (i, line) in enumerate(lines)
        #     println("$i: $line")
        # end

        variables = ["nbNodes",
                    "nbVehi",
                    "nbServi",
                    "r",
                    "DS",
                    "a",
                    "x",
                    "y",
                    "d",
                    "p",
                    "mind",
                    "maxd",
                    "e",
                    "l"
        ]

        line_var = Int64[]
        for var in variables
            if !isempty(findall(occursin.(var, lines)))
                println("line: $(findfirst(occursin.("$var=", lines))) has $(var)=")
                push!(line_var, findfirst(occursin.("$var=", lines)))
            end
        end

        # parse number of nodes, vehicles, services
        num_node = parse(Int64, split(replace(lines[1], ";" => ""), "=")[2]) - 1
        num_vehi = parse(Int64, split(replace(lines[2], ";" => ""), "=")[2])
        num_serv = parse(Int64, split(replace(lines[3], ";" => ""), "=")[2])



        # r
        r = []
        rr = lines[(line_var[4] + 1):(line_var[5]-2)]
        rr = replace.(rr, "[" => "")
        rr = replace.(rr, "]" => "")
        for (i, j) in enumerate(rr)
            j = String.(split(j, ","))
            push!(r, parse.(Int64, j))
        end
        r = [r[i][j] for i in 1:num_node+1, j in 1:num_serv]
        
        # a
        a = []
        @show aa = lines[(line_var[6]):(line_var[7]-2)]
        @show aa = replace.(aa, "a=" => "")
        @show aa = replace.(aa, "[" => "")
        @show aa = replace.(aa, "]" => "")
        for (i, j) in enumerate(aa)
            j = String.(split(j, ","))
            push!(a, parse.(Int64, j))
        end
        @show a = [a[i][j] for i in 1:num_vehi, j in 1:num_serv]


        # d
        d = []
        dd = lines[(line_var[9]+1):(line_var[10]-2)]
        dd = replace.(dd, "a=" => "")
        dd = replace.(dd, "[" => "")
        dd = replace.(dd, "]" => "")
        for (i, j) in enumerate(dd)
            j = String.(split(j, ","))
            push!(d, parse.(Float64, j))
        end
        d = [d[i][j] for i in 1:num_node, j in 1:num_node]

        mind = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[11]], ";" => ""), "[" => ""), "]" => ""), "$(variables[11])=" => ""), ",")))
        maxd = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[12]], ";" => ""), "[" => ""), "]" => ""), "$(variables[12])=" => ""), ",")))
        e = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[13]], ";" => ""), "[" => ""), "]" => ""), "$(variables[13])=" => ""), ",")))
        l = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[14]], ";" => ""), "[" => ""), "]" => ""), "$(variables[14])=" => ""), ",")))
        
        x = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[7]], ";" => ""), "{" => ""), "}" => ""), "//$(variables[7])=" => ""), ",")))
        y = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[8]], ";" => ""), "{" => ""), "}" => ""), "//$(variables[8])=" => ""), ",")))
        DS = parse.(Int64, String.(split(replace(replace(replace(replace(lines[line_var[5]], ";" => ""), "{" => ""), "}" => ""), "$(variables[5])=" => ""), ",")))

        # p
        pp = parse.(Float64, String.(split(replace(replace(replace(lines[line_var[10]+num_vehi+2], ";" => ""), "[" => ""), "]" => ""), ",")))
        service_time = pp[1]
        p = zeros(Float64, num_vehi, num_serv)
        for i in 2:num_node
            p = cat(p, service_time*ones(Float64, num_vehi, num_serv), dims=3)
        end


        # find Synchronization and Precedence constraints
        serv_r = find_service_request(r)
        serv_a = find_compat_vehicle_node(a, r)
        SYN = find_SYN(serv_r, mind, maxd)
        PRE = find_PRE(serv_r, mind, maxd)
        if isnothing(locc)
            return Ins("ins$(num_node-1)-$num_ins", num_node, num_vehi, num_serv, serv_a, serv_r, mind, maxd, a, r, d, p, e, l, PRE, SYN)
        else
            return x, y
        end
    end
end


function save_particle(particle::Sol; initial=false)
    instance_name = particle.ins.name
    location = "$(location_simulation(instance_name, initial=initial))"
    num = length(glob("$instance_name*.jld2", location))
    save_object("$(location_simulation(instance_name, initial=initial))/$instance_name-$(num+1).jld2", particle)
end


function save_information_run(name, best_par, t, iter)
    nothing
end


function df_conclusion_table()
    CSV.File(joinpath(dir(), "data", "table", "our.csv")) |> DataFrame
end


function df_pretty_table()
    df = df_conclusion_table()
    header = (["Instance", "Obj2014", "ObjOur", "#iteration", "time", "time"], 
    ["", "", "", "", "second", "min"])
    h1 = Highlighter((df, i, j) -> (df[i, 6] > 0.5), foreground = :blue)
    pretty_table(df, show_row_number=true, header=header, highlighters=h1)
end


function report()
    stopper = Ref(false)
    task = @async while !stopper[]
        sent_email_report("report2014", "")
        sleep(21600)
    end
    return task, stopper
end


function save_particle(particle::Particle, instance_name::String; initial=false)
    location = "$(location_simulation(instance_name, initial=initial))"
    num = length(glob("$instance_name*.jld2", location))
    save_object("$(location_simulation(instance_name, initial=initial))/$instance_name-$(num+1).jld2", particle)
end


function save_ins(ins::Ins)
    opath = joinpath(dir(), "data", "simulations", ins.name)
    mkpath(opath)
    try save_object(joinpath(opath, "$(ins.name).jld2"), ins) catch e; save_ins_txt(ins) end
    save_ins_txt(ins)
end


function dir()
    splitdir(splitdir(Base.find_package("VehicleRoutingSynPre"))[1])[1]
end


function dir(p...)
    joinpath(dir(), p...)
end


function save_ins_txt(ins::Ins)
    opath = joinpath(dir(), "data", "simulations", ins.name)
    mkpath(opath)
    println("........Writing text file........")
    io = open(joinpath(opath, "$(ins.name).txt"), "w")
    write(io, "Name: $(ins.name)\n")
    write(io, "#nodes: $(ins.num_node)\n")
    write(io, "#vehicles: $(ins.num_vehi)\n")
    write(io, "#services: $(ins.num_serv)\n")
    write(io, "a: $(ins.a)\n")
    write(io, "r: $(ins.r)\n")
    write(io, "e: $(ins.e)\n")
    write(io, "l: $(ins.l)\n")
    write(io, "mind: $(ins.mind)\n")
    write(io, "maxd: $(ins.maxd)\n")
    write(io, "d: $(ins.d)\n")
    write(io, "PRE: $(ins.PRE)\n")
    write(io, "SYN: $(ins.SYN)\n")
    close(io)
    println("........Completed........")
end