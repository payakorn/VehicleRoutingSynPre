const ins_names = [
    "ins$i-$j" for i in [10, 25, 50, 75, 100, 200, 300] for j in 1:10
]


function create_dict_for_opt(model)
    nothing
end

function find_opt(ins_name...)
    for ins in ins_name
        find_opt(ins)
    end
end


function find_opt()
    find_opt("ins10-1")
end


function find_opt(number_of_node::Integer)
    for ins in ("ins$number_of_node-$j" for j in 1:10)
        @info "finding optimal solution for instance: $ins"
        find_opt(ins)
    end
end


function find_opt(ins_name; timelim=3600)

    location = dir("data", "HHCRSP_opt", "$ins_name.yml")

    if isfile(location)

        data = YAML.load_file(location)

        if data["termination_status"] != "OPTIMAL" && data["solve_time"] <= timelim
            nothing
        else
            return nothing
        end
    end

    @info "optimizing $(ins_name)"
    # @load "data/raw_HHCRSP/ins10-1.jld2"

    hello_dict = load("data/raw_HHCRSP/$ins_name.jld2")

    e = hello_dict["e"]
    r = hello_dict["r"]
    num_node = hello_dict["num_node"]
    a = hello_dict["a"]
    DS = hello_dict["DS"]
    d = hello_dict["d"]
    mind = hello_dict["mind"]
    maxd = hello_dict["maxd"]
    xx = hello_dict["xx"]
    l = hello_dict["l"]
    num_serv = hello_dict["num_serv"]
    num_vehi = hello_dict["num_vehi"]
    yy = hello_dict["yy"]
    p = hello_dict["p"]

    # a = ones(num_vehi, num_serv)
    # r = ones(num_node, num_serv)

    # load parameters
    # num_node = 11
    # num_vehi = 3
    # num_serv = 6
    M = num_node * 1000

    # create set of indices
    N = 1:(num_node)
    N_c = 2:(num_node)
    K = 1:num_vehi
    S = 1:num_serv

    # generate set of i, j in N with i != j
    IJ = Iterators.filter(x -> x[1] != x[2], Iterators.product(N, N))
    SS = Iterators.filter(x -> x[1] != x[2], Iterators.product(S, S))
    KK = Iterators.filter(x -> x[1] != x[2], Iterators.product(K, K))

    r_syn = deepcopy(r)
    Q = []
    SYN = []
    SYN_num = ones(num_node, num_serv)
    for i in N_c
        pos = findall(x -> x == 1.0, r[i, :])
        if mind[i] == 0 && maxd[i] == 0
            push!(Q, (i, length(pos) - 1))
            push!(SYN, (i, pos[1]))
            r_syn[i, pos[1]] = 2.0
            r_syn[i, pos[2]] = 0.0
            r[i, pos[2]] = 0.0
            SYN_num[i, pos[1]] = 2.0
        else
            push!(Q, (i, length(pos)))
        end
    end

    # create PRE set
    PRE = []
    PRE_node = []
    for i in N_c
        fx = findall(x -> x == 1, r[i, :])
        if length(fx) > 1
            for j in 2:length(fx)
                push!(PRE, (i, fx[j-1], fx[j], mind[i], maxd[i]))
                push!(PRE_node, i)
            end
        end
    end

    # for (i, x1, x2, min_d, max_d) in PRE
    #     Iterators.filter!(x -> (x[1], x[2]) == (x1, x2), SS)
    # end

    # model
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", timelim)

    # variables
    @variable(model, x[i=N, j=N, k=K; i != j], Bin)
    # @variable(model, e[i]<=t[i=N, k=K]<=l[i])
    @variable(model, t[i=N, k=K] >= e[i])
    @variable(model, ts[i=N, k=K, s=S] >= 0)
    @variable(model, y[i=N_c, k=K, s=S], Bin)
    # @variable(model, z[j=N_c, s1=S, s2=S]<=r[j, s2], Bin) # s1 is position of job, s2 is service
    @variable(model, z[j=N_c, s1=S, s2=S], Bin) # s1 is position of job, s2 is service
    @variable(model, zz[i=N_c, k=K, s=S] >= 0)
    @variable(model, Tmax >= 0)
    # constraints

    # 1
    for (i, j) in IJ
        if j > 1
            for k in K
                @constraint(model, x[i, j, k] <= sum(y[j, k, s] for s in S))
            end
        end
    end

    for k in K
        for j in N_c
            @constraint(model, sum(y[j, k, s] for s in S) <= M * sum(x[i, j, k] for i in N if i != j))
        end
    end

    # 2, 3
    for k in K
        @constraint(model, sum(x[1, j, k] for j in N if 1 != j) == 1)
        @constraint(model, sum(x[i, 1, k] for i in N if 1 != i) == 1)
    end

    # 4
    for k in K
        for j in N_c
            @constraint(model, sum(x[i, j, k] for i in N if i != j) - sum(x[j, l, k] for l in N if j != l) == 0.0)
        end
    end

    # subtour
    # @variable(model, 1 <= u[i=N_c, k=K] <= num_node-1)
    # for (i, j) in IJ
    #     if i > 1 && j > 1
    #         for k in K
    #             @constraint(model, u[i, k] - u[j, k] + 1 <= M*(1 - x[i, j, k]))
    #         end
    #     end
    # end

    # 5
    for s in S
        for j in N_c
            if r_syn[j, s] != 0.0
                @constraint(model, sum(a[k, s] * y[j, k, s] for k in K) == r_syn[j, s])
            else
                @constraint(model, sum(y[j, k, s] for k in K) == 0)
            end

            # @constraint(model, sum(y[j, k, s] for k in K) == r[j, s])
        end
    end

    # 6
    for k in K
        # fix(t[0,k], 0, force=true)
        for j in N_c
            @constraint(model, d[1, j] <= t[j, k] + M * (1 - x[1, j, k]))
        end
    end

    for (i, j) in IJ
        if i > 1
            for k in K
                # @constraint(model, t[i, k] + sum(p[k, s, i]*y[i, k, s] for s in S) + d[i, j] - M*(1-x[i, j, k]) <= t[j, k])
                for s in S
                    @constraint(model, ts[i, k, s] + p[k, s, i] + d[i, j] - M * (1 - x[i, j, k]) <= t[j, k])
                end
            end
        end
    end

    # for j in N
    #     for k in K
    #         @constraint(model, e[j]*sum(x[i, j, k] for i in N if i != j) <= t[j, k])
    #         @constraint(model, l[j]*sum(x[i, j, k] for i in N if i != j) >= t[j, k])
    #     end
    # end

    for j in N_c
        # test
        for k in K
            for s in S
                @constraint(model, e[j] * y[j, k, s] <= ts[j, k, s])
                @constraint(model, l[j] * y[j, k, s] >= ts[j, k, s] - zz[j, k, s])
                @constraint(model, -M * y[j, k, s] <= ts[j, k, s])
                @constraint(model, M * y[j, k, s] >= ts[j, k, s])
            end
        end
    end

    # Tmax
    for i in N_c
        for s in S
            for k in K
                @constraint(model, Tmax >= zz[i, k, s])
            end
        end
    end

    # 7
    for (i, s1, s2, min_d, max_d) in PRE
        # @constraint(model, sum(ts[i, k, s1] for k in K) + sum(p[k, s1, i]*y[i, k, s1] for k in K) <= sum(ts[i, k, s2] for k in K) + M*(2-sum(y[i, k, s1] for k in K)-sum(y[i, k, s2] for k in K)))
        @constraint(model, sum(ts[i, k, s1] for k in K) + min_d <= sum(ts[i, k, s2] for k in K) + M * (2 - sum(y[i, k, s1] for k in K) - sum(y[i, k, s2] for k in K)))
        @constraint(model, sum(ts[i, k, s2] for k in K) - max_d <= sum(ts[i, k, s1] for k in K) + M * (2 - sum(y[i, k, s1] for k in K) - sum(y[i, k, s2] for k in K)))
    end

    for i in N_c
        for k in K
            for s in S
                @constraint(model, t[i, k] <= ts[i, k, s] + M * (1 - y[i, k, s]))
            end
        end
    end

    # Synchronization
    # SYN = [(11, 3, 1, 2)]
    # for (i, s) in SYN
    #     @constraint(model, sum(ts[i, k1, s]) == sum(ts[i, k2, s]))
    # end
    # @constraint(model, ts[11, 1, 3] == ts[11, 2, 3])

    for (i, s) in SYN
        for (k1, k2) in KK
            @constraint(model, -M * (2 - y[i, k1, s] - y[i, k2, s]) <= ts[i, k1, s] - ts[i, k2, s])
            @constraint(model, ts[i, k1, s] - ts[i, k2, s] <= M * (2 - y[i, k1, s] - y[i, k2, s]))
        end
    end
    # new constraints z positions ()

    # for j in N_c
    #     position = findall(x -> x != 1.0, r[j, :])
    #     for i in 1:num_serv-length(position)
    #         for l in position
    #             fix(z[j, i, l], 0, force=true)
    #         end
    #     end

    #     for i in num_serv-length(position)+1:num_serv
    #         for l in S
    #             fix(z[j, i, l], 0, force=true)
    #         end
    #     end

    # end

    for (j, num_q) in Q
        for i in 1:num_q
            @constraint(model, sum(z[j, i, s2] for s2 in S) == 1)
        end
        for s in S
            @constraint(model, sum(z[j, i, s] for i in 1:num_q) == r[j, s])
        end
    end

    for (s1, s2) in SS
        for s in setdiff(S, 1) # position from 2 to S
            for j in setdiff(N_c, PRE_node)
                # @constraint(model, sum(ts[j, k, s1] for k in K) + p[2, s1, j] - M*(2 - z[j, s-1, s1] - z[j, s, s2]) <= sum(ts[j, k, s2] for k in K))
                @constraint(model, sum(ts[j, k, s1] for k in K) / SYN_num[j, s1] + p[2, s1, j] - M * (2 - z[j, s-1, s1] - z[j, s, s2]) <= sum(ts[j, k, s2] for k in K) / SYN_num[j, s2])
            end
        end
    end

    # objective 
    @objective(model, Min, 1 / 3 * sum(d[i, j] * x[i, j, k] for i in N for j in N for k in K if i != j) + 1 / 3 * Tmax + 1 / 3 * sum(zz[i, k, s] for i in N_c for k in K for s in S))

    # optimize
    optimize!(model)
    println(solution_summary(model; verbose=true))
    route = Dict()
    starttime = Dict()
    late = Dict()
    num_job = Dict()
    if model |> has_values
        for k in K
            for i in N
                for j in N
                    if i != j
                        if value.(x[i, j, k]) == 1.0
                            println("x($i, $j, $k)")
                        end
                    end
                end
            end
        end

        # for i in N
        #     for k in K
        #         for s in S
        #             if r[i, s] == 1
        #                 println("k=$k, r[$i, $s], y: $(value.(y[i, k, s]))")
        #             end
        #         end
        #     end
        # end

        # function print_value(X)
        #     for x in value.(X)
        #         println(values.(x))
        #     end
        # end

        for k in K
            route[k] = [1]
            starttime[k] = [0.0]
            late[k] = [0.0]
            num_job[k] = [0]

            job = 1
            for j in N_c
                if abs(value.(x[1, j, k]) - 1.0) <= 1e-6
                    job = deepcopy(j)
                    push!(route[k], job)
                    push!(starttime[k], value.(t[j, k]))
                    push!(late[k], l[j] - value.(t[j, k]))
                    push!(num_job[k], sum([value.(y[job, k, s]) for s in S]))
                    break
                end
            end

            iter = 1
            while job != 1 && iter <= num_node - 1
                iter += 1
                for j in setdiff(N, job)
                    if abs(value.(x[job, j, k]) - 1.0) <= 1e-20
                        job = deepcopy(j)
                        push!(route[k], job)
                        push!(starttime[k], value.(t[j, k]))
                        push!(late[k], l[j] - value.(t[j, k]))
                        if job != 1
                            push!(num_job[k], sum([value.(y[job, k, s]) for s in S]))
                        end
                        break
                    end
                end
            end
        end

        resultDict = Dict(
            "name" => ins_name,
            "route" => route,
            "starttime" => starttime,
            "num_job" => num_job,
            "solver" => JuMP.solver_name(model),
            "late" => late,
            "Tmax" => value.(Tmax),
            "solve_time" => JuMP.solve_time(model),
            "obj_value" => JuMP.objective_value(model),
            "relative_gap" => JuMP.relative_gap(model),
            "output_text" => "$(solution_summary(model))",
            "termination_status" => termination_status(model),
        )
    else

        resultDict = Dict(
            "name" => ins_name,
            "route" => route,
            "starttime" => starttime,
            "num_job" => num_job,
            "solver" => JuMP.solver_name(model),
            "late" => late,
            "Tmax" => Inf,
            "solve_time" => JuMP.solve_time(model),
            "obj_value" => Inf,
            "relative_gap" => 1,
            "output_text" => "$(solution_summary(model))",
            "termination_status" => termination_status(model),
        )
    end

    # the example to use the dictionary for the variables
    # resultDict = Dict(k => !(value.(v) isa JuMP.Containers.DenseAxisArray) ? value.(v) : value.(v) |> Array for (k, v) in object_dictionary(model) if v isa AbstractArray{VariableRef})
    # resultDict = Dict(k => value.(v) for (k, v) in object_dictionary(model) if v isa AbstractArray{VariableRef})
    save_opt_solution_to_YAML(resultDict)

    return resultDict
end

function save_opt_solution_to_YAML(results::Dict)
    YAML.write_file(dir("data", "HHCRSP_opt", "$( results["name"] ).yml"), results)
end

function load_opt_solution_YAML(ins_name::String)
    location = dir("data", "HHCRSP_opt", "$ins_name.yml")
    return YAML.load_file(location)
end
