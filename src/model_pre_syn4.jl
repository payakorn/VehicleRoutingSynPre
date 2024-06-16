using JuMP, Gurobi, JLD2

@load "data/raw_HHCRSP/ins10-2.jld2"

# a = ones(num_vehi, num_serv)
# r = ones(num_node, num_serv)

# load parameters
# num_node = 11
# num_vehi = 3
# num_serv = 6
M = num_node*1000
mind = zeros(num_node)
maxd = M*ones(num_node)

# create set of indices
N = 1:(num_node)
N_c = 2:(num_node)
K = 1:num_vehi
S = 1:num_serv

# generate set of i, j in N with i != j
IJ = Iterators.filter(x -> x[1] != x[2], Iterators.product(N, N))
SS = Iterators.filter(x -> x[1] != x[2], Iterators.product(S, S))
KK = Iterators.filter(x -> x[1] != x[2], Iterators.product(K, K))

# create PRE set (include Synchronization)
PRE = Dict()
for i in N_c
    PRE[i] = Int64[]
    xx = findall(x -> x == 1, r[i, :])
    append!(PRE[i], xx)
end

# model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 12000)
# set_optimizer_attribute(model, "Presolve", 0) # for Gurobi

# variables
@variable(model, x[i=N, j=N, k=K; i!=j], Bin)
# @variable(model, e[i]<=t[i=N, k=K]<=l[i])
@variable(model, t[i=N, k=K] >= e[i])
@variable(model, ts[i=N, k=K, s=S] >= 0)
@variable(model, y[i=N_c, k=K, s=S], Bin)
@variable(model, zz[i=N_c, s=S] >= 0)
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
        @constraint(model, sum(y[j, k, s] for s in S) <= M*sum(x[i, j, k] for i in N if i != j))
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

# 5
for s in S
    for j in N_c
        if r[j, s] != 0.0
            @constraint(model, sum(a[k, s]*y[j, k, s] for k in K) == r[j, s])
        else
            @constraint(model, sum(y[j, k, s] for k in K) == 0)
        end
    end
end

# 6
for k in K
    for j in N_c
        @constraint(model, d[1, j] <= t[j, k]+ M*(1-x[1, j, k]))
    end
end

for (i, j) in IJ
    if i > 1
        for k in K
            # @constraint(model, t[i, k] + sum(p[k, s, i]*y[i, k, s] for s in S) + d[i, j] - M*(1-x[i, j, k]) <= t[j, k])
            for s in S
                @constraint(model, ts[i, k, s] + p[k, s, i] + d[i, j] - M*(1-x[i, j, k]) <= t[j, k])
            end
        end
    end
end

for j in N_c
    # test
    for k in K
        for s in S
            # @constraint(model, e[j]*y[j, k, s] <= ts[j, k, s])
            @constraint(model, l[j]*y[j, k, s] >= ts[j, k, s] - zz[j, s])
            @constraint(model, -M*y[j, k, s] <= ts[j, k, s])
            @constraint(model, M*y[j, k, s] >= ts[j, k, s])
        end
    end
end

# Tmax
for i in N_c
    for s in S
        for k in K
            @constraint(model, Tmax >= zz[i, s])
        end
    end
end

# 7
# for (i, s1, s2, min_d, max_d) in PRE
for i in N_c
    len_pre = length(PRE[i])
    if len_pre > 1
        for s in 1:(len_pre-1)
            (s1, s2) = (PRE[i][s], PRE[i][s+1])
            if mind[i] == 0 && maxd[i] == 0
                for k in K
                    @constraint(model, sum(y[i, k, s] for s in [s1, s2]) <= 1)
                end
            end
            @constraint(model, sum(ts[i, k, s1] for k in K) + mind[i] <= sum(ts[i, k, s2] for k in K) + M*(2-sum(y[i, k, s1] for k in K)-sum(y[i, k, s2] for k in K)))
            @constraint(model, sum(ts[i, k, s2] for k in K) - maxd[i] <= sum(ts[i, k, s1] for k in K) + M*(2-sum(y[i, k, s1] for k in K)-sum(y[i, k, s2] for k in K)))
        end
    end
end


for i in N_c
    for k in K
        for s in S
            @constraint(model, t[i, k] <= ts[i, k, s] + M*(1-y[i, k, s]))
        end
    end
end

# Synchronization
# SYN = [(11, 3, 1, 2)]
# for (i, s) in SYN
#     @constraint(model, sum(ts[i, k1, s]) == sum(ts[i, k2, s]))
# end
# @constraint(model, ts[11, 1, 3] == ts[11, 2, 3])

# for (i, s) in SYN
#     for (k1, k2) in KK
#         @constraint(model, -M*(2-y[i, k1, s]-y[i, k2, s]) <= ts[i, k1, s] - ts[i, k2, s])
#         @constraint(model, ts[i, k1, s] - ts[i, k2, s] <= M*(2-y[i, k1, s]-y[i, k2, s]))
#     end
# end


# for (s1, s2) in SS
#     for s in setdiff(S, 1)
#         for j in N_c
#             # @constraint(model, sum(ts[j, k, s1] for k in K) + p[2, s1, j] - M*(2 - z[j, s-1, s1] - z[j, s, s2]) <= sum(ts[j, k, s2] for k in K))
#             @constraint(model, sum(ts[j, k, s1] for k in K)/SYN_num[j, s1] + p[2, s1, j] - M*(2 - z[j, s-1, s1] - z[j, s, s2]) <= sum(ts[j, k, s2] for k in K)/SYN_num[j, s2])
#         end
#     end
# end

# objective function
@objective(model, Min, 1/3*sum(d[i, j]*x[i, j, k] for i in N for j in N for k in K if i != j) + 1/3*Tmax + 1/3*sum(zz[i, s] for i in N_c for s in S))


# optimize
optimize!(model)


route = Dict()
starttime = Dict()
late = Dict()
num_job = Dict()
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
    while job != 1 && iter <= num_node-1
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

resultDict = Dict(k => value.(v) for (k, v) in object_dictionary(model))