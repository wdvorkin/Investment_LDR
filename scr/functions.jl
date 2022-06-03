## auxiliary functions used in data manipulation and model specification
ns(l) = Int(net[:n_s][l])
nr(l) = Int(net[:n_r][l])
Y_(t,w,h,dict) = hcat(dict["Y"][t][w][h]...)
P_(t,w,h,dict) = hcat(dict["P"][t][w][h]...)
Y̅_(t,dict) = hcat(dict["Y̅"][t]...)
Θ̅_(t,dict) = hcat(dict["Θ̅"][t]...)
Φ̅_(t,dict) = hcat(dict["Φ̅"][t]...)
function  Ψ(x)
    if set[:cc_ref] == "Normal"
        return quantile(Normal(0,1),1-x)
    else
        return sqrt((1-x)/x)
    end
end
chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]
function my_chunk(arr,n)
    h = []
    for w in 1:set[:W]
        push!(h,chunk(arr, 168)[w][length(chunk(arr, 168)[w])-(set[:H]-1):length(chunk(arr, 168)[w])])
    end
    return h
end
ind_can(string_) = findall(x->x==true, [occursin(string_, data[:gen_c][:tech][i]) for i in 1:data[:dims][:J]])
ind_exs(string_) = findall(x->x==true, [occursin(string_, data[:gen_e][:tech][i]) for i in 1:data[:dims][:G]])
function remove_col_and_row(B,refbus)
    @assert size(B,1) == size(B,2)
    n = size(B,1)
    return B[1:n .!= refbus, 1:n .!= refbus]
end
function build_B̆(B̂inv,refbus)
    Nb = size(B̂inv,1)+1
    B̆ = zeros(Nb,Nb)
    for i in 1:Nb, j in 1:Nb
        if i < refbus && j < refbus
            B̆[i,j] = B̂inv[i,j]
        end
        if i > refbus && j > refbus
            B̆[i,j] = B̂inv[i-1,j-1]
        end
        if i > refbus && j < refbus
            B̆[i,j] = B̂inv[i-1,j]
        end
        if i < refbus && j > refbus
            B̆[i,j] = B̂inv[i,j-1]
        end
    end
    return B̆
end
function partioning(br_num,br,N_s)
    x = 1:N_s
    reshape(x, Int(N_s/br_num),br_num)[:,br]
end
function gen_rand_ξ(t)
    ξ̂ = ones(data[:n_t][data[:T]],1)
    ξ̂[2:end,:] = rand(MvNormal(ξ[:μ],ξ[:Σ]),1)
    return S(t)*ξ̂
end
function scenario_tree(data)
    N_s = (data[:BF]^(data[:T]-1))
    Tree_Load = zeros(data[:D],N_s,data[:T])
    Tree_CAPEX = zeros(data[:J],N_s,data[:T])
    Tree_OPEX_ex = zeros(data[:G],N_s,data[:T])
    Tree_OPEX_ca = zeros(data[:J],N_s,data[:T])
    for t in 1:data[:T]
        if t == 1
            Tree_Load[:,1:N_s,t] .= L(1)
            Tree_CAPEX[:,1:N_s,t] .= Q_g(1)
            Tree_OPEX_ex[:,1:N_s,t] .= V(1)
            Tree_OPEX_ca[:,1:N_s,t] .= B(1)
        end
        # loads
        if t >= 2
            br_num = data[:BF]^(t-1)
            for br in 1:br_num
                ξ̂ = gen_rand_ξ(t)
                Tree_Load[:,partioning(br_num,br,N_s),t] .= L(t)*ξ̂
                Tree_CAPEX[:,partioning(br_num,br,N_s),t] .= Q_g(t)*ξ̂
            end
        end
        # capex
        if t >= 2
            br_num = data[:BF]^(t-1)
            for br in 1:br_num
                ξ̂ = gen_rand_ξ(t)
                Tree_CAPEX[:,partioning(br_num,br,N_s),t] .= Q_g(t)*ξ̂
            end
        end
        # opex
        if t >= 2
            br_num = data[:BF]^(t-1)
            for br in 1:br_num
                ξ̂ = gen_rand_ξ(t)
                Tree_OPEX_ex[:,partioning(br_num,br,N_s),t] .= V(t)*ξ̂
                Tree_OPEX_ca[:,partioning(br_num,br,N_s),t] .= B(t)*ξ̂
            end
        end
    end
    return Dict(:L => Tree_Load, :Q_g => Tree_CAPEX, :V => Tree_OPEX_ex, :B =>Tree_OPEX_ca )
end
function uniform_sample(n,set,μ)
    a = zeros(n); b = zeros(n)
    a = μ .- sqrt(12)/2*sqrt(set[:σ])
    b = sqrt(12)/2*sqrt(set[:σ]) .+ μ
    ξ̂ = ones(n,set[:Ns])
    for i in 1:n-1
        ξ̂[i+1,:] = rand(Uniform(a[i],b[i]),set[:Ns])
    end
    return ξ̂
end
function laplace_sample(n,set,μ,Σ)
    b = zeros(n-1)
    for i in 1:n-1
        b[i] = sqrt(Σ[i,i]/2)
    end
    ξ̂ = ones(n,set[:Ns])
    for i in 1:n-1
        ξ̂[i+1,:] = rand(Laplace(μ[i],b[i]),set[:Ns])
    end
    return ξ̂
end
function logistic_sample(n,set,μ,Σ)
    θ = zeros(n-1)
    for i in 1:n-1
        θ[i] = sqrt(3*Σ[i,i])/π
    end
    ξ̂ = ones(n,set[:Ns])
    for i in 1:n-1
        ξ̂[i+1,:] = rand(Logistic(μ[i],θ[i]),set[:Ns])
    end
    return ξ̂
end

## function to extract data
function load_data(caseID,set)

    tech_data    = CSV.read("$(caseID)/technology.csv"     ,        DataFrame; header=1, skipto = 3)

    gen_e_data   = CSV.read("$(caseID)/existing_units.csv",         DataFrame; header=1, skipto = 3)

    gen_c_data   = CSV.read("$(caseID)/candidate_units.csv",        DataFrame; header=1, skipto = 3)

    # power_factor = CSV.read("$(caseID)/renewable_power_factor.csv", DataFrame; header=1, skipto = 2, limit = set[:H])
    power_factor = CSV.read("$(caseID)/renewable_power_factor.csv", DataFrame; header=1, skipto = 2)

    trans_data   = CSV.read("$(caseID)/transmission.csv",           DataFrame; header=1)

    demand_data  = CSV.read("$(caseID)/load.csv",           DataFrame; header=1)

    # load_factor  = CSV.read("$(caseID)/load_factor.csv", DataFrame; header=1, skipto = 2, limit = set[:H])
    load_factor  = CSV.read("$(caseID)/load_factor.csv", DataFrame; header=1, skipto = 2)

    weights      = CSV.read("$(caseID)/oper_cond_weights.csv",           DataFrame; header=0, limit = set[:W])

    retirement   = CSV.read("$(caseID)/retirement.csv",           DataFrame; header=1)


    # problem dimentions
    N = maximum(gen_e_data[!,:node])
    E = size(trans_data,1)
    G = size(gen_e_data,1)
    J = size(gen_c_data,1)
    D = size(demand_data,1)
    W = set[:W]; H = set[:H]
    K = 4
    # save to dictionary
    dims = Dict(:N => N,:E => E,:G =>G ,:J => J,:D => D, :K => K)


    # Load existing generation data
    p̅ = zeros(G,set[:T]); p̲ = zeros(G); r̅ = zeros(G); # capacity limits
    cost_fuel = zeros(G); cost_var = zeros(G); cost_fix = zeros(G); # cost data
    M = zeros(N,G); k = ones(G,W,H) # location and power factor data
    renew_flag = zeros(G); em_rate = zeros(G); h_rate = zeros(G); tech = []# technology data data
    for i in 1:G
        for t in 1:set[:T]
            p̅[i,t] = round(gen_e_data[i,:p_max] .* retirement[i,t+2],digits=2)
        end
        p̲[i] = gen_e_data[i,:p_min]
        r̅[i]         = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:ramping_capacity]
        cost_fuel[i] = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:fuel_price]
        cost_var[i]  = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:cost_variable]
        cost_fix[i]  = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:cost_fixed]
        em_rate[i]   = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:emission_rate]
        h_rate[i]    = tech_data[findall(tech_data.type .== gen_e_data[i,:technology])[1],:heat_rate]
        M[gen_e_data[i,:node],i] = 1
        push!(tech,gen_e_data[i,:technology])
        (occursin("Wind", gen_e_data[i,:technology]) || occursin("PV", gen_e_data[i,:technology]) || occursin("Hydro", gen_e_data[i,:technology])) ? renew_flag[i] = 1 : NaN
        if renew_flag[i] == 1
            for w in 1:W
                # h = vec([H*(week-1)+1:H*week for week in 1:W])
                h = my_chunk(1:2352,set[:H])
                k[i,w,:] = power_factor[h[w],gen_e_data[i,:technology]]
            end
        end
    end

    # save to dictionary
    gen_e = Dict(:p̅ => p̅, :p̲ => p̲, :r̅ => r̅, :k => k,
    :c_fuel => cost_fuel, :c_var => cost_var, :c_fix => cost_fix,
    :e => em_rate, :h => h_rate, :M => M, :renew_flag => renew_flag, :tech => tech)

    # Load candidate generation data
    y̅ = zeros(J); y̲ = zeros(J); r̅ = zeros(J); # capacity limits
    cost_fuel = zeros(J); cost_var = zeros(J); cost_fix = zeros(J); cost_inv = zeros(J); # cost data
    M = zeros(N,J); k = ones(J,W,H); # location and power factor data
    renew_flag = zeros(J); em_rate = zeros(J); h_rate = zeros(J); life_time = zeros(J); tech = [] # technology data
    for i in 1:J
        y̅[i] = gen_c_data[i,:y_max] ; y̲[i] = gen_c_data[i,:y_min]
        r̅[i]         = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:ramping_capacity]
        cost_inv[i]  = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:cost_inv_overnight]
        cost_fuel[i] = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:fuel_price]
        cost_var[i]  = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:cost_variable]
        cost_fix[i]  = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:cost_fixed]
        em_rate[i]   = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:emission_rate]
        h_rate[i]    = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:heat_rate]
        life_time[i] = tech_data[findall(tech_data.type .== gen_c_data[i,:technology])[1],:life_time]
        M[gen_c_data[i,:node],i] = 1
        push!(tech,gen_c_data[i,:technology])
        (occursin("Wind", gen_c_data[i,:technology]) || occursin("PV", gen_c_data[i,:technology])) ? renew_flag[i] = 1 : NaN
        if renew_flag[i] == 1
            for w in 1:W
                # h = vec([H*(week-1)+1:H*week for week in 1:W])
                h = my_chunk(1:2352,set[:H])
                k[i,w,:] = power_factor[h[w],gen_c_data[i,:technology]]
            end
        end
    end
    # save to dictionary
    gen_c = Dict(:y̅ => y̅, :y̲ => y̲, :r̅ => r̅, :k => k,
    :y̅_max => [y̅[i] == -1 ? 500000.0 : y̅[i]  for i in 1:J],
    :y̅_pv => [occursin("PV", tech[j]) == true ? 1 : 0 for j in 1:J],
    :y̅_wd => [occursin("Wind", tech[j]) == true ? 1 : 0 for j in 1:J],
    :c_fuel => cost_fuel, :c_var => cost_var, :c_fix => cost_fix, :c_inv => cost_inv,
    :e => em_rate, :h => h_rate, :M => M, :renew_flag => renew_flag, :life_time => life_time, :tech => tech)

    # Load network data
    β = trans_data[!,:beta]
    f̅ = trans_data[!,:f_max]
    n_s = trans_data[!,:n_s]
    n_r = trans_data[!,:n_r]
    # Compute PTDF matrix
    B_line = zeros(E,N); B̃_bus = zeros(N,N); B_ = zeros(N,N)
    for n in 1:N
        for l in 1:E
            if n_s[l] == n
                B_[n,n] += β[l]
                B_line[l,n] = β[l]
            end
            if n_r[l] == n
                B_[n,n] += β[l]
                B_line[l,n] = -β[l]
            end
        end
    end
    for l in 1:E
        B_[Int(n_s[l]),Int(n_r[l])] = - β[l]
        B_[Int(n_r[l]),Int(n_s[l])] = - β[l]
    end
    ref_node = 3
    B̃_bus = remove_col_and_row(B_,ref_node)
    B̃_bus = inv(B̃_bus)
    B̃_bus = build_B̆(B̃_bus,ref_node)
    PTDF = round.(B_line * B̃_bus, digits=5)
    # save to dictionary
    net = Dict(:f̅ => f̅, :n_s => n_s, :n_r => n_r, :F => PTDF)

    # Load demand data
    l = demand_data[!,:max_load]
    M = zeros(N,D); k = zeros(D,set[:T],W,H);
    for i in 1:D
        for w in 1:W
            # h = vec([H*(week-1)+1:H*week for week in 1:W])
            h = my_chunk(1:2352,set[:H])
            for t in 1:set[:T]
                load_factor__  = CSV.read("$(caseID)/load_factor_$(t).csv", DataFrame; header=1, skipto = 2)
                k[i,t,w,:] = load_factor__[h[w],"$(demand_data[i,:index])"]
            end
        end
        for n in 1:N
            demand_data[!,:node][i] == n ? M[n,i] = 1 : NaN
        end
    end
    # save to dictionary
    load = Dict(:l => l, :M => M, :k => k, :w => weights[!,:Column1]./168)

    # storage data
    M = zeros(N,K); ϑ_o = zeros(K); cost_inv_energy = zeros(K); cost_inv_charge = zeros(K);
    M = diagm(ones(4)); cost_inv_energy[:] .= 2*190000; cost_inv_charge[:] .= 198000; LT = 15 .*ones(K);
    η⁻ = 0.95; η⁺ = 0.95
    # save to dictionary
    stor = Dict(:M => M, :ϑ_o => ϑ_o, :cost_inv_energy => cost_inv_energy, :cost_inv_charge => cost_inv_charge, :life_time => LT, :η⁻ => 0.95, :η⁺ => 0.95)

    return Dict(:gen_e => gen_e, :gen_c => gen_c, :net => net, :load => load, :stor => stor, :dims => dims)
end

## functions to model uncertain planning data
function V(i)
    V = zeros(data[:dims][:G],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:G]
        V[j,1] = (data[:gen_e][:c_fuel] .* data[:gen_e][:h] .+ data[:gen_e][:c_var])[j]
        if occursin("Coal", data[:gen_e][:tech][j]) == true
            t > 1 ? V[j,set[:n_t][t]-1] = V[j,1]*set[:r_cl][t] : NaN
        else
            t > 1 ? V[j,set[:n_t][t]-1] = V[j,1]*set[:r_ng][t] : NaN
        end
    end
    return V
end
function B(i)
    B = zeros(data[:dims][:J],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:J]
        B[j,1] = (data[:gen_c][:c_fuel] .* data[:gen_c][:h] .+ data[:gen_c][:c_var])[j]
        if occursin("Coal", data[:gen_c][:tech][j]) == true
            t > 1 ? B[j,set[:n_t][t]-1] = B[j,1]*set[:r_cl][t] : NaN
        else
            t > 1 ? B[j,set[:n_t][t]-1] = B[j,1]*set[:r_ng][t] : NaN
        end
    end
    return B
end
function L(i)
    L = zeros(data[:dims][:D],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:D]
        L[j,1] = data[:load][:l][j]
        t > 1 ? L[j,set[:n_t][t]] = L[j,1]*set[:r_l][t] : NaN
    end
    return L
end
function E̅(i)
    ### Option 0
    E = zeros(1,set[:n_t][i])
    E[1,1] = set[:e̅][i]
    return E
end
function Q(i)
    data_ = deepcopy(data)
    AIC = data_[:gen_c][:c_inv].*set[:int_rate].*(1+set[:int_rate]).^(data_[:gen_c][:life_time])./((1+set[:int_rate]).^(data_[:gen_c][:life_time]) .- 1)
    M = zeros(data[:dims][:J])
    C = zeros(data[:dims][:J])
    for j in 1:data[:dims][:J]
        M[j] = min(data_[:gen_c][:life_time][j],set[:year_per_period]*(set[:T]-i+1))
    end
    for j in 1:data[:dims][:J]
        C[j] = sum(AIC[j]/(1+set[:int_rate])^(k-1) for k in 1:M[j])
    end
    # build matrix Q
    Q_ = zeros(data[:dims][:J],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:J]
         Q_[j,1] = C[j]
         if occursin("Wind", data[:gen_c][:tech][j]) == true
             t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_wind][t] : NaN
         elseif occursin("PV", data[:gen_c][:tech][j]) == true
             t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_pv][t] : NaN
         elseif occursin("Nuclear", data[:gen_c][:tech][j]) == true
             t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_nucl][t] : NaN
         elseif occursin("GT", data[:gen_c][:tech][j]) == true || occursin("NG", data[:gen_c][:tech][j]) == true
             t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_ng][t] : NaN
         end
    end
    return Q_
end
function Q_ϑ̅(i) # energy storage capacity
    data_ = deepcopy(data)
    AIC = data_[:stor][:cost_inv_energy].*set[:int_rate].*(1+set[:int_rate]).^(data_[:stor][:life_time])./((1+set[:int_rate]).^(data_[:stor][:life_time]) .- 1)
    M = zeros(data[:dims][:K])
    C = zeros(data[:dims][:K])
    for j in 1:data[:dims][:K]
        M[j] = min(data_[:stor][:life_time][j],set[:year_per_period]*(set[:T]-i+1))
    end
    for j in 1:data[:dims][:K]
        C[j] = sum(AIC[j]/(1+set[:int_rate])^(k-1) for k in 1:M[j])
    end
    # build matrix Q
    Q_ = zeros(data[:dims][:K],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:K]
         Q_[j,1] = C[j]
         t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_s_e][t] : NaN
    end
    return Q_
end
function Q_φ̅(i) # charging storage capacity
    data_ = deepcopy(data)
    AIC = data_[:stor][:cost_inv_charge].*set[:int_rate].*(1+set[:int_rate]).^(data_[:stor][:life_time])./((1+set[:int_rate]).^(data_[:stor][:life_time]) .- 1)
    M = zeros(data[:dims][:K])
    C = zeros(data[:dims][:K])
    for j in 1:data[:dims][:K]
        M[j] = min(data_[:stor][:life_time][j],set[:year_per_period]*(set[:T]-i+1))
    end
    for j in 1:data[:dims][:K]
        C[j] = sum(AIC[j]/(1+set[:int_rate])^(k-1) for k in 1:M[j])
    end
    # build matrix Q
    Q_ = zeros(data[:dims][:K],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:K]
         Q_[j,1] = C[j]
         t > 1 ? Q_[j,set[:n_t][t]-2] = Q_[j,1] * set[:r_q_s_p][t] : NaN
    end
    return Q_
end
function O_c(i)
    O = zeros(data[:dims][:J],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:J]
        O[j,1] = (data[:gen_c][:c_fix])[j]
        if occursin("Wind", data[:gen_c][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1]*set[:r_o_wind][t] : NaN
        elseif occursin("PV", data[:gen_c][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_pv][t] : NaN
        elseif occursin("Nuclear", data[:gen_c][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_nucl][t] : NaN
        elseif occursin("GT", data[:gen_c][:tech][j]) == true || occursin("NG", data[:gen_c][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_ng][t] : NaN
        end
    end
    return O
end
function O_e(i)
    O = zeros(data[:dims][:G],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:G]
        O[j,1] = (data[:gen_e][:c_fix])[j]
        if occursin("Wind", data[:gen_e][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1]*set[:r_o_wind][t] : NaN
        elseif occursin("PV", data[:gen_e][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_pv][t] : NaN
        elseif occursin("Nuclear", data[:gen_e][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_nucl][t] : NaN
        elseif occursin("GT", data[:gen_e][:tech][j]) == true || occursin("NG", data[:gen_e][:tech][j]) == true
            t > 1 ? O[j,set[:n_t][t]-2] = O[j,1] * set[:r_o_ng][t] : NaN
        end
    end
    return O
end
function R_e(i)
    data_ = deepcopy(data)
    c_inv = data_[:gen_c][:c_inv]
    R = zeros(data[:dims][:G],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:G]
        if occursin("Wind", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Wind")[1]] * 0.01
        end
        if occursin("PV", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("PV")[1]] * 0.01
        end
        if occursin("Nuclear", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Nuclear")[1]] * 0.01
        end
        if occursin("GT", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("GT")[1]] * 0.01
        end
        if occursin("NG", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("NG")[1]] * 0.01
        end
        if occursin("Coal", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("NG")[1]] * 0.01
        end
        if occursin("Hydro", data[:gen_e][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Nuclear")[1]] * 0.01
        end

        if occursin("Wind", data[:gen_e][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1]*set[:r_o_wind][t] : NaN
        elseif occursin("PV", data[:gen_e][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_pv][t] : NaN
        elseif occursin("Nuclear", data[:gen_e][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_nucl][t] : NaN
        elseif occursin("GT", data[:gen_e][:tech][j]) == true || occursin("NG", data[:gen_e][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_ng][t] : NaN
        end
    end
    return R
end
function R_c(i)
    data_ = deepcopy(data)
    c_inv = data_[:gen_c][:c_inv]
    R = zeros(data[:dims][:J],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:J]
        if occursin("Wind", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Wind")[1]] * 0.01
        end
        if occursin("PV", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("PV")[1]] * 0.01
        end
        if occursin("Nuclear", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Nuclear")[1]] * 0.01
        end
        if occursin("GT", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("GT")[1]] * 0.01
        end
        if occursin("NG", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("NG")[1]] * 0.01
        end
        if occursin("Coal", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("NG")[1]] * 0.01
        end
        if occursin("Hydro", data[:gen_c][:tech][j]) == true
            R[j,1] = c_inv[ind_can("Nuclear")[1]] * 0.01
        end

        if occursin("Wind", data[:gen_c][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1]*set[:r_o_wind][t] : NaN
        elseif occursin("PV", data[:gen_c][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_pv][t] : NaN
        elseif occursin("Nuclear", data[:gen_c][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_nucl][t] : NaN
        elseif occursin("GT", data[:gen_c][:tech][j]) == true || occursin("NG", data[:gen_c][:tech][j]) == true
            t > 1 ? R[j,set[:n_t][t]-2] = R[j,1] * set[:r_o_ng][t] : NaN
        end
    end
    return R
end
function O_ϑ̅(i) # fix cost of operation and maintaince for energy storage capacity
    O = zeros(data[:dims][:K],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:K]
        O[j,1] = (Q_ϑ̅(i)*S(i)*ξ[:μ̂])[j]*0.025*2/3 # 2230
        # (Q_ϑ̅(i)*S(i)*ξ[:μ̂])[j]*4/5
        t > 1 ? O[j,set[:n_t][t]-2] = O[j,1]*set[:r_o_s][t] : NaN
    end
    return O
end
function O_φ̅(i) # fix cost of operation and maintaince for energy storage dis/charging capacity
    O = zeros(data[:dims][:K],set[:n_t][i])
    for t in 1:i, j in 1:data[:dims][:K]
        O[j,1] = (Q_ϑ̅(i)*S(i)*ξ[:μ̂])[j]*0.025*1/3 # 750
        # (Q_φ̅(i)*S(i)*ξ[:μ̂])[j]*1/5
        t > 1 ? O[j,set[:n_t][t]-2] = O[j,1]*set[:r_o_s][t] : NaN
    end
    return O
end
function S(t)
    S = zeros(set[:n_t][t],set[:n_t][set[:T]])
    S[diagind(S)] .= 1
    return S
end

## function to generate planning uncertainty data (mean, covariance and samples)
function error_stat(set)
    n = set[:n_t][set[:T]]
    # mean and covariacne
    μ = ones(n-1)
    Σ = zeros(n-1,n-1)
    for i in 1:n-1
        Σ[i,i] = set[:σ]
    end
    # first two moments of altered dimentions
    μ̂ = vcat(1,μ)
    # E[ξξᵀ] = Cov + μμᵀ
    Σ̂ = zeros(n,n)
    Σ̂[1,1] = 1; Σ̂[2:end,1] .= μ; Σ̂[1,2:end] .= μ;
    Σ̂[2:end,2:end] .=  Σ .+ μ*μ'
    # factorization of E[ξξᵀ]
    Σ̅ = cholesky(Σ̂).L .+ zeros(n,n)
    Σ̅[:,1] .= 0
    ξ̂ = ones(n,set[:Ns])
    if set[:dist] == "Normal"
        ξ̂[2:end,:] = rand(MvNormal(μ,Σ),set[:Ns])
    elseif set[:dist] == "Uniform"
        ξ̂ = uniform_sample(n,set,μ)
    elseif set[:dist] == "Laplace"
        ξ̂ = laplace_sample(n,set,μ,Σ)
    elseif set[:dist] == "Logistic"
        ξ̂ = logistic_sample(n,set,μ,Σ)
    end

    return Dict(:μ => μ, :Σ => Σ, :μ̂ => μ̂, :Σ̂ => Σ̂, :Σ̅ => Σ̅, :ξ̂ => ξ̂)
end

## deterministic planning models
function determenistic_expansion_prim(set,data,ξ)
    # model definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # model variables
    @variable(model, ϑ̅[1:data[:dims][:K],1:set[:T]] >= 0)
    @variable(model, φ̅[1:data[:dims][:K],1:set[:T]] >= 0)
    @variable(model, y̅[1:data[:dims][:J],1:set[:T]] >= 0)
    @variable(model, ϑ[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, φ⁺[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, φ⁻[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, y[1:data[:dims][:J],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, p[1:data[:dims][:G],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    # model objective
    @objective(model, Min,
                    sum(# investment and fixed costs
                        + (Q(t)*S(t)*ξ[:μ̂] + sum(set[:year_per_period]*(O_c(τ)*S(τ)*ξ[:μ̂]) for τ in t:set[:T]))'*y̅[:,t]
                        + (Q_ϑ̅(t)*S(t)*ξ[:μ̂] + sum(set[:year_per_period]*(O_ϑ̅(τ)*S(τ)*ξ[:μ̂]) for τ in t:set[:T]))'*ϑ̅[:,t]
                        + (Q_φ̅(t)*S(t)*ξ[:μ̂] + sum(set[:year_per_period]*(O_φ̅(τ)*S(τ)*ξ[:μ̂]) for τ in t:set[:T]))'*φ̅[:,t]
                        # operating fuel costs
                        + set[:year_per_period]*52/sum(data[:load][:w])*sum(data[:load][:w][w]*
                        (
                          168/set[:H]*sum((B(t)*S(t)*ξ[:μ̂])'*y[:,w,h,t] for h in 1:set[:H]) +
                          168/set[:H]*sum((V(t)*S(t)*ξ[:μ̂])'*p[:,w,h,t] for h in 1:set[:H])
                        )
                    for w in 1:set[:W])
                    for t in 1:set[:T])
    )
    # power flow eqations
    @constraint(model, λ_b[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  ones(data[:dims][:N])'*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:μ̂])) .== 0)
    @constraint(model, λ_f̅[w=1:set[:W],h=1:set[:H],t=1:set[:T]], -data[:net][:F]*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:μ̂])) .>= -data[:net][:f̅])
    @constraint(model, λ_f̲[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  data[:net][:F]*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:μ̂])) .>= -data[:net][:f̅])
    # generation limits
    @constraint(model, λ_p[w=1:set[:W],h=1:set[:H],t=1:set[:T]], -p[:,w,h,t] .>= -data[:gen_e][:k][:,w,h] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_y[w=1:set[:W],h=1:set[:H],t=1:set[:T]], data[:gen_c][:k][:,w,h] .* sum(y̅[:,τ] for τ in 1:t) - y[:,w,h,t] .>= 0)
    @constraint(model, λ_r_e⁻[w=1:set[:W],h=2:set[:H],t=1:set[:T]], p[:,w,h,t] .- p[:,w,h-1,t] .>= - data[:gen_e][:r̅] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_r_e⁺[w=1:set[:W],h=2:set[:H],t=1:set[:T]], p[:,w,h-1,t] .- p[:,w,h,t] .>= - data[:gen_e][:r̅] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_r_c⁻[w=1:set[:W],h=2:set[:H],t=1:set[:T]], data[:gen_c][:r̅] .* sum(y̅[:,τ] for τ in 1:t) .+ y[:,w,h,t] .- y[:,w,h-1,t] .>= 0)
    @constraint(model, λ_r_c⁺[w=1:set[:W],h=2:set[:H],t=1:set[:T]], data[:gen_c][:r̅] .* sum(y̅[:,τ] for τ in 1:t) .+ y[:,w,h-1,t] .- y[:,w,h,t] .>= 0)
    # storage limits and equations
    @constraint(model, λ_s[w=1:set[:W],h=2:set[:H],t=1:set[:T]], ϑ[:,w,h,t] .- ϑ[:,w,h-1,t] .- φ⁺[:,w,h,t]*data[:stor][:η⁺] .+ φ⁻[:,w,h,t]/data[:stor][:η⁻] .== 0)
    @constraint(model, λ_s_1[w=1:set[:W],t=1:set[:T]], ϑ[:,w,1,t] .- φ⁺[:,w,1,t]*data[:stor][:η⁺] .+ φ⁻[:,w,1,t]/data[:stor][:η⁻] .== 0)
    @constraint(model, λ_φ⁺[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(φ̅[:,τ] for τ=1:t) .- φ⁺[:,w,h,t] .>= 0)
    @constraint(model, λ_φ⁻[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(φ̅[:,τ] for τ=1:t) .- φ⁻[:,w,h,t] .>= 0)
    @constraint(model, λ_φ[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  sum(φ̅[:,τ] for τ=1:t) .- φ⁺[:,w,h,t] .- φ⁻[:,w,h,t] .>= 0)
    @constraint(model, λ_ϑ[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(ϑ̅[:,τ] for τ=1:t) .- ϑ[:,w,h,t] .>= 0)
    # emission limits
    @constraint(model, λ_e[t=1:set[:T]],  - sum(data[:load][:w][w]*168/set[:H]*sum(sum(data[:gen_e][:e].*data[:gen_e][:h].*p[:,w,h,t])
    + sum(data[:gen_c][:e].*data[:gen_c][:h].*y[:,w,h,t]) for h in 1:set[:H]) for w in 1:set[:W])./10^6 >= - (E̅(t)*S(t)*ξ[:μ̂])[1])
    # investment limits
    @constraint(model, λ_y̅,  - sum(y̅[:,t] for t in 1:set[:T]) .>= - data[:gen_c][:y̅_max])
    @constraint(model, λ_wd[t=1:set[:T]],  - data[:gen_c][:y̅_wd]' * y̅[:,t] >= - set[:y̅_wd][t])
    @constraint(model, λ_pv[t=1:set[:T]],  - data[:gen_c][:y̅_pv]' * y̅[:,t] >= - set[:y̅_pv][t])
    # solve model
    optimize!(model)
    @info("done solving primal det. capacity expansion: $(termination_status(model))")

    sol = Dict(
    "ϑ̅" => JuMP.value.(ϑ̅),
    "φ̅" => JuMP.value.(φ̅),
    "y̅" => JuMP.value.(y̅),
    "y" => JuMP.value.(y),
    "p" => JuMP.value.(p),
    "φ⁺" => JuMP.value.(φ⁺),
    "φ⁻" => JuMP.value.(φ⁻),
    "ϑ" => JuMP.value.(ϑ),
    "obj" => JuMP.objective_value.(model)
    )

    # save solution dictionaries to json files
    save_json_determenistic = JSON.json(sol)
    open("$(outdir)determenistic.json","w") do f
        write(f, save_json_determenistic)
    end

    # save summary
    open("$(outdir)summary.txt","a") do io
       println(io,"========= in-sample solution =========")
       println(io,"det_cost = ", sol["obj"])
    end

    return sol
end
function determenistic_expansion_dual(set,data,ξ)
    # model definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # opf constraints dual variables
    @variable(model, λ_b[1:set[:T],1:set[:W],1:set[:H]])
    @variable(model, λ_f̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:E]] >= 0)
    @variable(model, λ_f̲[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:E]] >= 0)
    # generation constraints dual variables
    @variable(model, λ_p[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, λ_y[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:J]] >= 0)
    @variable(model, λ_r_e⁻[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, λ_r_e⁺[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, λ_r_c⁻[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:J]] >= 0)
    @variable(model, λ_r_c⁺[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:J]] >= 0)
    # storage constraints dual variables
    @variable(model, λ_s[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:N]])
    @variable(model, λ_s_1[1:set[:T],1:set[:W],1:data[:dims][:N]])
    @variable(model, λ_φ⁺[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:N]] >= 0)
    @variable(model, λ_φ⁻[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:N]] >= 0)
    @variable(model, λ_φ[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:N]] >= 0)
    @variable(model, λ_ϑ[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:N]] >= 0)
    # emission constraints dual variables
    @variable(model, λ_e[1:set[:T]] >= 0)
    # investment constraints dual variables
    @variable(model, λ_y̅[1:data[:dims][:J]] >= 0)
    @variable(model, λ_wd[1:set[:T]] >= 0)
    @variable(model, λ_pv[1:set[:T]] >= 0)
    # model objective
    @objective(model, Max,
    sum(sum(sum(
        ones(data[:dims][:N])'*diagm(data[:load][:k][:,t,w,h])*L(t)*S(t)*ξ[:μ̂]*λ_b[t,w,h]
        -(data[:net][:F]*data[:load][:M]*diagm(data[:load][:k][:,t,w,h])*L(t)*S(t)*ξ[:μ̂] .+ data[:net][:f̅])'*λ_f̅[t,w,h,:]
        +(data[:net][:F]*data[:load][:M]*diagm(data[:load][:k][:,t,w,h])*L(t)*S(t)*ξ[:μ̂] .- data[:net][:f̅])'*λ_f̲[t,w,h,:]
        -(diagm(data[:gen_e][:k][:,w,h])*data[:gen_e][:p̅][:,t])'*λ_p[t,w,h,:]
        for h in 1:set[:H]) for w in 1:set[:W])
        -set[:e̅][t]*λ_e[t]
        -set[:y̅_pv][t]*λ_pv[t]
        -set[:y̅_wd][t]*λ_wd[t]
        for t in 1:set[:T])
        -data[:gen_c][:y̅_max]'*λ_y̅
        -sum(
         (diagm(data[:gen_e][:r̅])*data[:gen_e][:p̅][:,t])'*λ_r_e⁻[t,w,h,:]
        +(diagm(data[:gen_e][:r̅])*data[:gen_e][:p̅][:,t])'*λ_r_e⁺[t,w,h,:]
        for h in 2:set[:H] for w in 1:set[:W] for t in 1:set[:T])
    )
    # dual constraint associated with y̅
    @constraint(model, y̅[t=1:set[:T]],
        (Q(t)*S(t)*ξ[:μ̂] .+ sum(set[:year_per_period]*(O_c(τ)*S(τ)*ξ[:μ̂]) for τ in t:set[:T]))
            .>= sum(   sum(diagm(data[:gen_c][:k][:,w,h])'*λ_y[τ,w,h,:] for h in 1:set[:H])
                    .+ sum(diagm(data[:gen_c][:r̅])'*λ_r_c⁻[τ,w,h,:] + diagm(data[:gen_c][:r̅])'*λ_r_c⁺[τ,w,h,:] for h in 2:set[:H])
                    for τ in t:set[:T] for w in 1:set[:W])
                        .- data[:gen_c][:y̅_pv] * λ_pv[t]
                            .- data[:gen_c][:y̅_wd] * λ_wd[t]
                                .- λ_y̅)
    # dual constraint associated with ϑ̅
    @constraint(model, ϑ̅[t=1:set[:T]],
        (Q_ϑ̅(t)*S(t)*ξ[:μ̂] .+ sum(set[:year_per_period]*(O_ϑ̅(τ)*S(τ)*ξ[:μ̂]) for τ in t:set[:T]))
            .>= sum(λ_ϑ[τ,w,h,:] for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))
    # dual constraint associated with φ̅
    @constraint(model, φ̅[t=1:set[:T]],
        (Q_φ̅(t)*S(t)*ξ[:μ̂] .+ sum(set[:year_per_period]*O_φ̅(τ)*S(τ)*ξ[:μ̂] for τ in t:set[:T]))
            .>= sum(λ_φ⁺[τ,w,h,:] + λ_φ⁻[τ,w,h,:] + λ_φ[τ,w,h,:] for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))
    # dual constraint associated with p
    @constraint(model, p[t=1:set[:T],w=1:set[:W],h=1:set[:H]],
        set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*V(t)*S(t)*ξ[:μ̂]
            .>= λ_b[t,w,h] * (ones(data[:dims][:N])'*data[:gen_e][:M])'
                .+ (data[:net][:F]*data[:gen_e][:M])'*(λ_f̲[t,w,h,:] .- λ_f̅[t,w,h,:])
                    .- λ_p[t,w,h,:]
                        .+ (h >= 2 ? λ_r_e⁻[t,w,h,:] .- λ_r_e⁺[t,w,h,:] : 0)
                            .+ (h < set[:H] ? λ_r_e⁺[t,w,h+1,:] .- λ_r_e⁻[t,w,h+1,:] : 0)
                                .- data[:load][:w][w]*168/set[:H] * λ_e[t] * (data[:gen_e][:e].*data[:gen_e][:h])./10^6)
    # dual constraint associated with y
    @constraint(model, y[t=1:set[:T],w=1:set[:W],h=1:set[:H]],
        set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*B(t)*S(t)*ξ[:μ̂]
            .>= λ_b[t,w,h] * (ones(data[:dims][:N])'*data[:gen_c][:M])'
                .+ (data[:net][:F]*data[:gen_c][:M])'*(λ_f̲[t,w,h,:] .- λ_f̅[t,w,h,:])
                    .- λ_y[t,w,h,:]
                        .+ (h >= 2 ? λ_r_c⁻[t,w,h,:] .- λ_r_c⁺[t,w,h,:] : 0)
                            .+ (h < set[:H] ? λ_r_c⁺[t,w,h+1,:] .- λ_r_c⁻[t,w,h+1,:] : 0)
                                .- data[:load][:w][w]*168/set[:H] * λ_e[t] * (data[:gen_c][:e].*data[:gen_c][:h])./10^6)
    # dual constraint associated with φ⁻
    @constraint(model, φ⁻[t=1:set[:T],w=1:set[:W],h=1:set[:H]],
        0 .>=
            + λ_b[t,w,h] * (ones(data[:dims][:N])'*data[:stor][:M])'
                .+ (data[:net][:F]*data[:stor][:M])'*(λ_f̲[t,w,h,:] .- λ_f̅[t,w,h,:])
                    .+ (h >= 2 ? λ_s[t,w,h,:] ./ data[:stor][:η⁻] : 0)
                        .+ (h == 1 ? λ_s_1[t,w,:] ./ data[:stor][:η⁻] : 0)
                            .- λ_φ⁻[t,w,h,:]
                                .- λ_φ[t,w,h,:])
    # dual constraint associated with φ⁺
    @constraint(model, φ⁺[t=1:set[:T],w=1:set[:W],h=1:set[:H]],
        0 .>=
            - λ_b[t,w,h] * (ones(data[:dims][:N])'*data[:stor][:M])'
                .- (data[:net][:F]*data[:stor][:M])'*(λ_f̲[t,w,h,:] .- λ_f̅[t,w,h,:])
                    .- (h >= 2 ? λ_s[t,w,h,:] * data[:stor][:η⁺] : 0)
                        .- (h == 1 ? λ_s_1[t,w,:] * data[:stor][:η⁺] : 0)
                            .- λ_φ⁺[t,w,h,:]
                                .- λ_φ[t,w,h,:])
    # dual constraint associated with ϑ
    @constraint(model, ϑ[t=1:set[:T],w=1:set[:W],h=1:set[:H]],
        0 .>= .+ (h == 1 ? λ_s_1[t,w,:] : 0)
                .+ (h >= 2 ? λ_s[t,w,h,:] : 0)
                    .- (h < set[:H] ? λ_s[t,w,h+1,:] : 0)
                        .- λ_ϑ[t,w,h,:])
    # solve model
    optimize!(model)
    @info("done solving dual det. capacity expansion: $(termination_status(model))")
    # save solution
    sol = Dict(
    "y̅" => [JuMP.dual.(y̅[t]) for t in 1:set[:T]],
    "φ̅" => [JuMP.dual.(φ̅[t]) for t in 1:set[:T]],
    "ϑ̅" => [JuMP.dual.(ϑ̅[t]) for t in 1:set[:T]],
    "obj" => JuMP.objective_value.(model),
    )
    return sol
end

## stochastic planning models
function stochastic_expansion_prim(data,set,ξ)
    # model definition
    @info("start building stochastic expansion model")
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # auxilary variables for variance penalization
    s_Y̅ =  [@variable(model, [1:data[:dims][:J]]) for t = 1:set[:T]]
    # investment decisions
    Y̅ =  [@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for t = 1:set[:T]]
    Θ̅ =  [@variable(model, [1:data[:dims][:K], 1:set[:n_t][t]]) for t = 1:set[:T]]
    Φ̅ =  [@variable(model, [1:data[:dims][:K], 1:set[:n_t][t]]) for t = 1:set[:T]]
    # operational decisions
    Y =[[[@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    Y[T][W][H][J,:]
    P =[[[@variable(model, [1:data[:dims][:G], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    P[T][W][H][J,:]
    Θ =[[[@variable(model, [1:data[:dims][:K], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    Θ[T][W][H][J,:]
    Φ⁻ =[[[@variable(model, [1:data[:dims][:K], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    Φ[T][W][H][J,:]
    Φ⁺ =[[[@variable(model, [1:data[:dims][:K], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    Φ[T][W][H][J,:]
    # auxiliary variables
    Z =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]] #    Z[T][W][H][N,:]
    @variable(model, z_f̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:E]] >= 0)
    @variable(model, x_f̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:E]] >= 0)
    @variable(model, z_p̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, x_p̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, z_r̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @variable(model, x_r̅[1:set[:T],1:set[:W],1:set[:H],1:data[:dims][:G]] >= 0)
    @info("done initializing variables")
    # objective function
    @objective(model, Min,
    sum(
    # investment costs
    + tr(S(t)*ξ[:Σ̂]*S(t)'*Q(t)'*Y̅[t])
    + tr(S(t)*ξ[:Σ̂]*S(t)'*Q_ϑ̅(t)'*Θ̅[t])
    + tr(S(t)*ξ[:Σ̂]*S(t)'*Q_φ̅(t)'*Φ̅[t])
    # fixed costs
    + set[:year_per_period]*(O_c(t)*S(t)*ξ[:μ̂])'*sum(Y̅[τ]*S(τ) for τ in 1:t)*ξ[:μ̂]
    + set[:year_per_period]*(O_ϑ̅(t)*S(t)*ξ[:μ̂])'*sum(Θ̅[τ]*S(τ) for τ in 1:t)*ξ[:μ̂]
    + set[:year_per_period]*(O_φ̅(t)*S(t)*ξ[:μ̂])'*sum(Φ̅[τ]*S(τ) for τ in 1:t)*ξ[:μ̂]
    # fuel costs
    + set[:year_per_period]*52/sum(data[:load][:w])*sum(data[:load][:w][w]*(
    168/set[:H]*sum(tr(S(t)*ξ[:Σ̂]*S(t)'*B(t)'*Y[t][w][h]) for h in 1:set[:H]) +
    168/set[:H]*sum(tr(S(t)*ξ[:Σ̂]*S(t)'*V(t)'*P[t][w][h]) for h in 1:set[:H])
    )
    for w in 1:set[:W])
    for t in 1:set[:T])
    )
    @info("done building objective function")
    # variance constraints
    @constraint(model, var_y̅[t=2:set[:T], i=1:data[:dims][:J]], [set[:α_y]*(Y̅[t]*S(t)*ξ[:μ̂])[i];ξ[:Σ̅]*(Y̅[t]*S(t))[i,:]] in SecondOrderCone())
    @constraint(model, var_ϑ̅[t=2:set[:T], i=1:data[:dims][:K]], [set[:α_ϑ]*(Θ̅[t]*S(t)*ξ[:μ̂])[i];ξ[:Σ̅]*(Θ̅[t]*S(t))[i,:]] in SecondOrderCone())
    @constraint(model, var_φ̅[t=2:set[:T], i=1:data[:dims][:K]], [set[:α_φ]*(Φ̅[t]*S(t)*ξ[:μ̂])[i];ξ[:Σ̅]*(Φ̅[t]*S(t))[i,:]] in SecondOrderCone())
    # power conservation equations
    @constraint(model, aux_con_1[t=1:set[:T],w=1:set[:W],h=1:set[:H]], Z[t][w][h] .== data[:gen_e][:M]*P[t][w][h]
        .+ data[:gen_c][:M]*Y[t][w][h]
        .+ data[:stor][:M] * Φ⁻[t][w][h]
        .- data[:load][:M]*diagm(data[:load][:k][:,t,w,h])*L(t)
        .- data[:stor][:M] * Φ⁺[t][w][h])
    @constraint(model, bal[t=1:set[:T],w=1:set[:W],h=1:set[:H]], ones(data[:dims][:N])' * Z[t][w][h] .== 0)
    @info("done building power conservation equations")
    # power flow limits
    if set[:cc_ref] != "DRO-DS"
        @constraint(model,   μ̅[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            [(data[:net][:f̅] .- data[:net][:F]*Z[t][w][h]*S(t)*ξ[:μ̂])[i];
                Ψ(set[:ε̅_f]/(2*data[:dims][:E]))*ξ[:Σ̅]*(.-data[:net][:F]*Z[t][w][h]*S(t))[i,:]] in SecondOrderCone())
        @constraint(model,   μ̲[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            [(data[:net][:F]*Z[t][w][h]*S(t)*ξ[:μ̂] .+ data[:net][:f̅])[i];
                Ψ(set[:ε̅_f]/(2*data[:dims][:E]))*ξ[:Σ̅]*(data[:net][:F]*Z[t][w][h]*S(t))[i,:]] in SecondOrderCone())
    else
        @constraint(model,   con_f̅_1[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            [sqrt(set[:ε̅_f]/(data[:dims][:E]))*(data[:net][:f̅][i] - x_f̅[t,w,h,i]);
                vcat(ξ[:Σ̅]*(data[:net][:F]*Z[t][w][h]*S(t))[i,:],z_f̅[t,w,h,i])]
                    in SecondOrderCone())
        @constraint(model,   con_f̅_2[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            (data[:net][:F]*Z[t][w][h]*S(t)*ξ[:μ̂])[i] <= z_f̅[t,w,h,i] + x_f̅[t,w,h,i])
        @constraint(model,   con_f̅_3[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            (data[:net][:F]*Z[t][w][h]*S(t)*ξ[:μ̂])[i] >= - z_f̅[t,w,h,i] - x_f̅[t,w,h,i])
        @constraint(model,   con_f̅_4[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
            data[:net][:f̅][i] >= x_f̅[t,w,h,i])
    end
    @info("done enforcing power flow limits")
    # existing power generation limits
    for i in 1:data[:dims][:G]
        if set[:cc_ref] != "DRO-DS"
            for t in 1:set[:T], w in 1:set[:W], h in 1:set[:H]
                @constraint(model, [data[:gen_e][:k][i,w,h].*data[:gen_e][:p̅][i,t] .- P[t][w][h][i,:]'*S(t)*ξ[:μ̂] ; Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(P[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
                @constraint(model, [(P[t][w][h]*S(t)*ξ[:μ̂])[i] ; Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(P[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
            end
            for t in 1:set[:T], w in 1:set[:W], h in 2:set[:H]
                @constraint(model, [P[t][w][h][i,:]'*S(t)*ξ[:μ̂] - P[t][w][h-1][i,:]'*S(t)*ξ[:μ̂] + data[:gen_e][:r̅][i] * data[:gen_e][:p̅][i]; Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J])) * ξ[:Σ̅] * (P[t][w][h][i,:]'*S(t) - P[t][w][h-1][i,:]'*S(t))'] in SecondOrderCone())
                @constraint(model, [P[t][w][h-1][i,:]'*S(t)*ξ[:μ̂] - P[t][w][h][i,:]'*S(t)*ξ[:μ̂] + data[:gen_e][:r̅][i] * data[:gen_e][:p̅][i]; Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J])) * ξ[:Σ̅] * (P[t][w][h-1][i,:]'*S(t) - P[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
            end
        else
            for t in 1:set[:T], w in 1:set[:W], h in 1:set[:H]
                @constraint(model, [sqrt(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J])*2)*(0.5*data[:gen_e][:k][i,w,h].*data[:gen_e][:p̅][i,t] - x_p̅[t,w,h,i]);
                                        vcat(ξ[:Σ̅]*(P[t][w][h][i,:]'*S(t))',z_p̅[t,w,h,i])] in SecondOrderCone())
                @constraint(model, P[t][w][h][i,:]'*S(t)*ξ[:μ̂] - 0.5*data[:gen_e][:k][i,w,h].*data[:gen_e][:p̅][i,t] <=  z_p̅[t,w,h,i] + x_p̅[t,w,h,i])
                @constraint(model, P[t][w][h][i,:]'*S(t)*ξ[:μ̂] - 0.5*data[:gen_e][:k][i,w,h].*data[:gen_e][:p̅][i,t] >= -z_p̅[t,w,h,i] - x_p̅[t,w,h,i])
                @constraint(model, 0.5*data[:gen_e][:k][i,w,h].*data[:gen_e][:p̅][i,t] >= x_p̅[t,w,h,i])
            end
            for t in 1:set[:T], w in 1:set[:W], h in 2:set[:H]
                @constraint(model, [sqrt(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J])*2)*(0.5*(data[:gen_e][:r̅][i] + data[:gen_e][:r̅][i]).*data[:gen_e][:p̅][i,t] - x_r̅[t,w,h,i]);
                                        vcat(ξ[:Σ̅]*((P[t][w][h][i,:]-P[t][w][h-1][i,:])'*S(t))',z_r̅[t,w,h,i])] in SecondOrderCone())
                @constraint(model, vec((P[t][w][h][i,:]-P[t][w][h-1][i,:])'*S(t))'*ξ[:μ̂] - 0.5*(data[:gen_e][:r̅][i] + data[:gen_e][:r̅][i]).*data[:gen_e][:p̅][i,t] <= z_r̅[t,w,h,i] + x_r̅[t,w,h,i])
                @constraint(model, vec((P[t][w][h][i,:]-P[t][w][h-1][i,:])'*S(t))'*ξ[:μ̂] - 0.5*(data[:gen_e][:r̅][i] + data[:gen_e][:r̅][i]).*data[:gen_e][:p̅][i,t] >= - z_r̅[t,w,h,i] - x_r̅[t,w,h,i])
                @constraint(model, 0.5*(data[:gen_e][:r̅][i] + data[:gen_e][:r̅][i]).*data[:gen_e][:p̅][i,t] >= x_r̅[t,w,h,i])
            end
        end
        @info("done with enforcing $(i)-th existing power generator limits")
    end
    # candidate power generation limits
    for i in 1:data[:dims][:J]
        for t in 1:set[:T], w in 1:set[:W], h in 1:set[:H]
            @constraint(model, [(data[:gen_c][:k][i,w,h] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) - Y[t][w][h][i,:]'*S(t))*ξ[:μ̂];
                Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(data[:gen_c][:k][i,w,h] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) - Y[t][w][h][i,:]'*S(t))']
                    in SecondOrderCone())
            @constraint(model, [Y[t][w][h][i,:]'*S(t)*ξ[:μ̂];Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(Y[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
        end
        for t in 1:set[:T], w in 1:set[:W], h in 2:set[:H]
            @constraint(model,[(data[:gen_c][:r̅][i] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) + Y[t][w][h][i,:]'*S(t) - Y[t][w][h-1][i,:]'*S(t))*ξ[:μ̂];
                Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(data[:gen_c][:r̅][i] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) + Y[t][w][h][i,:]'*S(t) - Y[t][w][h-1][i,:]'*S(t))']
                    in SecondOrderCone())
            @constraint(model,[(data[:gen_c][:r̅][i] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) + Y[t][w][h-1][i,:]'*S(t) - Y[t][w][h][i,:]'*S(t))*ξ[:μ̂];
                Ψ(set[:ε̅_g]/(2*data[:dims][:G] + 2*data[:dims][:J]))*ξ[:Σ̅]*(data[:gen_c][:r̅][i] * sum(Y̅[τ][i,:]'*S(τ) for τ in 1:t) + Y[t][w][h-1][i,:]'*S(t) - Y[t][w][h][i,:]'*S(t))']
                    in SecondOrderCone())
        end
        @info("done with enforcing $(i)-th candidate power generator limits")
    end
    # storage equations
    @constraint(model, SoC_t[t=1:set[:T],w=1:set[:W],h=2:set[:H]], Θ[t][w][h] .- Θ[t][w][h-1] .- Φ⁺[t][w][h]*data[:stor][:η⁺] .+ Φ⁻[t][w][h]/data[:stor][:η⁻]  .== 0)
    @constraint(model, SoC_1[t=1:set[:T],w=1:set[:W]            ], Θ[t][w][1] .- Φ⁺[t][w][1]*data[:stor][:η⁺] .+ Φ⁻[t][w][1]/data[:stor][:η⁻] .== 0)
    @info("done building storage")
    # storage limits
    @constraint(model, lim_φ̅_pl[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [(sum(Φ̅[τ][i,:]'*S(τ) for τ=1:t) - Φ⁺[t][w][h][i,:]'*S(t))*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(sum(Φ̅[τ][i,:]'*S(τ) for τ=1:t) - Φ⁺[t][w][h][i,:]'*S(t))']
            in SecondOrderCone())
    @constraint(model, lim_φ̅_mn[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [(sum(Φ̅[τ][i,:]'*S(τ) for τ=1:t) - Φ⁻[t][w][h][i,:]'*S(t))*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(sum(Φ̅[τ][i,:]'*S(τ) for τ=1:t) - Φ⁻[t][w][h][i,:]'*S(t))']
            in SecondOrderCone())
    @constraint(model, lim_φ̲_pl[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [Φ⁺[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(Φ⁺[t][w][h][i,:]'*S(t))']
            in SecondOrderCone())
    @constraint(model, lim_φ̲_mn[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [Φ⁻[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(Φ⁻[t][w][h][i,:]'*S(t))']
            in SecondOrderCone())
    @constraint(model, lim_ϑ̅[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [(sum(Θ̅[τ][i,:]'*S(τ) for τ=1:t) - Θ[t][w][h][i,:]'*S(t))*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(sum(Θ̅[τ][i,:]'*S(τ) for τ=1:t) - Θ[t][w][h][i,:]'*S(t))']
            in SecondOrderCone())
    @constraint(model, lim_ϑ̲[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:K]],
        [Θ[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_s]/(6*data[:dims][:K]))*ξ[:Σ̅]*(Θ[t][w][h][i,:]'*S(t))']
             in SecondOrderCone())
    @info("done with enforcing storage limits")
    # capacity expansion limits
    @constraint(model, y̲̅[t=1:set[:T], i=1:data[:dims][:J]],
        [Y̅[t][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*(Y̅[t][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, φ̲̅[t=1:set[:T], i=1:data[:dims][:K]],
        [Φ̅[t][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*(Φ̅[t][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, ϑ̲̅[t=1:set[:T], i=1:data[:dims][:K]],
        [Θ̅[t][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*(Θ̅[t][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, y̅̅[i=1:data[:dims][:J];data[:gen_c][:y̅][i] != -1],
        [data[:gen_c][:y̅][i] .- sum(Y̅[t][i,:]'*S(t) for t in 1:set[:T])*ξ[:μ̂]; Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*sum(Y̅[t][i,:]'*S(t) for t in 1:set[:T])'] in SecondOrderCone())
    @constraint(model, τ_w[t=1:set[:T]],
        [set[:y̅_wd][t] - sum(Y̅[t][j,:] for j in 1:data[:dims][:J] if occursin("Wind", data[:gen_c][:tech][j]) == true)'*S(t)*ξ[:μ̂];
            Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*(sum(Y̅[t][j,:] for j in 1:data[:dims][:J] if occursin("Wind", data[:gen_c][:tech][j]) == true)'*S(t))'] in SecondOrderCone())
    @constraint(model, τ_p[t=1:set[:T]],
        [set[:y̅_pv][t] - sum(Y̅[t][j,:] for j in 1:data[:dims][:J] if occursin("PV",   data[:gen_c][:tech][j]) == true)'*S(t)*ξ[:μ̂];
            Ψ(set[:ε̅_i]/(2*data[:dims][:K] + 2*data[:dims][:K] + 2))*ξ[:Σ̅]*(sum(Y̅[t][j,:] for j in 1:data[:dims][:J] if occursin("PV",   data[:gen_c][:tech][j]) == true)'*S(t))'] in SecondOrderCone())
    @info("done with enforcing capacity expansion limits")
    # emission limits
    @constraint(model, e_max[t=1:set[:T]],
        [(E̅(t)*S(t)*ξ[:μ̂])[1] .- sum(data[:load][:w][w]*168/set[:H]*sum(sum(data[:gen_e][:e][i]*data[:gen_e][:h][i]*P[t][w][h][i,:]'*S(t) for i in 1:data[:dims][:G]) + sum(data[:gen_c][:e][i]*data[:gen_c][:h][i]*Y[t][w][h][i,:]'*S(t) for i in 1:data[:dims][:J]) for h in 1:set[:H]) for w in 1:set[:W])*ξ[:μ̂]./10^6;
            Ψ(set[:ε̅_e])*ξ[:Σ̅]*((E̅(t)*S(t))[:] .- sum(data[:load][:w][w]*168/set[:H]*sum(sum(data[:gen_e][:e][i]*data[:gen_e][:h][i]*P[t][w][h][i,:]'*S(t) for i in 1:data[:dims][:G]) + sum(data[:gen_c][:e][i]*data[:gen_c][:h][i]*Y[t][w][h][i,:]'*S(t) for i in 1:data[:dims][:J]) for h in 1:set[:H]) for w in 1:set[:W])'./10^6)] in SecondOrderCone());
    @info("done with enforcing emission limits")
    # solve model
    optimize!(model)
    @info("done solving stochastic capacity expansion: $(termination_status(model))")
    # prepare "average" results
    ϑ̅ = zeros(data[:dims][:K],set[:T])
    y̅ = zeros(data[:dims][:J],set[:T])
    y = zeros(data[:dims][:J],set[:W],set[:H],set[:T])
    p = zeros(data[:dims][:J],set[:W],set[:H],set[:T])
    for j in 1:data[:dims][:J], t in 1:set[:T]
        y̅[j,t] = JuMP.value.(Y̅[t])[j,:]'*S(t)*ξ[:μ̂]
        for w in 1:set[:W], h in 1:set[:H]
            y[j,w,h,t] = JuMP.value.(Y[t][w][h])[j,:]'*S(t)*ξ[:μ̂]
        end
    end
    for j in 1:data[:dims][:G], t in 1:set[:T], w in 1:set[:W], h in 1:set[:H]
        p[j,w,h,t] = JuMP.value.(P[t][w][h])[j,:]'*S(t)*ξ[:μ̂]
    end
    for j in 1:data[:dims][:K], t in 1:set[:T]
        ϑ̅[j,t] = JuMP.value.(Θ̅[t])[j,:]'*S(t)*ξ[:μ̂]
    end
    cost = sum(
            # investment costs
            + tr(S(t)*ξ[:Σ̂]*S(t)'*Q(t)'*JuMP.value.(Y̅[t]))
            + tr(S(t)*ξ[:Σ̂]*S(t)'*Q_ϑ̅(t)'*JuMP.value.(Θ̅[t]))
            + tr(S(t)*ξ[:Σ̂]*S(t)'*Q_φ̅(t)'*JuMP.value.(Φ̅[t]))
            # fixed costs
            + set[:year_per_period]*(O_c(t)*S(t)*ξ[:μ̂])'*sum(JuMP.value.(Y̅[τ])*S(τ) for τ in 1:t)*ξ[:μ̂]
            + set[:year_per_period]*(O_ϑ̅(t)*S(t)*ξ[:μ̂])'*sum(JuMP.value.(Θ̅[τ])*S(τ) for τ in 1:t)*ξ[:μ̂]
            + set[:year_per_period]*(O_φ̅(t)*S(t)*ξ[:μ̂])'*sum(JuMP.value.(Φ̅[τ])*S(τ) for τ in 1:t)*ξ[:μ̂]
            # fuel costs
            + set[:year_per_period]*52/sum(data[:load][:w])*sum(data[:load][:w][w]*(
            168/set[:H]*sum(tr(S(t)*ξ[:Σ̂]*S(t)'*B(t)'*JuMP.value.(Y[t][w][h])) for h in 1:set[:H]) +
            168/set[:H]*sum(tr(S(t)*ξ[:Σ̂]*S(t)'*V(t)'*JuMP.value.(P[t][w][h])) for h in 1:set[:H])
            )
            for w in 1:set[:W])
            for t in 1:set[:T])
    Y̅_std = zeros(data[:dims][:J],set[:T])
    Θ̅_std = zeros(data[:dims][:K],set[:T])
    Φ̅_std = zeros(data[:dims][:K],set[:T])
    for t in 1:set[:T], i in 1:data[:dims][:J]
        Y̅_std[i,t] = norm(ξ[:Σ̅]*(JuMP.value.(Y̅[t])*S(t))[i,:],2)
    end
    for t in 1:set[:T], i in 1:data[:dims][:K]
        Θ̅_std[i,t] = norm(ξ[:Σ̅]*(JuMP.value.(Θ̅[t])*S(t))[i,:],2)
        Φ̅_std[i,t] = norm(ξ[:Σ̅]*(JuMP.value.(Φ̅[t])*S(t))[i,:],2)
    end

    sol_ = Dict(
    "obj" => JuMP.objective_value.(model),
    "cost" => cost,
    "Y̅" => [JuMP.value.(Y̅[t]) for t in 1:set[:T]],
    "Y" => [[[JuMP.value.(Y[t][w][h]) for h in 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]],
    "P" => [[[JuMP.value.(P[t][w][h]) for h in 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]],
    # "s_Y̅" => [JuMP.value.(s_Y̅[t]) for t in 1:set[:T]],
    "Y̅_std" => Y̅_std, "Θ̅_std" => Θ̅_std, "Φ̅_std" => Φ̅_std,
    "Θ̅" => [JuMP.value.(Θ̅[t]) for t in 1:set[:T]],
    "Φ̅" => [JuMP.value.(Φ̅[t]) for t in 1:set[:T]],
    "y̅" => y̅,
    "y" => y,
    "p" => p,
    "ϑ̅" => ϑ̅,
    )


    save_json_stochastic = JSON.json(sol_)
    open("$(outdir)stochastic.json","w") do f
        write(f, save_json_stochastic)
    end
    # save summary
    open("$(outdir)summary.txt","a") do io
       println(io,"sto_cost = ", sol_["cost"])
       println(io,"tot_std_Y̅ = ", sum(Y̅_std))
       println(io,"tot_std_Θ̅ = ", sum(Θ̅_std))
       println(io,"tot_std_Φ̅ = ", sum(Φ̅_std))
       println(io,"========= out-of-sample solution =========")
    end
    return sol_
end
function stochastic_expansion_dual(data,set,ξ)
    # model definition
    @info("start building stochastic dual expansion model")
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # opf constraints dual variables
    Λ_b =[[[@variable(model, [1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_f̅ =[[[@variable(model, [1:data[:dims][:E], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_f̲ =[[[@variable(model, [1:data[:dims][:E], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    # generation constraints dual variables
    Λ_p =[[[@variable(model, [1:data[:dims][:G], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_y =[[[@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_r_e⁻ =[[[@variable(model, [1:data[:dims][:G], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_r_e⁺ =[[[@variable(model, [1:data[:dims][:G], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_r_c⁻ =[[[@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_r_c⁺ =[[[@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    # storage constraints dual variables
    Λ_s   =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_s_1 =[[ @variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_φ⁺  =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_φ⁻  =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_φ   =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    Λ_ϑ   =[[[@variable(model, [1:data[:dims][:N], 1:set[:n_t][t]]) for h = 1:set[:H]] for w in 1:set[:W]] for t in 1:set[:T]]
    # emission constraints dual variables
    Λ_e   =[@variable(model, [1:set[:n_t][t]]) for t in 1:set[:T]]
    # investment constraints dual variables
    Λ_y̅   =[@variable(model, [1:data[:dims][:J], 1:set[:n_t][t]]) for t in set[:T]]
    Λ_wd  =[@variable(model, [1:set[:n_t][t]]) for t in 1:set[:T]]
    Λ_pv  =[@variable(model, [1:set[:n_t][t]]) for t in 1:set[:T]]

    num_variables(model)
    @info("done dual variable definitions")
    # model objective
    @objective(model, Max,
    sum(sum(sum(
        tr(S(t)*ξ[:Σ̂]*S(t)'*(ones(data[:dims][:N])'*diagm(data[:load][:k][:,t,w,h])*L(t))'*Λ_b[t][w][h]')
        -tr(S(t)*ξ[:Σ̂]*S(t)'*(data[:net][:F]*data[:load][:M]*diagm(data[:load][:k][:,t,w,h])*L(t))'*Λ_f̅[t][w][h])
        -(data[:net][:f̅])'*Λ_f̅[t][w][h]*S(t)*ξ[:μ̂]
        +tr(S(t)*ξ[:Σ̂]*S(t)'*(data[:net][:F]*data[:load][:M]*diagm(data[:load][:k][:,t,w,h])*L(t))'*Λ_f̲[t][w][h])
        -(data[:net][:f̅])'*Λ_f̲[t][w][h]*S(t)*ξ[:μ̂]

        -(diagm(data[:gen_e][:k][:,w,h])*data[:gen_e][:p̅][:,t])'*Λ_p[t][w][h]*S(t)*ξ[:μ̂]
        for h in 1:set[:H]) for w in 1:set[:W])
        -set[:e̅][t]*Λ_e[t]'*S(t)*ξ[:μ̂]
        -set[:y̅_pv][t]*Λ_pv[t]'*S(t)*ξ[:μ̂]
        -set[:y̅_wd][t]*Λ_wd[t]'*S(t)*ξ[:μ̂]
        for t in 1:set[:T])
        -data[:gen_c][:y̅_max]'*Λ_y̅[1]*ξ[:μ̂]
        -sum(
         (diagm(data[:gen_e][:r̅])*data[:gen_e][:p̅][:,t])'*Λ_r_e⁻[t][w][h]*S(t)*ξ[:μ̂]
        +(diagm(data[:gen_e][:r̅])*data[:gen_e][:p̅][:,t])'*Λ_r_e⁺[t][w][h]*S(t)*ξ[:μ̂]
        for h in 2:set[:H] for w in 1:set[:W] for t in 1:set[:T])
    )
    @info("done building objective function")
    # dual constraint associated with y̅
    @constraint(model, y̅[t=1:set[:T],i=1:data[:dims][:J]],
        [((Q(t)*S(t) .+ sum(set[:year_per_period]*(O_c(τ)*S(τ)) for τ in t:set[:T]))
        .- sum(   sum(diagm(data[:gen_c][:k][:,w,h])'*Λ_y[τ][w][h]*S(τ) for h in 1:set[:H])
                .+ sum(diagm(data[:gen_c][:r̅])'*Λ_r_c⁻[τ][w][h]*S(τ) + diagm(data[:gen_c][:r̅])'*Λ_r_c⁺[τ][w][h]*S(τ) for h in 2:set[:H])
                for τ in t:set[:T] for w in 1:set[:W])
                    .+ data[:gen_c][:y̅_pv] * Λ_pv[t]'*S(t)
                        .+ data[:gen_c][:y̅_pv] * Λ_pv[t]'*S(t)
                            .+ data[:gen_c][:y̅_wd] * Λ_wd[t]'*S(t)
                                .+ Λ_y̅[1])[i,:]'*ξ[:μ̂]
        ;Ψ(set[:ε̅])*ξ[:Σ̅]*
        ((Q(t)*S(t) .+ sum(set[:year_per_period]*(O_c(τ)*S(τ)) for τ in t:set[:T]))
        .- sum(   sum(diagm(data[:gen_c][:k][:,w,h])'*Λ_y[τ][w][h]*S(τ) for h in 1:set[:H])
                .+ sum(diagm(data[:gen_c][:r̅])'*Λ_r_c⁻[τ][w][h]*S(τ) + diagm(data[:gen_c][:r̅])'*Λ_r_c⁺[τ][w][h]*S(τ) for h in 2:set[:H])
                for τ in t:set[:T] for w in 1:set[:W])
                    .+ data[:gen_c][:y̅_pv] * Λ_pv[t]'*S(t)
                        .+ data[:gen_c][:y̅_pv] * Λ_pv[t]'*S(t)
                            .+ data[:gen_c][:y̅_wd] * Λ_wd[t]'*S(t)
                                .+ Λ_y̅[1])[i,:]
        ] in SecondOrderCone())
    @info("done building y̅ constraints")
    # dual constraint associated with ϑ̅
    @constraint(model, ϑ̅[t=1:set[:T],i=1:data[:dims][:N]],
    [
    (Q_ϑ̅(t)*S(t)
        .+ sum(set[:year_per_period]*(O_ϑ̅(τ)*S(τ)) for τ in t:set[:T])
            .- sum(Λ_ϑ[τ][w][h]*S(τ) for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))[i,:]'*ξ[:μ̂]
    ;
    Ψ(set[:ε̅])*ξ[:Σ̅]*(Q_ϑ̅(t)*S(t)
        .+ sum(set[:year_per_period]*(O_ϑ̅(τ)*S(τ)) for τ in t:set[:T])
            .- sum(Λ_ϑ[τ][w][h]*S(τ) for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))[i,:]
    ] in SecondOrderCone())
    @info("done building ϑ̅ constraints")
    # # dual constraint associated with φ̅
    @constraint(model, φ̅[t=1:set[:T],i=1:data[:dims][:N]],
    [
    (Q_φ̅(t)*S(t) .+ sum(set[:year_per_period]*O_φ̅(τ)*S(τ) for τ in t:set[:T])
        .- sum(Λ_φ⁺[τ][w][h]*S(τ) + Λ_φ⁻[τ][w][h]*S(τ) + Λ_φ[τ][w][h]*S(τ) for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))[i,:]'*ξ[:μ̂]
    ;
    Ψ(set[:ε̅])*ξ[:Σ̅]*(Q_φ̅(t)*S(t) .+ sum(set[:year_per_period]*O_φ̅(τ)*S(τ) for τ in t:set[:T])
        .- sum(Λ_φ⁺[τ][w][h]*S(τ) + Λ_φ⁻[τ][w][h]*S(τ) + Λ_φ[τ][w][h]*S(τ) for τ in t:set[:T] for w in 1:set[:W] for h in 1:set[:H]))[i,:]
    ] in SecondOrderCone())
    @info("done building φ̅ constraints")
    # dual constraint associated with p
    @constraint(model, p[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:G]],
        [(set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*V(t)
            .- (ones(data[:dims][:N])'*data[:gen_e][:M])'*Λ_b[t][w][h]'
                .- (data[:net][:F]*data[:gen_e][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                    .+ Λ_p[t][w][h]
                        .- (h >= 2 ? Λ_r_e⁻[t][w][h] .- Λ_r_e⁺[t][w][h] : 0)
                            .- (h < set[:H] ? Λ_r_e⁺[t][w][h+1] .- Λ_r_e⁻[t][w][h+1] : 0)
                                .+ data[:load][:w][w]*168/set[:H] * (data[:gen_e][:e].*data[:gen_e][:h])./10^6 * Λ_e[t]')[i,:]'*S(t)*ξ[:μ̂]
        ;Ψ(set[:ε̅])*ξ[:Σ̅]*
        ((set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*V(t)
            .- (ones(data[:dims][:N])'*data[:gen_e][:M])'*Λ_b[t][w][h]'
                .- (data[:net][:F]*data[:gen_e][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                    .+ Λ_p[t][w][h]
                        .- (h >= 2 ? Λ_r_e⁻[t][w][h] .- Λ_r_e⁺[t][w][h] : 0)
                            .- (h < set[:H] ? Λ_r_e⁺[t][w][h+1] .- Λ_r_e⁻[t][w][h+1] : 0)
                                .+ data[:load][:w][w]*168/set[:H] * (data[:gen_e][:e].*data[:gen_e][:h])./10^6 * Λ_e[t]')[i,:]'*S(t))'
        ] in SecondOrderCone())
    @info("done building p constraints")
    # dual constraint associated with y
    @constraint(model, y[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:J]],
        [(set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*B(t)
            .- (ones(data[:dims][:N])'*data[:gen_c][:M])'*Λ_b[t][w][h]'
                .- (data[:net][:F]*data[:gen_c][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                    .+ Λ_y[t][w][h]
                        .- (h >= 2 ? Λ_r_c⁻[t][w][h] .- Λ_r_c⁺[t][w][h] : 0)
                            .- (h < set[:H] ? Λ_r_c⁺[t][w][h+1] .- Λ_r_c⁻[t][w][h+1] : 0)
                                .+ data[:load][:w][w]*168/set[:H] * (data[:gen_c][:e].*data[:gen_c][:h])./10^6 * Λ_e[t]')[i,:]'*S(t)*ξ[:μ̂]
        ;Ψ(set[:ε̅])*ξ[:Σ̅]*
        ((set[:year_per_period]*52/sum(data[:load][:w])*data[:load][:w][w]*168/set[:H]*B(t)
            .- (ones(data[:dims][:N])'*data[:gen_c][:M])'*Λ_b[t][w][h]'
                .- (data[:net][:F]*data[:gen_c][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                    .+ Λ_y[t][w][h]
                        .- (h >= 2 ? Λ_r_c⁻[t][w][h] .- Λ_r_c⁺[t][w][h] : 0)
                            .- (h < set[:H] ? Λ_r_c⁺[t][w][h+1] .- Λ_r_c⁻[t][w][h+1] : 0)
                                .+ data[:load][:w][w]*168/set[:H] * (data[:gen_c][:e].*data[:gen_c][:h])./10^6 * Λ_e[t]')[i,:]'*S(t))'
        ] in SecondOrderCone())
    @info("done building y constraints")
    # dual constraint associated with φ⁻
    @constraint(model, φ⁻[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
    [
    ((Λ_φ[t][w][h] .+ Λ_φ⁻[t][w][h] .- (h == 1 ? Λ_s_1[t][w] ./ data[:stor][:η⁻] : 0)
        .- (h >= 2 ? Λ_s[t][w][h] ./ data[:stor][:η⁻] : 0)
            .- (data[:net][:F]*data[:stor][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                .- (ones(data[:dims][:N])'*data[:stor][:M])' * Λ_b[t][w][h]')*S(t))[i,:]'*ξ[:μ̂]
    ;
    Ψ(set[:ε̅])*ξ[:Σ̅]*((Λ_φ[t][w][h] .+ Λ_φ⁻[t][w][h] .- (h == 1 ? Λ_s_1[t][w] ./ data[:stor][:η⁻] : 0)
        .- (h >= 2 ? Λ_s[t][w][h] ./ data[:stor][:η⁻] : 0)
            .- (data[:net][:F]*data[:stor][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
                .- (ones(data[:dims][:N])'*data[:stor][:M])' * Λ_b[t][w][h]')*S(t))[i,:]
    ] in SecondOrderCone())
    @info("done building φ⁻ constraints")
    # dual constraint associated with φ⁺
    @constraint(model, φ⁺[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
    [(((ones(data[:dims][:N])'*data[:stor][:M])' * Λ_b[t][w][h]'
        .+ (data[:net][:F]*data[:stor][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
            .+ (h >= 2 ? Λ_s[t][w][h] * data[:stor][:η⁺] : 0)
                .+ (h == 1 ? Λ_s_1[t][w] * data[:stor][:η⁺] : 0)
                    .+ Λ_φ⁺[t][w][h]
                        .+ Λ_φ[t][w][h])*S(t))[i,:]'*ξ[:μ̂]
    ;
    Ψ(set[:ε̅])*ξ[:Σ̅]*(((ones(data[:dims][:N])'*data[:stor][:M])' * Λ_b[t][w][h]'
        .+ (data[:net][:F]*data[:stor][:M])'*(Λ_f̲[t][w][h] .- Λ_f̅[t][w][h])
            .+ (h >= 2 ? Λ_s[t][w][h] * data[:stor][:η⁺] : 0)
                .+ (h == 1 ? Λ_s_1[t][w] * data[:stor][:η⁺] : 0)
                    .+ Λ_φ⁺[t][w][h]
                        .+ Λ_φ[t][w][h])*S(t))[i,:]
    ] in SecondOrderCone())
    @info("done building φ⁺ constraints")
    # dual constraint associated with ϑ
    @constraint(model, ϑ[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
    [
    ((Λ_ϑ[t][w][h] .+ (h < set[:H] ? Λ_s[t][w][h+1] : 0)
            .- (h >= 2 ? Λ_s[t][w][h] : 0) .- (h == 1 ? Λ_s_1[t][w] : 0))*S(t))[i,:]'*ξ[:μ̂]
    ;
    Ψ(set[:ε̅])*ξ[:Σ̅]*((Λ_ϑ[t][w][h] .+ (h < set[:H] ? Λ_s[t][w][h+1] : 0)
            .- (h >= 2 ? Λ_s[t][w][h] : 0) .- (h == 1 ? Λ_s_1[t][w] : 0))*S(t))[i,:]
    ] in SecondOrderCone())
    @info("done building ϑ constraints")
    # non-negativity cosntraints
    @constraint(model, λ_f̅[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
        [Λ_f̅[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_f̅[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_f̲[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:E]],
        [Λ_f̲[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_f̲[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_p[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:G]],
        [Λ_p[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_p[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_y[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:J]],
        [Λ_y[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_y[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_r_e⁻[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:G]],
        [Λ_r_e⁻[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_r_e⁻[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_r_e⁺[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:G]],
        [Λ_r_e⁺[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_r_e⁺[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_r_c⁻[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:J]],
        [Λ_r_c⁻[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_r_c⁻[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_r_c⁺[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:J]],
        [Λ_r_c⁺[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_r_c⁺[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_φ[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
        [Λ_φ[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_φ[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_φ⁺[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
        [Λ_φ⁺[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_φ⁺[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_φ⁻[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
        [Λ_φ⁻[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_φ⁻[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_ϑ[t=1:set[:T],w=1:set[:W],h=1:set[:H],i=1:data[:dims][:N]],
        [Λ_ϑ[t][w][h][i,:]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_ϑ[t][w][h][i,:]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_e[t=1:set[:T]],
        [Λ_e[t]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_e[t]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_y̅[j=1:data[:dims][:J]],
        [Λ_y̅[1][j,:]'*S(set[:T])*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_y̅[1][j,:]'*S(set[:T]))'] in SecondOrderCone())
    @constraint(model, λ_wd[t=1:set[:T]],
        [Λ_wd[t]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_wd[t]'*S(t))'] in SecondOrderCone())
    @constraint(model, λ_pv[t=1:set[:T]],
        [Λ_pv[t]'*S(t)*ξ[:μ̂]; Ψ(set[:ε̅])*ξ[:Σ̅]*(Λ_pv[t]'*S(t))'] in SecondOrderCone())
    @info("done building non-negativity constraints")
    # solve model
    optimize!(model)
    @info("done solving dual sto. capacity expansion: $(termination_status(model))")
    return sol_ = Dict("obj" => JuMP.objective_value.(model))
end
function SAA_gen_exp(data,set,ξ,ST)
    Ns = (data[:BF]^(data[:T]-1))
    # model definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    # model variables
    @variable(model, y̅[1:data[:J],1:data[:T],s=1:Ns] >= 0)
    @variable(model, y[1:data[:J],1:data[:H],1:data[:T],s=1:Ns] >= 0)
    @variable(model, p[1:data[:G],1:data[:H],1:data[:T],s=1:Ns] >= 0)
    # model objective
    @objective(model, Min,
                    1/Ns * sum(ST[:Q_g][:,s,t]'*y̅[:,t,s]
                    + sum(ST[:B][:,s,t]'*y[:,h,t,s] for h in 1:data[:H])
                    + sum(ST[:V][:,s,t]'*p[:,h,t,s] for h in 1:data[:H])
                    for t in 1:data[:T] for s in 1:Ns)
    )
    # power flow eqations
    @constraint(model, con_μ[h=1:data[:H],t=1:data[:T],s=1:Ns], ones(data[:N])'*(data[:M_e] * p[:,h,t,s] + data[:M_c] * y[:,h,t,s] - data[:M_l] * (data[:k_l][:,h] .* ST[:L][:,s,t])) .== 0)
    @constraint(model, con_ζ̅[h=1:data[:H],t=1:data[:T],s=1:Ns],  data[:F]*(data[:M_e] * p[:,h,t,s] + data[:M_c] * y[:,h,t,s] - data[:M_l] * (data[:k_l][:,h] .* ST[:L][:,s,t])) .<= data[:f̅])
    @constraint(model, con_ζ̲[h=1:data[:H],t=1:data[:T],s=1:Ns], - data[:f̅] .<= data[:F]*(data[:M_e] * p[:,h,t,s] + data[:M_c] * y[:,h,t,s] - data[:M_l] * (data[:k_l][:,h] .* ST[:L][:,s,t])))
    # generation limits
    @constraint(model, con_κ_e[h=1:data[:H],t=1:data[:T],s=1:Ns], p[:,h,t,s] .<= data[:k_e][:,h] .* data[:p̅])
    @constraint(model, con_κ_c[h=1:data[:H],t=1:data[:T],s=1:Ns], y[:,h,t,s] .<= data[:k_c][:,h] .* sum(y̅[:,τ,s] for τ in 1:t))
    @constraint(model, con_ϱ_e⁻[h=2:data[:H],t=1:data[:T],s=1:Ns], p[:,h,t,s] .- p[:,h-1,t,s] .>= - data[:r_e] .* data[:p̅])
    @constraint(model, con_ϱ_e⁺[h=2:data[:H],t=1:data[:T],s=1:Ns], p[:,h-1,t,s] .- p[:,h,t,s] .>= - data[:r_e] .* data[:p̅])
    @constraint(model, con_ϱ_c⁻[h=2:data[:H],t=1:data[:T],s=1:Ns], y[:,h,t,s] .- y[:,h-1,t,s] .+ data[:r_c] .* sum(y̅[:,τ,s] for τ in 1:t) .>= 0)
    @constraint(model, con_ϱ_c⁺[h=2:data[:H],t=1:data[:T],s=1:Ns], y[:,h-1,t,s] .- y[:,h,t,s] .+ data[:r_c] .* sum(y̅[:,τ,s] for τ in 1:t) .>= 0)
    # non-anticipativity constraints
    @constraint(model, non_ant[t=1:data[:T],s=1:Ns,s_=1:Ns; ST[:L][1,s,t] .== ST[:L][1,s_,t] && s != s_], y̅[:,t,s] .== y̅[:,t,s_])
    # solve model
    optimize!(model)
    @info("done SAA capacity expansion: $(termination_status(model))")
    sol = Dict(
    :obj => JuMP.objective_value.(model),
    :y̅ => JuMP.value.(y̅),
    )
    return sol
end

## functions for post processing
function single_scen_OPF(set,data,ξ,s,y̅_fix,ϑ̅_fix,φ̅_fix)
    # model definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # model variables
    @variable(model, ϑ[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, φ⁺[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, φ⁻[1:data[:dims][:K],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, y[1:data[:dims][:J],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, p[1:data[:dims][:G],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    @variable(model, l[1:data[:dims][:D],1:set[:W],1:set[:H],1:set[:T]] >= 0)
    # model objective
    @objective(model, Min,
                    sum(
                        + (Q(t)*S(t)*ξ[:ξ̂][:,s] + sum(set[:year_per_period]*(O_c(τ)*S(τ)*ξ[:ξ̂][:,s]) for τ in t:set[:T]))'*y̅_fix[:,t]
                        + (Q_ϑ̅(t)*S(t)*ξ[:ξ̂][:,s] + sum(set[:year_per_period]*(O_ϑ̅(τ)*S(τ)*ξ[:ξ̂][:,s]) for τ in t:set[:T]))'*ϑ̅_fix[:,t]
                        + (Q_φ̅(t)*S(t)*ξ[:ξ̂][:,s] + sum(set[:year_per_period]*(O_φ̅(τ)*S(τ)*ξ[:ξ̂][:,s]) for τ in t:set[:T]))'*φ̅_fix[:,t]
                        + set[:year_per_period]*52/sum(data[:load][:w])*sum(data[:load][:w][w]*
                        (
                          168/set[:H]*sum((B(t)*S(t)*ξ[:ξ̂][:,s])'*y[:,w,h,t] for h in 1:set[:H]) +
                          168/set[:H]*sum((V(t)*S(t)*ξ[:ξ̂][:,s])'*p[:,w,h,t] for h in 1:set[:H]) +
                          168/set[:H]*sum(9e3*ones(data[:dims][:D])'*l[:,w,h,t] for h in 1:set[:H])
                        )
                    for w in 1:set[:W])
                    for t in 1:set[:T])
    )
    # power flow equations
    @constraint(model, λ_b[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  ones(data[:dims][:N])'*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:ξ̂][:,s] .- l[:,w,h,t])) .== 0)
    @constraint(model, λ_f̅[w=1:set[:W],h=1:set[:H],t=1:set[:T]], -data[:net][:F]*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:ξ̂][:,s] .- l[:,w,h,t])) .>= -data[:net][:f̅])
    @constraint(model, λ_f̲[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  data[:net][:F]*(data[:gen_e][:M] * p[:,w,h,t] + data[:gen_c][:M] * y[:,w,h,t] + data[:stor][:M] * φ⁻[:,w,h,t] - data[:stor][:M] * φ⁺[:,w,h,t] - data[:load][:M] * (data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:ξ̂][:,s] .- l[:,w,h,t])) .>= -data[:net][:f̅])
    # generation limits
    @constraint(model, λ_p[w=1:set[:W],h=1:set[:H],t=1:set[:T]], -p[:,w,h,t] .>= -data[:gen_e][:k][:,w,h] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_y[w=1:set[:W],h=1:set[:H],t=1:set[:T]], data[:gen_c][:k][:,w,h] .* sum(y̅_fix[:,τ] for τ in 1:t) - y[:,w,h,t] .>= 0)
    @constraint(model, λ_r_e⁻[w=1:set[:W],h=2:set[:H],t=1:set[:T]], p[:,w,h,t] .- p[:,w,h-1,t] .>= - data[:gen_e][:r̅] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_r_e⁺[w=1:set[:W],h=2:set[:H],t=1:set[:T]], p[:,w,h-1,t] .- p[:,w,h,t] .>= - data[:gen_e][:r̅] .* data[:gen_e][:p̅][:,t])
    @constraint(model, λ_r_c⁻[w=1:set[:W],h=2:set[:H],t=1:set[:T]], data[:gen_c][:r̅] .* sum(y̅_fix[:,τ] for τ in 1:t) .+ y[:,w,h,t] .- y[:,w,h-1,t] .>= 0)
    @constraint(model, λ_r_c⁺[w=1:set[:W],h=2:set[:H],t=1:set[:T]], data[:gen_c][:r̅] .* sum(y̅_fix[:,τ] for τ in 1:t) .+ y[:,w,h-1,t] .- y[:,w,h,t] .>= 0)
    # storage limits and equations
    @constraint(model, λ_s[w=1:set[:W],h=2:set[:H],t=1:set[:T]], ϑ[:,w,h,t] .- ϑ[:,w,h-1,t] .- φ⁺[:,w,h,t]*data[:stor][:η⁺] .+ φ⁻[:,w,h,t]/data[:stor][:η⁻] .== 0)
    @constraint(model, λ_s_1[w=1:set[:W],t=1:set[:T]], ϑ[:,w,1,t] .- φ⁺[:,w,1,t]*data[:stor][:η⁺] .+ φ⁻[:,w,1,t]/data[:stor][:η⁻] .== 0)
    @constraint(model, λ_φ⁺[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(φ̅_fix[:,τ] for τ=1:t) .- φ⁺[:,w,h,t] .>= 0)
    @constraint(model, λ_φ⁻[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(φ̅_fix[:,τ] for τ=1:t) .- φ⁻[:,w,h,t] .>= 0)
    @constraint(model, λ_φ[w=1:set[:W],h=1:set[:H],t=1:set[:T]],  sum(φ̅_fix[:,τ] for τ=1:t) .- φ⁺[:,w,h,t] .- φ⁻[:,w,h,t] .>= 0)
    @constraint(model, λ_ϑ[w=1:set[:W],h=1:set[:H],t=1:set[:T]], sum(ϑ̅_fix[:,τ] for τ=1:t) .- ϑ[:,w,h,t] .>= 0)
    # emission limits
    @constraint(model, λ_e[t=1:set[:T]],  - sum(data[:load][:w][w]*168/set[:H]*sum(sum(data[:gen_e][:e].*data[:gen_e][:h].*p[:,w,h,t])
    + sum(data[:gen_c][:e].*data[:gen_c][:h].*y[:,w,h,t]) for h in 1:set[:H]) for w in 1:set[:W])./10^6 >= - (E̅(t)*S(t)*ξ[:ξ̂][:,s])[1])
    # solve model
    optimize!(model)
    @info("done solving a sample of OPF problem: $(termination_status(model))")
    sol = Dict(
    "status" => termination_status(model),
    "obj" => JuMP.objective_value.(model),
    "shed" => sum(JuMP.value.(l)))
    return sol
end
function post_processing(outdir,data,set,ξ)
    # open solution from json files as dictionaries
    sol_sto = Dict()
    open("$(outdir)stochastic.json", "r") do f
        # global sol_sto
        dicttxt = read(f,String)     # file information to string
        sol_sto=JSON.parse(dicttxt)  # parse and transform data
    end
    sol_det = Dict()
    open("$(outdir)determenistic.json", "r") do f
        # local sol_det
        dicttxt = read(f,String)     # file information to string
        sol_det=JSON.parse(dicttxt)  # parse and transform data
    end
    # stochastic solution
    cost_sto = zeros(set[:Ns]);
    inf_sto = zeros(set[:Ns]);
    shed_sto = zeros(set[:Ns]);
    price_sto = zeros(set[:T],set[:W],set[:H],data[:dims][:N],set[:Ns])
    _Y̅_ = zeros(data[:dims][:J]*set[:T],set[:Ns])
    _Θ̅_ = zeros(data[:dims][:K]*set[:T],set[:Ns])
    _Φ̅_ = zeros(data[:dims][:K]*set[:T],set[:Ns])
    for scenario in 1:set[:Ns]
        y̅_fix = max.(0,hcat([Y̅_(t,sol_sto)*S(t)*ξ[:ξ̂][:,scenario] for t in 1:set[:T]]...))
        ϑ̅_fix = max.(0,hcat([Θ̅_(t,sol_sto)*S(t)*ξ[:ξ̂][:,scenario] for t in 1:set[:T]]...))
        φ̅_fix = max.(0,hcat([Φ̅_(t,sol_sto)*S(t)*ξ[:ξ̂][:,scenario] for t in 1:set[:T]]...))
        sol = single_scen_OPF(set,data,ξ,scenario,y̅_fix,ϑ̅_fix,φ̅_fix)
        cost_sto[scenario] = sol["obj"]
        sol["shed"] >= 1 ? inf_sto[scenario] = 1 : NaN
        shed_sto[scenario] = sol["shed"]
        display(shed_sto[scenario])
        _Y̅_[:,scenario] = reshape(y̅_fix, (data[:dims][:J]*set[:T],1))
        _Θ̅_[:,scenario] = reshape(ϑ̅_fix, (data[:dims][:K]*set[:T],1))
        _Φ̅_[:,scenario] = reshape(φ̅_fix, (data[:dims][:K]*set[:T],1))
    end
    # determenistic solution
    cost_det = zeros(set[:Ns]);
    inf_det = zeros(set[:Ns]);
    shed_det = zeros(set[:Ns]);
    price_det = zeros(set[:T],set[:W],set[:H],data[:dims][:N],set[:Ns])
    for scenario in 1:set[:Ns]
        y̅_fix = hcat(sol_det["y̅"]...)
        ϑ̅_fix = hcat(sol_det["ϑ̅"]...)
        φ̅_fix = hcat(sol_det["φ̅"]...)
        sol = single_scen_OPF(set,data,ξ,scenario,y̅_fix,ϑ̅_fix,φ̅_fix)
        cost_det[scenario] = sol["obj"]
        sol["shed"] >= 1 ? inf_det[scenario] = 1 : NaN
        shed_det[scenario] = sol["shed"]
        display(shed_det[scenario])
    end
    # save summary
    open("$(outdir)summary.txt","a") do io
        println(io,"det_ofs_cost = ", mean(cost_det))
        println(io,"sto_ofs_cost = ", mean(cost_sto))
        println(io,"sto_ofs_stdY̅ = ", sum(std(_Y̅_,dims=2)))
        println(io,"sto_ofs_stdΘ̅ = ", sum(std(_Θ̅_,dims=2)))
        println(io,"sto_ofs_stdΦ̅ = ", sum(std(_Φ̅_,dims=2)))
        println(io,"det_ofs_shed = ", mean(shed_det))
        println(io,"sto_ofs_shed = ", mean(shed_sto))
        println(io,"det_ofs_viol = ", mean(sum(inf_det)/set[:Ns]*100))
        println(io,"sto_ofs_viol = ", mean(sum(inf_sto)/set[:Ns]*100))
    end
end
function out_of_sample(Y̅,s)
    # model definition
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # model variables
    @variable(model, y̅[1:data[:J],1:data[:T]] >= 0)
    @variable(model, y[1:data[:J],1:data[:H],1:data[:T]] >= 0)
    @variable(model, p[1:data[:G],1:data[:H],1:data[:T]] >= 0)
    @variable(model, c_e[1:data[:G],1:data[:H],1:data[:T]])
    @variable(model, c_c[1:data[:J],1:data[:H],1:data[:T]])
    # model objective
    @objective(model, Min,
                    sum((Q_g(t)*S(t)*ξ[:ξ̂][:,s])'*y̅[:,t]
                    + sum((B(t)*S(t)*ξ[:ξ̂][:,s])'*y[:,h,t] + sum(c_c[i,h,t] for i in 1:data[:J]) for h in 1:data[:H])
                    + sum((V(t)*S(t)*ξ[:ξ̂][:,s])'*p[:,h,t] + sum(c_e[i,h,t] for i in 1:data[:G]) for h in 1:data[:H])
                    for t in 1:data[:T])
    )
    @constraint(model, cost_e[j=1:data[:G],h=1:data[:H],t=1:data[:T]], [1/2;c_e[j,h,t];sqrt(diag(data[:H̃])[j])*p[j,h,t]] in RotatedSecondOrderCone())
    @constraint(model, cost_c[j=1:data[:J],h=1:data[:H],t=1:data[:T]], [1/2;c_c[j,h,t];sqrt(diag(data[:C])[j])*y[j,h,t]] in RotatedSecondOrderCone())
    # fix investment decision
    @constraint(model, y_fx[t=1:data[:T]], y̅[:,t] .== Y̅[:,t])
    # # power flow eqations
    @constraint(model, con_μ[h=1:data[:H],t=1:data[:T]],  ones(data[:N])'*(data[:M_e] * p[:,h,t] + data[:M_c] * y[:,h,t] - data[:M_l] * (data[:k_l][:,h] .* L(t)*S(t)*ξ[:ξ̂][:,s])) .== 0)
    @constraint(model, con_ζ̅[h=1:data[:H],t=1:data[:T]],  data[:F]*(data[:M_e] * p[:,h,t] + data[:M_c] * y[:,h,t] - data[:M_l] * (data[:k_l][:,h] .* L(t)*S(t)*ξ[:ξ̂][:,s])) .<= data[:f̅])
    @constraint(model, con_ζ̲[h=1:data[:H],t=1:data[:T]], -data[:f̅] .<= data[:F]*(data[:M_e] * p[:,h,t] + data[:M_c] * y[:,h,t] - data[:M_l] * (data[:k_l][:,h] .* L(t)*S(t)*ξ[:ξ̂][:,s])))
    # generation limits
    @constraint(model, con_κ_e[h=1:data[:H],t=1:data[:T]], p[:,h,t] .<= data[:k_e][:,h] .* data[:p̅])
    @constraint(model, con_κ_c[h=1:data[:H],t=1:data[:T]], y[:,h,t] .<= data[:k_c][:,h] .* sum(y̅[:,τ] for τ in 1:t))
    @constraint(model, con_ϱ_e⁻[h=2:data[:H],t=1:data[:T]], p[:,h,t] .- p[:,h-1,t] .>= - data[:r_e] .* data[:p̅])
    @constraint(model, con_ϱ_e⁺[h=2:data[:H],t=1:data[:T]], p[:,h-1,t] .- p[:,h,t] .>= - data[:r_e] .* data[:p̅])
    @constraint(model, con_ϱ_c⁻[h=2:data[:H],t=1:data[:T]], y[:,h,t] .- y[:,h-1,t] .+ data[:r_c] .* sum(y̅[:,τ] for τ in 1:t) .>= 0)
    @constraint(model, con_ϱ_c⁺[h=2:data[:H],t=1:data[:T]], y[:,h-1,t] .- y[:,h,t] .+ data[:r_c] .* sum(y̅[:,τ] for τ in 1:t) .>= 0)
    optimize!(model)
    @info("ofs status in scenatio $(s): $(termination_status(model))")
    return Dict(:status => termination_status(model), :obj => JuMP.objective_value(model))
end
function polot_assumptions(set)
    assumptions = plot(legend=:topleft,ylabel="percentage change", title = "Model assumptions")
    plot!([sum(set[:r_l][1:i]) for i in 1:set[:T]], label="load",lw=2,linestyle=:dash)
    plot!([sum(set[:r_q_s_e][1:i]) for i in 1:set[:T]], label="storage_capex",lw=2)
    plot!([sum(set[:r_q_wind][1:i]) for i in 1:set[:T]], label="wind_capex",lw=2)
    plot!([sum(set[:r_q_pv][1:i]) for i in 1:set[:T]], label="pv_capex",lw=2)
    plot!([sum(set[:r_q_nucl][1:i]) for i in 1:set[:T]], label="nucl_capex",lw=2)
    plot!([sum(set[:r_q_ng][1:i]) for i in 1:set[:T]], label="ng_capex",lw=2)
    plot!([sum(set[:r_ng][1:i]) for i in 1:set[:T]], label="ng_price",lw=2, linestyle=:dashdotdot)
    plot!([sum(set[:r_cl][1:i]) for i in 1:set[:T]], label="coal_price",lw=2, linestyle=:dashdotdot)
    plot!(xticks=(1:set[:T], ["2025","2030","2035","2040","2045","2050"]))
end
function plot_generation(sol)
    wind_gen_can = zeros(set[:T]); sola_gen_can = zeros(set[:T]); ccgt_gen_can = zeros(set[:T]); ccgt_ccs_can = zeros(set[:T])
    nucl_gen_can = zeros(set[:T]); hydr_gen_can = zeros(set[:T]); coal_gen_can = zeros(set[:T])

    wind_gen_exs = zeros(set[:T]); sola_gen_exs = zeros(set[:T]); ccgt_gen_exs = zeros(set[:T]); ccgt_ccs_exs = zeros(set[:T])
    nucl_gen_exs = zeros(set[:T]); hydr_gen_exs = zeros(set[:T]); coal_gen_exs = zeros(set[:T])

    for t in 1:set[:T]
        wind_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("Wind"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        sola_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("PV"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        ccgt_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("GT"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        ccgt_ccs_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("NG"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        nucl_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("Nucl"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        hydr_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("Hydro"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        coal_gen_can[t] = sum(168/set[:H]*sum(sol["y"][ind_can("Coal"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6

        wind_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("Wind"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        sola_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("PV"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        ccgt_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("GT"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        ccgt_ccs_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("NG"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        nucl_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("Nucl"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        hydr_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("Hydro"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
        coal_gen_exs[t] = sum(168/set[:H]*sum(sol["p"][ind_exs("Coal"),w,:,t])*data[:load][:w][w] for w in 1:set[:W])/10^6
    end
    generation = groupedbar([wind_gen_exs wind_gen_can sola_gen_exs sola_gen_can ccgt_gen_exs ccgt_gen_can ccgt_ccs_exs ccgt_ccs_can hydr_gen_exs hydr_gen_can coal_gen_exs coal_gen_can nucl_gen_exs nucl_gen_can],
            title = "Annual generation [TWh]",
            bar_position = :stack,
            legend=:topleft,
            bar_width=0.7,
            xticks=(1:set[:T], ["2025","2030","2035","2040","2045","2050"]),
            seriescolor = [:deepskyblue2 :deepskyblue2 :yellow3 :yellow3 :maroon :maroon :navajowhite3 :navajowhite3 :aqua :aqua :gray :gray :salmon :salmon],
            label=["wind" false "solar" false "ccgt" false "ccgt_ccs" false "hydro" false "coal" false "nuclear" false],
            fillalpha = [0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0])
    return generation
end
function plot_capacity(sol,plottitle)
    wind_cap_can = zeros(set[:T]); sola_cap_can = zeros(set[:T]); ccgt_cap_can = zeros(set[:T]); ccgt_ccs_can = zeros(set[:T])
    nucl_cap_can = zeros(set[:T]); hydr_cap_can = zeros(set[:T]); coal_cap_can = zeros(set[:T])
    stor_cap_can = zeros(set[:T]);

    wind_cap_exs = zeros(set[:T]); sola_cap_exs = zeros(set[:T]); ccgt_cap_exs = zeros(set[:T]); ccgt_ccs_exs = zeros(set[:T])
    nucl_cap_exs = zeros(set[:T]); hydr_cap_exs = zeros(set[:T]); coal_cap_exs = zeros(set[:T])
    stor_cap_exs = zeros(set[:T]);

    for t in 1:set[:T]
        wind_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("Wind"),t])/1000
        sola_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("PV"),t])/1000
        ccgt_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("GT"),t])/1000
        ccgt_ccs_exs[t] = sum(data[:gen_e][:p̅][ind_exs("NG"),t])/1000
        nucl_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("Nuc"),t])/1000
        hydr_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("Hydro"),t])/1000
        coal_cap_exs[t] = sum(data[:gen_e][:p̅][ind_exs("Coal"),t])/1000

        wind_cap_can[t] = sum(sol["y̅"][ind_can("Wind"),1:t])/1000
        sola_cap_can[t] = sum(sol["y̅"][ind_can("PV"),1:t])/1000
        ccgt_cap_can[t] = sum(sol["y̅"][ind_can("GT"),1:t])/1000
        ccgt_ccs_can[t] = sum(sol["y̅"][ind_can("NG"),1:t])/1000
        nucl_cap_can[t] = sum(sol["y̅"][ind_can("Nuc"),1:t])/1000
        hydr_cap_can[t] = sum(sol["y̅"][ind_can("Hydro"),1:t])/1000
        stor_cap_can[t] = sol["ϑ̅"] == 0 ? 0 : sum(sol["ϑ̅"][:,1:t])/1000
    end
    capacity = groupedbar([stor_cap_exs stor_cap_can wind_cap_exs wind_cap_can sola_cap_exs sola_cap_can ccgt_cap_exs ccgt_cap_can ccgt_ccs_exs ccgt_ccs_can hydr_cap_exs hydr_cap_can coal_cap_exs coal_cap_can nucl_cap_exs nucl_cap_can],
            title = plottitle,
            framestyle = :box,
            bar_position = :stack,
            legend=:topleft,
            bar_width=0.7,
            ylim=(0,800),
            xticks=(1:set[:T], ["2025","2030","2035","2040","2045"]),
            seriescolor = [:forestgreen :forestgreen :deepskyblue2 :deepskyblue2 :yellow3 :yellow3 :maroon :maroon :navajowhite3 :navajowhite3 :aqua :aqua :gray :gray :salmon :salmon],
            label=["storage" false "wind" false "solar" false "ccgt" false "ccgt_ccs" false "hydro" false "coal" false "nuclear" false],
            # fillstyle = [:x :/ :x :/ :x :/ :x :/ :x :/ :x :/ :x :/ :x :/]
            fillalpha = [0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0 0.75 1.0]
            )
    return capacity
end
function plot_week(t,w,sol)
    # load
    load = [sum(data[:load][:k][:,t,w,h] .* L(t)*S(t)*ξ[:μ̂]) for h in 1:set[:H]]

    # generation
    wind = sum(sol["y"][ind_can("Wind"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("Wind"),w,:,t], dims=1)
    sola = sum(sol["y"][ind_can("PV"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("PV"),w,:,t], dims=1)
    nucl = sum(sol["y"][ind_can("Nuc"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("Nuc"),w,:,t], dims=1)
    ccgt = sum(sol["y"][ind_can("GT"),w,:,t], dims=1) .+ sum(sol["y"][ind_can("NG"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("GT"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("NG"),w,:,t], dims=1)
    hydr = sum(sol["y"][ind_can("Hydro"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("Hydro"),w,:,t], dims=1)
    coal = sum(sol["y"][ind_can("Coal"),w,:,t], dims=1) .+ sum(sol["p"][ind_exs("Coal"),w,:,t], dims=1)

    # storage
    stor = sum(sol[:φ][:,w,:,t], dims=1)
    charg_ind = findall(x->x>=0,stor)
    disch_ind = findall(x->x<0,stor)
    stor_charg = deepcopy(stor); stor_charg[disch_ind] .= 0
    stor_disch = deepcopy(stor); stor_disch[charg_ind] .= 0
    Y = [nucl;coal;hydr;ccgt;sola;wind;abs.(stor_disch)]'./1000

    plo=plot(legend=false)
    areaplot!(1:set[:H], Y,
    seriescolor = [:salmon :gray :aqua :maroon :yellow3 :deepskyblue2 :forestgreen],
    label=["nuclear" "coal" "hydro" "ccgt" "solar" "wind" "storage"],legend=:outerright, legend_columns=-2)
    areaplot!(1:set[:H], [-stor_charg;zeros(set[:H])']'./1000,
    seriescolor = [:forestgreen],
    label=false,legend=:outerright, legend_columns=-1)

    plot!(load./1000, c=:red,lw=2,label="load")

    plot!(title="stage $(t) week $(w)")

    return plo
end
