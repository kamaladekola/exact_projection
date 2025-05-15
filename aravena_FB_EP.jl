module FB_EP
# include("Data.jl")
using JuMP, Gurobi
import ..Data: ZONES, NODES, zone_of, GENERATORS, COST, DEMAND, BRANCHES, PTDF, TCONNECT

export FB_ep

function FB_ep()

    m = Model(Gurobi.Optimizer)
    set_silent(m)

    # dayhead dispatch
    @variable(m, 0 <= v_dayahead[i=1:length(GENERATORS)] <= 1)

    # redispatch
    # @variable(m, 0 <= v_redispatch[i=1:length(GENERATORS)] <= 1)
    @variable(m, -1 <= v_redispatch[i=1:length(GENERATORS)] <= 1) # new
    # line flows
    @variable(m, f[ℓ=1:length(BRANCHES)])

    # net position
    @variable(m, net_pos[z in ZONES])

    # cross‐border exchanges
    @variable(m, e[t in TCONNECT])
    for (ℓ,(i,j,F)) in enumerate(BRANCHES)
        # if (i,j) is A–B, f[ℓ] == e[(:A,:B)], etc.
        if (zone_of[i],zone_of[j]) in TCONNECT
            @constraint(m, f[ℓ] ==  e[(zone_of[i],zone_of[j])])
        elseif (zone_of[j],zone_of[i]) in TCONNECT
            @constraint(m, f[ℓ] == -e[(zone_of[j],zone_of[i])])
        end
        # thermal limits on every line
        @constraint(m, -F <= f[ℓ] <= F)
    end
    # netposition = sum over exchanges
    @constraint(m, [z in ZONES],
    net_pos[z] ==
        sum(e[t] for t in TCONNECT if t[1] == z)   # exports
      - sum(e[t] for t in TCONNECT if t[2] == z))  # imports

    # zone balance
    for z in ZONES
        @constraint(m, sum(GENERATORS[i][2]*v_dayahead[i] for i in 1:length(GENERATORS) if zone_of[GENERATORS[i][1]] == z) - net_pos[z]
            == sum(DEMAND[n] for n in NODES if zone_of[n]==z)
        )

        @constraint(m, sum(GENERATORS[i][2] * v_redispatch[i] for i in eachindex(GENERATORS) if zone_of[GENERATORS[i][1]] == z) == 0) # new
    end

    # nodal balance
    for ℓ in eachindex(BRANCHES)
        @constraint(m, f[ℓ] == sum(PTDF[ℓ,j] * (sum(GENERATORS[i][2]*(v_dayahead[i] + v_redispatch[i]) for i in eachindex(GENERATORS)
        if GENERATORS[i][1] == NODES[j]) - DEMAND[NODES[j]])
        for j in eachindex(NODES))
        )
    end
    
 
    @constraint(m, [i in eachindex(GENERATORS)], 0 <= v_dayahead[i] + v_redispatch[i] <= 1) # new


    @objective(m, Min, sum( COST[i] * GENERATORS[i][2] * v_dayahead[i] for i in eachindex(GENERATORS)))

    optimize!(m)
    
    # results
    if termination_status(m) in (MOI.OPTIMAL, MOI.FEASIBLE_POINT)
        println("→ Objective value: ", objective_value(m))
        println("\n→ Generator dispatch v[i]:")
        for (i,(node,cap)) in enumerate(GENERATORS)
            # println("   v[$i] @ $(node)  = ", value(v[i]))
            println("   v_dayhead[$i] @ $(node)  = ", value(v_dayahead[i]))
            println("   v_redispatch[$i] @ $(node)  = ", value(v_redispatch[i]))
        end
    
        println("\n→ Line flows f[ℓ]:")
        for (ℓ,(i,j,F)) in enumerate(BRANCHES)
            println("   f[$ℓ] ($(i)→$(j)) = ", value(f[ℓ]))
        end
    
        println("\n→ Net positions net_pos[z]:")
        for z in ZONES
            println("   net_pos[$z] = ", value(net_pos[z]))
        end
    
        println("\n→ Cross-border exchanges e[t]:")
        for t in TCONNECT
            println("   e[$(t[1])→$(t[2])] = ", value(e[t]))
        end

        # collect in variables
        obj          = objective_value(m)
        dayahead  = [value(v_dayahead[i])   for i in eachindex(GENERATORS)]
        redispatch = [value(v_redispatch[i]) for i in eachindex(GENERATORS)]
        line_flows           = [value(f[ℓ])            for ℓ in eachindex(BRANCHES)]
        netpos     = Dict(z => value(net_pos[z]) for z in ZONES)
        exchange           = Dict(t => value(e[t])       for t in TCONNECT)

        return m, obj, dayahead, redispatch, line_flows, netpos, exchange
       
    else
        println("Model did not solve to optimality. Status: ", termination_status(m))
        return m, nothing, nothing, nothing, nothing, nothing, nothing
    end
end

# FB_ep() # for debugging
end # module FB_EP