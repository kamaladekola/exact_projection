
module get_ATC

using JuMP, Gurobi
import ..Data: ZONES, NODES, zone_of, GENERATORS, COST, DEMAND, BRANCHES, PTDF, TCONNECT, BRANCH_LIMIT

export get_atc

function get_atc()
    m = Model(Gurobi.Optimizer)
    # m = Model(Ipopt.Optimizer)
    set_silent(m)

    # ATC variables - two for  each interconnector

    @variable(m, atc_plus[t in TCONNECT])
    @variable(m, atc_minus[t in TCONNECT])


    # 2^3 vertices:
    num_INT = length(TCONNECT)
    signs = collect(Iterators.product(fill([-1,1], num_INT)...))
    verts = length(signs)

    @variable(m, 0 <= v_dayahead[v=1:verts, i=1:length(GENERATORS)] <= 1)     # utilisation factor of generator i; dayhead
    @variable(m, -1  <= v_redispatch[v=1:verts, i=1:length(GENERATORS)] <= 1)   # redispatch

    @variable(m, f[v=1:verts, ℓ=1:length(BRANCHES)])  # realised flow (MW) on physical line ℓ at vertex v.
    @variable(m, net_pos[v=1:verts, z in ZONES])      # zonal net position (MW) at vertex v
    @variable(m, e[v=1:verts, t in TCONNECT])         # cross-border exchange (MW) at vertex v

    @constraint(m, [t in TCONNECT], -atc_minus[t] <= atc_plus[t])

    for v in 1:verts, i in eachindex(GENERATORS)
        @constraint(m, v_dayahead[v,i] + v_redispatch[v,i] >= 0)
        @constraint(m, v_dayahead[v,i] + v_redispatch[v,i] <= 1)
    end

    # For each vertex
    for (v,s) in enumerate(signs)
        @constraint(m, [t in TCONNECT], e[v,t] == -atc_minus[t])
        @constraint(m, [t in TCONNECT], e[v,t] == atc_plus[t])

        for (k,t) in enumerate(TCONNECT)
            @constraint(m, e[v,t] == (s[k] == 1 ?  atc_plus[t] : -atc_minus[t]))
            cap = sum(ℓ[3] for ℓ in BRANCHES if (zone_of[ℓ[1]],zone_of[ℓ[2]]) == t)
            @constraint(m, -cap <= e[v,t] <= cap)
        end

        # nodal balance
        for (ℓ,(_,_,F)) in enumerate(BRANCHES)
            @constraint(m, f[v,ℓ] == sum(PTDF[ℓ,j] * (sum(GENERATORS[i][2] * (v_dayahead[v,i] + v_redispatch[v,i]) for i in eachindex(GENERATORS)
            if GENERATORS[i][1] == NODES[j]) - DEMAND[NODES[j]])
            for j in eachindex(NODES))
            )

            @constraint(m, -F <= f[v,ℓ] <= F) # Thermal limits
        end

        for t in TCONNECT
            @constraint(m, e[v,t] == sum(f[v,ℓ] for (ℓ,(i,j,_)) in enumerate(BRANCHES)
                            if (zone_of[i], zone_of[j]) == t)
                    - sum(f[v,ℓ] for (ℓ,(i,j,_)) in enumerate(BRANCHES)
                            if (zone_of[j], zone_of[i]) == t))
        end


        # zone balance
        for z in ZONES
                    # netposition = sum over exchanges
            @constraint(m, net_pos[v,z] == sum(e[v,t] for t in TCONNECT if t[1] == z) -  # exports
                                        sum(e[v,t] for t in TCONNECT if t[2] == z))  # imports

            @constraint(m, sum(GENERATORS[i][2]*v_dayahead[v,i] for i in eachindex(GENERATORS) if zone_of[GENERATORS[i][1]] == z) - net_pos[v,z]
                == sum(DEMAND[n] for n in NODES if zone_of[n]==z)
            )
                                        
            @constraint(m, sum(GENERATORS[i][2] * v_redispatch[v,i] for i in eachindex(GENERATORS) if zone_of[GENERATORS[i][1]] == z) == 0)                        
        end



    end

    @objective(m, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT))
    # @NLobjective(m, Max, prod(atc_plus[t] + atc_minus[t] for t in TCONNECT))


    optimize!(m)


    ATCplus  = value.(atc_plus)
    ATCminus = value.(atc_minus)
    
    println("=== ATC-EP results ===")
    for t in TCONNECT
      println("  ", t,
              "  ", round(ATCplus[t], digits=1),
              "   ", round(ATCminus[t], digits=1))
    end

     # same order as TCONNECT
        atc_vector = [ (abs(ATCplus[t]), -1* abs(ATCminus[t])) for t in TCONNECT ]
        # Return the ATC values
        return atc_vector
end
# get_atc() # for debugging

end