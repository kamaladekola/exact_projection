
module Domain

using LazySets, Plots
import ..Data: PTDF, zonal_PTDF, BRANCH_LIMIT

export plot_domains

"""
    fb_polyhedron()

Return the flow-based feasibility polyhedron in (NP_A, NP_B) space.
"""

function fb_polyhedron()
    H = HalfSpace[]

    for ℓ in axes(zonal_PTDF, 1)             
        ram = BRANCH_LIMIT[ℓ]

        pA = zonal_PTDF[ℓ,1]                  # zone A
        pB = zonal_PTDF[ℓ,2]                  # zone B
        pC = zonal_PTDF[ℓ,3]                  # zone C

        Acoef = pA - pC                      # (PTDF_A – PTDF_C)
        Bcoef = pB - pC                      # (PTDF_B – PTDF_C)

        push!(H, HalfSpace([ Acoef,  Bcoef],  ram))
        push!(H, HalfSpace([-Acoef, -Bcoef],  ram))
    end
    
    return HPolyhedron(H)
end


"""
    atc_polyhedron(atc)

Return the 2‑D polyhedron of feasible net positions that respect the
bilateral ATC limits of the interconnectors

"""
function atc_polyhedron(atc::Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}})

    # ATC limits
    export_AC, import_AC = atc[(:A,:C)]   # A ↔ C
    export_CB, import_CB = atc[(:C,:B)]   # C ↔ B
    export_BA, import_BA = atc[(:B,:A)]   # B ↔ A

    #   NP_A =  e_AC  −  e_BA
    #   NP_B =  e_BA  −  e_CB
    #   NP_C = −(NP_A + NP_B)

    H = HalfSpace[]

    # -----  NP_A = e_AC − e_BA -----
    # max export from A: e_AC ≤ export_AC
    # max import into A: e_BA ≥ import_BA
    # ⇒ NP_A ≤ export_AC - import_BA
    # ⇒ NP_A ≥ import_AC - export_BA
    push!(H, HalfSpace([+1.0,  0.0],  export_AC - import_BA))  # upper bound
    push!(H, HalfSpace([-1.0,  0.0], -import_AC + export_BA))  # lower bound

    # ----- NP_B = e_BA − e_CB -----
    # max export from B: e_BA ≤ export_BA
    # max import into B: e_CB ≥ import_CB
    push!(H, HalfSpace([ 0.0, +1.0],  export_BA - import_CB))  # upper bound
    push!(H, HalfSpace([ 0.0, -1.0], -import_BA + export_CB))  # lower bound

    # ----- bounds that couple NP_A and NP_B (border C↔B)
    # e_CB ≤ export_CB   and   e_CB ≥ import_CB
    # e_CB = -NP_A - NP_B + e_AC
    # ⇒ NP_A + NP_B ≤ export_AC - import_CB
    # push!(H, HalfSpace([+1.0, +1.0],  export_AC - import_CB))  # NP_A+NP_B upper
    # push!(H, HalfSpace([-1.0, -1.0], -import_AC + export_CB))  # NP_A+NP_B lower

    return HPolyhedron(H)
end

"""
    plot_domains(P_fb, atc_dict, clearing, atc_clearing)

Overlay flow-based and ATC domains, and market clearing points.
"""
function plot_domains(
    atc_dict::Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}},
    clearing::Tuple{Float64,Float64},
    atc_clearing::Tuple{Float64,Float64},
)
    P_fb = fb_polyhedron()
    # Build ATC region
    P_atc = atc_polyhedron(atc_dict)


    plt = plot(P_fb; lw=1.5, alpha=0.25,
               label="flow-based domain",
               xlabel="net position A [MW]",
               ylabel="net position B [MW]",
               ratio=:equal, legend=:topright)

    plot!(P_atc; lw=1.5, alpha=0.25, label="ATC domain")
    scatter!([clearing[1]], [clearing[2]];
            m=:star5, ms=8, mc=:royalblue, msc=:white, msw=1.5,
            label="FBMC")
    scatter!([atc_clearing[1]], [atc_clearing[2]];
            m=:diamond,      # marker shape
            ms=5,           # marker size 
            mc=:tomato,      # marker color
            msc=:white,     # marker stroke color (outline)
            msw=1.5,        # marker stroke width
            label="ATCMC"  # legend label
            )

    return plt
end

end # module Domain