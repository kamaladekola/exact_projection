module Data

export ZONES, NODES, zone_of, GENERATORS, COST, DEMAND,
       BRANCHES, PTDF, PTDF_zone, TCONNECT, BRANCH_LIMIT

const ZONES  = [:A, :B, :C]
const NODES  = [:n1, :n2, :n3, :n4]
const zone_of = Dict(:n1=>:A, :n2=>:A, :n3=>:B, :n4=>:C)

# generator capacities
const GENERATORS = [
    (:n1, 300.0),  # gen at node n1 (zone A)
    (:n2, 200.0),  # n2 (A)
    (:n3, 200.0),  # n3 (B)
    (:n4, 200.0)   # n4 (C)
]

const COST = [5.0, 70.0, 15.0, 40.0]    # €/MWh for gens 1–4


# fixed demand at each node
const DEMAND = Dict(:n1=>0.0, :n2=>300.0, :n3=>0.0, :n4=>300.0)

# branches = (from, to, thermal limit F [MW])
const BRANCHES = [
    (:n1, :n2, 100.0),   # ℓ₁₂   (A–A) 
    (:n2, :n4, 200.0),   # ℓ₂₄   (A–C)
    (:n4, :n3, 200.0),   # ℓ₄₃   (C-B)
    (:n3, :n1, 200.0)    # ℓ₃₁   (B–A)
]

# node‐to‐line PTDF[ℓ,n]
const PTDF = [
     0.5  -0.25   0.25   0.0;   # l1: 1–2
     0.5   0.75   0.25   0.0;   # l2: 2–4
    -0.5  -0.25  -0.75   0.0;   # l3: 4–3
    -0.5  -0.25   0.25   0.0    # l4: 3–1
]

# Zone-to-line PTDFs
const zonal_PTDF = [
     0.5   0.25   0.0;   # l1: 1–2
     0.5   0.25   0.0;   # l2: 2–4
    -0.5  -0.75   0.0;   # l3: 4–3
    -0.5   0.25   0.0    # l4: 3–1
]
# Interconnectors
const TCONNECT = [(:A,:C), (:C,:B), (:B,:A)]

const BRANCH_LIMIT = Dict(ℓ => F for (ℓ,(_,_,F)) in enumerate(BRANCHES))

end
