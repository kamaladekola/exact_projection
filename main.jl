include("Data.jl")
include("Domain.jl")
include("aravena_FB_EP.jl")
include("get_ATC.jl")
include("aravena_ATC_EP.jl")

using .Data, .FB_EP, .get_ATC, .ATC_EP
using .Domain: plot_domains

function main()
    println("→ Running FB EP…")
    fb_model, obj, dayahead, redispatch, line_flows, netpos, exchange = FB_ep()

    println("→ Computing ATC…")
    atc = get_atc()
    println("ATC: ", atc)

    atc_dict = Dict(zip(Data.TCONNECT, atc))
    println("ATC: ", atc_dict)

    println("→ Running ATC-EP…")
    atc_model, obj_atc, dayahead_atc, redispatch_atc, flows_atc, netpos_atc, exchange_atc = atc_ep(atc)
    
    println("→ plotting domain…")
    market_clearing = (netpos[:A], netpos[:B])
    atc_clearing = (netpos_atc[:A], netpos_atc[:B])
    plt = plot_domains(atc_dict, market_clearing, atc_clearing)
    display(plt)
    # savefig(plt, "flow_based_domain.png")

    println("done.")
end

main()