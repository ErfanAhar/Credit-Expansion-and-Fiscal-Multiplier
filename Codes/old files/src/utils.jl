
function grid_maker(pr)
    # b = exp.(exp.(range(0,2,length=pr.nb)) .- 1) .- (1-pr.min_b) 
     # a = exp.(exp.(range(0,2,length=pr.na)) .- 1) .- (1-pr.min_a) 
     # b = exp.(range(0,7,length=pr.nb))  .- (1-pr.min_b) 
     # a = exp.(range(0,7,length=pr.na))  .- (1-pr.min_a) 
     pivot = 0.25
     a = exp10.(range(log10(pivot), log10(pr.max_a - pr.min_a + pivot), length=pr.na) ) .+ pr.min_a .- pivot
     a[1]   = pr.min_a
     a[end] = pr.max_a
     b = exp10.(range(log10(pivot), log10(pr.max_b - pr.min_b + pivot), length=pr.nb) ) .+ pr.min_b .- pivot
     b[1]   = pr.min_b
     b[end] = pr.max_b
     k = exp10.(range(log10(pivot), log10(pr.max_k - pr.min_k + pivot), length=pr.nk) ) .+ pr.min_k .- pivot
     k[1]   = pr.min_k
     k[end] = pr.max_k
    z = rouwenhorst(pr.nz, pr.ρe, pr.σe, pr.μe).state_values
    z_p = rouwenhorst(pr.nz, pr.ρe, pr.σe, pr.μe).p
    return a,b,k,z,z_p
end

function grid_maker_new(pr)
    # b = exp.(exp.(range(0,2,length=pr.nb)) .- 1) .- (1-pr.min_b) 
     # a = exp.(exp.(range(0,2,length=pr.na)) .- 1) .- (1-pr.min_a) 
     # b = exp.(range(0,7,length=pr.nb))  .- (1-pr.min_b) 
     # a = exp.(range(0,7,length=pr.na))  .- (1-pr.min_a) 
     pivot = 0.25
     a = exp10.(range(log10(pivot), log10(pr.max_a - pr.min_a + pivot), length=pr.na) ) .+ pr.min_a .- pivot
     a[1]   = pr.min_a
     a[end] = pr.max_a
     b = exp10.(range(log10(pivot), log10(pr.max_b - pr.min_b + pivot), length=pr.nb) ) .+ pr.min_b .- pivot
     b[1]   = pr.min_b
     b[end] = pr.max_b
     k = exp10.(range(log10(pivot), log10(pr.max_k - pr.min_k + pivot), length=pr.nk) ) .+ pr.min_k .- pivot
     k[1]   = pr.min_k
     k[end] = pr.max_k
     
    function grid(min_g, max_g, n, θ)
        return range(0,1, length=n).^θ .* (max_g - min_g) .+ min_g
    end

    a = grid(pr.min_a, pr.max_a, pr.na, pr.aθ)
    b = grid(pr.min_b, pr.max_b, pr.nb, pr.bθ)
    k = grid(pr.min_k, pr.max_k, pr.nk, pr.bθ)

    z = rouwenhorst(pr.nz, pr.ρe, pr.σe, pr.μe).state_values
    z_p = rouwenhorst(pr.nz, pr.ρe, pr.σe, pr.μe).p
    return a,b,k,z,z_p
end
