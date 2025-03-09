using QuantEcon
include("src/utils.jl")


function stationary(Pi; pi_seed=nothing, tol=1E-27, maxit=100000)
    """Find invariant distribution of a Markov chain by iteration."""
    if pi_seed == nothing
        pi = ones(size(Pi)[1]) / size(Pi)[1]
    else
        pi = pi_seed
    end
    
    pi_new = copy(pi)
    converged=true
    for it in 1:maxit
        pi_new = Pi' * pi
        if maximum(abs.(pi_new .- pi)) < tol
            converged=false
            break
        end
        pi = copy(pi_new)
    end
    if converged
        throw("No convergence after $maxit forward iterations!")
    end
    pi = pi_new

    return pi
end

function setmin!(A, floor_value)
    @. A = max(A, floor_value)
end

function interpolate_f(x,y,xq) 
    if x[2] < x[1]
        x = reverse(x)
        y = reverse(y)
    end
    func = extrapolate( interpolate((x,), y, Gridded(Linear())), Line()) 
    return func.(xq)
end

function backward_core(Va_p, a_grid, y, pr)
    c_nextgrid = Va_p
    coh = (1 + pr.r) .* a_grid' .+ y
    
    # print( size(a) )
    a_endo = similar(coh)
    @threads for i in 1:size(c_nextgrid)[1]
        a_endo[i, :] = interpolate_f((c_nextgrid .+ a_grid')[i,:], a_grid,coh[i,:])
    end
    
    setmin!(a_endo, a_grid[1])
    c = coh - a_endo
    Va = (1 + pr.r) * c .^ ( -pr.σ)
    return a_endo, c
end

function enforce_constraints(a, abar, bbar, total)
    """Given a proposed split of 'total' into 'a' and residual 'b', enforce constraint
    that a>=abar and b>bbar (assume this is always possible, but prioritize b)"""
    a = max.(a, abar)
    b = total' .- a
    b = max.(b, bbar)
    a = total' .- b
    return a, b
end


function get_Psi_and_deriv(ap, a, pr)
    """Adjustment cost Psi(ap, a) and its derivatives with respect to
    first argument (ap) and second argument (a)"""
    a_with_return = (1 + pr.r) * a
    a_change = ap .- a_with_return
    abs_a_change = abs.(a_change)
    sign_change = sign.(a_change)

    adj_denominator = a_with_return .+ pr.χ0
    core_factor = (abs_a_change ./ adj_denominator) .^ (pr.χ2 - 1)

    Psi = (pr.χ1 / pr.χ2) .* abs_a_change .* core_factor
    Psi1 = pr.χ1 .* sign_change .* core_factor
    Psi2 = -(1 + pr.r) .* (Psi1 .+ (pr.χ2 - 1) .* Psi ./ adj_denominator)
    return Psi, Psi1, Psi2
end

function lhs_equals_rhs_interpolate(lhs::Array{Float64, 3}, rhs::Array{Float64, 2})::Tuple{Array{Int64, 3}, Array{Float64, 3}}
    pi        = similar(lhs)
    a_end_idx = Int.(floor.(lhs))
    for z in 1:size(lhs)[1]
        for i in 1:size(lhs)[2]
            for j in 1:size(rhs)[2]
                indx = searchsortedlast(lhs[z,i,:] .- rhs[:,j], 0, lt= >) # find idx of first a that is bigger than true a
                if indx ==0
                    a_end_idx[z,i,j] = max(indx,1)
                    pi[z,i,j] = 1
                elseif indx==size(rhs)[2]
                    a_end_idx[z,i,j] = min(size(rhs)[2]-1,indx)
                    pi[z,i,j] = 0
                else
                    err_upper = rhs[indx+1, j] - lhs[z,i,indx+1]
                    err_lower = rhs[indx, j]   - lhs[z,i,indx]
                    pi[z,i,j] = err_upper / (err_upper - err_lower)
                    a_end_idx[z,i,j] = indx
                end 
            end
        end
    end
    return a_end_idx, pi
end

function lhs_equals_rhs_interpolate_new(lhs, rhs)
    pi        = similar(lhs)
    a_end_idx = Int.(floor.(lhs))
    Φ = kron(Matrix(I,pr.nz,pr.nz),BasisMatrix(abasis,Direct()).vals[1])
    repeat(rhs', )
    
end


function apply_coord(x_i, x_pi, y)
    """Use representation xqi, xqpi to get yq at xq:
    yq = xqpi * y[xqi] + (1-xqpi) * y[xqi+1]

    Parameters
    ----------
    xqi  : array (nq), indices of lower bracketing gridpoints
    xqpi : array (nq), weights on lower bracketing gridpoints
    y  : array (n), data points

    Returns
    ----------
    yq : array (nq), interpolated points
    """
    yq = similar(x_pi)
    for iq in 1:length(x_i)
        y_low = y[x_i[iq]]
        y_high = y[x_i[iq]+1]
        yq[iq] = x_pi[iq]*y_low + (1-x_pi[iq])*y_high
    end
    return yq
end

function interpolate_coord(x, xq)
    """Get representation xqi, xqpi of xq interpolated against x:
    xq = xqpi * x[xqi] + (1-xqpi) * x[xqi+1]

    Parameters
    ----------
    x    : array (n), ascending data points
    xq   : array (nq), ascending query points

    Returns
    ----------
    i  : array (nq), indices of lower bracketing gridpoints
    pi : array (nq), weights on lower bracketing gridpoints
    """
    pi  = similar(x)
    idx = Int.(floor.(x))
    for i in 1:length(x)
        indx = searchsortedlast(x, xq[i], lt= <) # find idx of first a that is bigger than true a
        # print(indx,"_",i,"_____")
        if indx == 0 
            idx[i] = max(indx,1)
            pi[i] = 1
        elseif indx == length(x)
            idx[i] = min(length(x)-1,indx)
            pi[i] = 0
        else
            pi[i] = (x[indx+1] - xq[i]) / (x[indx+1] - x[indx])
            idx[i] = indx
        end
    end
    return idx, pi
end