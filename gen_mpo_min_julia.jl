using LinearAlgebra
using WignerSymbols
using ITensors
using ITensorMPOConstruction


↓ = false
↑ = true

function create_op_cache_vec(sites::Vector{<:Index})::OpCacheVec
    ## The identity must be the first operator.
    ## The next six operators are the only ones that are used directly in MPO construction.
    operatorNames = ["I", "Cdn", "Cup", "Cdagdn", "Cdagup", "Ndn", "Nup",
                     "Cup * Cdn", "Cup * Cdagdn", "Cup * Ndn",
                     "Cdagup * Cdn", "Cdagup * Cdagdn", "Cdagup * Ndn",
                     "Nup * Cdn", "Nup * Cdagdn", "Nup * Ndn"]
  
    return [[OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)]
end

opC(j::CartesianIndex, spin::Bool, mapping::Array{Int}) = OpID(2 + spin, mapping[wrap_around(size(mapping), j)])
opC(j::Int, spin::Bool) = OpID(2 + spin, j)

opCdag(j::CartesianIndex, spin::Bool, mapping::Array{Int}) = OpID(4 + spin, mapping[wrap_around(size(mapping), j)])
opCdag(j::Int, spin::Bool) = OpID(4 + spin, j)

function calc_int_ps(nm, ps_pot :: Vector{Float64})
    ## N = 2s+1, get V[i,j,k,l]
    int_el = zeros(nm, nm, nm)
    s = .5 * (nm - 1)
    for m1 in 1 : nm 
        m1r = m1 - s - 1
        for m2 in 1 : nm 
            m2r = m2 - s - 1
            for m3 = 1 : nm
                m3r = m3 - s - 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm)
                    continue
                end
                m4r = m4 - s - 1
                for l in 1 : length(ps_pot)
                    if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l)
                        break
                    end 
                    int_el[m1, m2, m3] += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
                end 
            end 
        end 
    end
    return int_el
end

function ITensors.space( :: SiteType"Electron" ; lz :: Int = 1, 
    conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_spin_parity = true,
    qnname_sz = "Sz", qnname_nf = "Nf", qnname_lz = "Lz", qnname_spin_parity = "Z2")
    ## Note that conserve nf and lz is always set on. Lz is TWICE the actual angular momentum
    if (conserve_sz || (!conserve_nf)) return 4 end
    if (conserve_spin_parity && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,      0), (qnname_spin_parity, 0, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz,     lz), (qnname_spin_parity, 1, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz,     lz), (qnname_spin_parity, 0, 2)) => 1
            QN((qnname_nf, 2, -1), (qnname_lz, 2 * lz), (qnname_spin_parity, 1, 2)) => 1
        ]
    elseif ((!conserve_spin_parity)&& conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,      0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz,     lz)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz,     lz)) => 1
            QN((qnname_nf, 2, -1), (qnname_lz, 2 * lz)) => 1
        ]        
    elseif (conserve_spin_parity && (!conserve_lz))
        return [
            QN((qnname_nf, 0, -1), (qnname_spin_parity, 0, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_spin_parity, 1, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_spin_parity, 0, 2)) => 1
            QN((qnname_nf, 2, -1), (qnname_spin_parity, 1, 2)) => 1
        ]
    else
        return [
            QN((qnname_spin_parity, 0, 2)) => 1
            QN((qnname_spin_parity, 1, 2)) => 1
            QN((qnname_spin_parity, 0, 2)) => 1
            QN((qnname_spin_parity, 1, 2)) => 1
        ]
    end
end

function generate_hmt_fuzzy_ising_4(nm :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16)
    os = OpIDSum{Float64}()
    @time "calc_int_ps" int_el = calc_int_ps(nm, ps_pot_u)
    for m1 in 1 : nm
        for m2 in 1 : nm
            for m3 in 1 : nm
                m4 = m1 + m2 - m3
                if (m4 <= 0 || m4 > nm) 
                    continue 
                end
                V = int_el[m1, m2, m3]

                push!(os, V * 2, opCdag(m1, ↑), opC(m4, ↑), opCdag(m2, ↓), opC(m3, ↓))
            end
        end
    end
    for m1 in 1 : nm
        push!(os, -fld_h, opCdag(m1, ↑), opC(m1, ↓))
        push!(os, -fld_h, opCdag(m1, ↓), opC(m1, ↓))
    end
    return os
end

nm = 24
sites = [siteind("Electron", lz = m, conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_spin_parity = false) for m :: Int in 1 : nm]
opCacheVec = create_op_cache_vec(sites)
os = generate_hmt_fuzzy_ising_4(nm)
hmt = MPO_new(os, sites, opCacheVec; basisOpCacheVec=opCacheVec, tol=4e-8)

nm = 100
sites = [siteind("Electron", lz = m, conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_spin_parity = false) for m :: Int in 1 : nm]
opCacheVec = create_op_cache_vec(sites)
println("nm = $nm")
@time "Constructing opsum" os = generate_hmt_fuzzy_ising_4(nm)
@time "constructing MPO" hmt = MPO_new(os, sites, opCacheVec; basisOpCacheVec=opCacheVec, tol=4e-8)

@show maxlinkdim(hmt)