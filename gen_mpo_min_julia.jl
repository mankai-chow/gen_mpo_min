using LinearAlgebra
using WignerSymbols
using ITensors

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
    os = OpSum()
    int_el = calc_int_ps(nm, ps_pot_u)
    for m1 in 1 : nm
        for m2 in 1 : nm
            for m3 in 1 : nm
                m4 = m1 + m2 - m3
                if (m4 <= 0 || m4 > nm) 
                    continue 
                end
                V = int_el[m1, m2, m3]
                os += V * 2, "Cdagup", m1, "Cup", m4, "Cdagdn", m2, "Cdn", m3
            end
        end
    end
    for m1 in 1 : nm
        os += -2 * fld_h, "Sx", m1
    end
    return os
end

nm = 24
sites = [siteind("Electron", lz = m, conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_spin_parity = false) for m :: Int in 1 : nm]

os = generate_hmt_fuzzy_ising_4(nm)
t0 = time()
hmt = MPO(os, sites ; cutoff = 1E-13)
t1 = time() - t0
@show t1
t0 = time()
hmt1 = MPO(os, sites ; cutoff = 1E-13)
t1 = time() - t0
@show t1
@show maxlinkdim(hmt)