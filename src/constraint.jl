function obj(g, pg)
    return g.cost1 * pg^2 + g.cost2 * pg + g.cost3 
end

function c1_polar(va)
    return va
end

function c2_polar(b, p_f, vm_f, vm_t, va_f, va_t)
    return p_f - b.c5 * vm_f^2 -
    b.c3 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c4 * (vm_f * vm_t * sin(va_f - va_t))
end

function c3_polar(b, q_f, vm_f, vm_t, va_f, va_t)
    return q_f +
    b.c6 * vm_f^2 +
    b.c4 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c3 * (vm_f * vm_t * sin(va_f - va_t))
end

function c4_polar(b, p_t, vm_f, vm_t, va_f, va_t)
    return p_t - b.c7 * vm_t^2 -
    b.c1 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c2 * (vm_t * vm_f * sin(va_t - va_f))
end

function c5_polar(b, q_t, vm_f, vm_t, va_f, va_t)
    return q_t +
    b.c8 * vm_t^2 +
    b.c2 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c1 * (vm_t * vm_f * sin(va_t - va_f))
end

function c6_polar(b, va_f, va_t)
    return va_f - va_t
end

function c7_polar(b, vm)
    return b.pd + b.gs * vm^2
end

function c8_polar(b, vm)
    return b.qd - b.bs * vm^2
end

#no coordinates specified
function c9_10(b, p,q)
    return p^2 + q^2 - b.rate_a_sq
end

#rectangular constraints
function c1_rect(vr, vim)
    return atan(vim/vr)
end

function c2_rect(b, p_f, vr_f, vr_t, vim_f, vim_t)
    return p_f - b.c5 * (vr_f^2+vim_f^2) -
    b.c3 * (vr_f*vr_t + vim_f*vim_t) -
    b.c4 * (vim_f*vr_t - vr_f*vim_t) 
end

function c3_rect(b, q_f, vr_f, vr_t, vim_f, vim_t)
    return q_f +
    b.c6 * (vr_f^2+vim_f^2) +
    b.c4 * (vr_f*vr_t + vim_f*vim_t) -
    b.c3 * (vim_f*vr_t - vr_f*vim_t)
end

function c4_rect(b, p_t, vr_f, vr_t, vim_f, vim_t)
    return  p_t - b.c7 * (vr_t^2+vim_t^2) -
    b.c1 * (vr_f*vr_t + vim_f*vim_t) -
    b.c2 * (vim_t*vr_f - vr_t*vim_f)
end

function c5_rect(b, q_t, vr_f, vr_t, vim_f, vim_t)
    return q_t +
        b.c8 * (vr_t^2+vim_t^2) +
        b.c2 * (vr_f*vr_t + vim_f*vim_t) -
        b.c1 * (vim_t*vr_f - vr_t*vim_f)
end

function c6_rect(b, vr_f, vr_t, vim_f, vim_t)
    return atan(vim_f/vr_f) - atan(vim_t/vr_t) 
end

function c7_rect(b, vr, vim)
    return b.pd + b.gs * (vr^2 + vim^2)
end

function c8_rect(b, vr, vim)
    return b.qd - b.bs * (vr^2 + vim^2)
end

function c11_rect(vr, vim)
    return vr^2 + vim^2
end

#only for mp
function c_12(pg0, pg1)
    return pg0 - pg1
end

#storage constraints
function c_13(s, pst, pstd, pstc, I2)
    return pst + pstd - pstc - s.Pexts - s.Zr*I2
end

function c13_smooth(s, pst, pstd, I2)
    return pst + pstd - s.Pexts - s.Zr*I2
end

function c_14(s, qst, qint, I2)
    return qst - qint - s.Qexts - s.Zim*I2
end

function c15_polar(pst, qst, vm, I2)
    return pst^2 + qst^2 - (vm^2)*I2
end

function c15_rect(pst, qst, vr, vim, I2)
    return pst^2 + qst^2 - (vr^2 + vim^2)*I2
end

function c16_17(s, E0, E1, pstc, pstd)
    return E0 - E1 - (s.etac*pstc - pstd/s.etad)
end

function c16_17_smooth(s, E0, E1, discharge_func::Function, pstd)
    return E0 - E1 + discharge_func(pstd, s.Srating)
end

function c_18(s, pst, qst)
    return pst^2 + qst^2 - s.Srating^2
end

function c_19(pstd, pstc)
    return pstd - pstc
end

function c19_smooth(pstd)
    return pstd
end

function c_20(pstd, pstc)
    return pstd*pstc
end