function constraint_polar_branch(b, p_f,p_t, q_f,q_t, vm_f, vm_t, va_f, va_t, i)


    c2 = p_f - b.c5 * vm_f^2 -
    b.c3 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c4 * (vm_f * vm_t * sin(va_f - va_t))

    c3 = q_f +
    b.c6 * vm_f^2 +
    b.c4 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c3 * (vm_f * vm_t * sin(va_f - va_t))

    c4 = p_t - b.c7 * vm_t^2 -
    b.c1 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c2 * (vm_t * vm_f * sin(va_t - va_f))

    c5 = q_t +
    b.c8 * vm_t^2 +
    b.c2 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c1 * (vm_t * vm_f * sin(va_t - va_f))
    
    c6 = va_f - va_t

    c9 = p_f^2 + q_f^2 - b.rate_a_sq

    c10 = p_t^2 + q_t^2 - b.rate_a_sq

    cs = [c2, c3, c4, c5, c6, c9, c10]
    return cs[i]
end

function constraint_polar_bus(b, vm, i)
    c7 = b.pd + b.gs * vm^2

    c8 = b.qd - b.bs * vm^2

    cs = [c7,c8]
    return cs[i]
end

function constraint_rect_branch(b, p_f,p_t, q_f,q_t, vr_f, vr_t, vim_f, vim_t, i)
    c2 = p_f - b.c5 * (vr_f^2+vim_f^2) -
            b.c3 * (vr_f*vr_t + vim_f*vim_t) -
            b.c4 * (vim_f*vr_t - vr_f*vim_t) 

    c3 = q_f +
        b.c6 * (vr_f^2+vim_f^2) +
        b.c4 * (vr_f*vr_t + vim_f*vim_t) -
        b.c3 * (vim_f*vr_t - vr_f*vim_t)
    

    c4 = p_t - b.c7 * (vr_t^2+vim_t^2) -
        b.c1 * (vr_f*vr_t + vim_f*vim_t) -
        b.c2 * (vim_t*vr_f - vr_t*vim_f)
    

    c5 = q_t +
        b.c8 * (vr_t^2+vim_t^2) +
        b.c2 * (vr_f*vr_t + vim_f*vim_t) -
        b.c1 * (vim_t*vr_f - vr_t*vim_f)

    c6 = atan(vim_f/vr_f) - atan(vim_t/vr_t) 
        
    c9 = p_f^2 + q_f^2 - b.rate_a_sq

    c10 = p_t^2 + q_t^2 - b.rate_a_sq
    
    cs = [c2, c3, c4, c5, c6, c9, c10]

    return cs[i]
end

function constraint_rect_bus(b, vr, vim, i)

    c7 = b.pd + b.gs * (vr^2 + vim^2)

    c8 = b.qd - b.bs * (vr^2 + vim^2)

    c11 = vr^2 + vim^2

    cs = [c7, c8, c11]

    return cs[i]
end