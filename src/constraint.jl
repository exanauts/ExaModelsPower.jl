function gen_cost(g, pg)
    return g.cost1 * pg^2 + g.cost2 * pg + g.cost3 
end

function c_ref_angle_polar(va)
    return va
end

function c_to_active_power_flow_polar(b, p_f, vm_f, vm_t, va_f, va_t)
    return p_f - b.c5 * vm_f^2 -
    b.c3 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c4 * (vm_f * vm_t * sin(va_f - va_t))
end

function c_to_reactive_power_flow_polar(b, q_f, vm_f, vm_t, va_f, va_t)
    return q_f +
    b.c6 * vm_f^2 +
    b.c4 * (vm_f * vm_t * cos(va_f - va_t)) -
    b.c3 * (vm_f * vm_t * sin(va_f - va_t))
end

function c_from_active_power_flow_polar(b, p_t, vm_f, vm_t, va_f, va_t)
    return p_t - b.c7 * vm_t^2 -
    b.c1 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c2 * (vm_t * vm_f * sin(va_t - va_f))
end

function c_from_reactive_power_flow_polar(b, q_t, vm_f, vm_t, va_f, va_t)
    return q_t +
    b.c8 * vm_t^2 +
    b.c2 * (vm_t * vm_f * cos(va_t - va_f)) -
    b.c1 * (vm_t * vm_f * sin(va_t - va_f))
end

function c_phase_angle_diff_polar(b, va_f, va_t)
    return va_f - va_t
end

function c_active_power_balance_demand_polar(b, vm)
    return b.pd + b.gs * vm^2
end

function c_reactive_power_balance_demand_polar(b, vm)
    return b.qd - b.bs * vm^2
end

#no coordinates specified
function c_thermal_limit(b, p,q)
    return p^2 + q^2 - b.rate_a_sq
end

#only for mp
function c_ramp(pg0, pg1)
    return pg0 - pg1
end

#rectangular constraints
function c_ref_angle_rect(vr, vim)
    return atan(vim/vr)
end

function c_to_active_power_flow_rect(b, p_f, vr_f, vr_t, vim_f, vim_t)
    return p_f - b.c5 * (vr_f^2+vim_f^2) -
    b.c3 * (vr_f*vr_t + vim_f*vim_t) -
    b.c4 * (vim_f*vr_t - vr_f*vim_t) 
end

function c_to_reactive_power_flow_rect(b, q_f, vr_f, vr_t, vim_f, vim_t)
    return q_f +
    b.c6 * (vr_f^2+vim_f^2) +
    b.c4 * (vr_f*vr_t + vim_f*vim_t) -
    b.c3 * (vim_f*vr_t - vr_f*vim_t)
end

function c_from_active_power_flow_rect(b, p_t, vr_f, vr_t, vim_f, vim_t)
    return  p_t - b.c7 * (vr_t^2+vim_t^2) -
    b.c1 * (vr_f*vr_t + vim_f*vim_t) -
    b.c2 * (vim_t*vr_f - vr_t*vim_f)
end

function c_from_reactive_power_flow_rect(b, q_t, vr_f, vr_t, vim_f, vim_t)
    return q_t +
        b.c8 * (vr_t^2+vim_t^2) +
        b.c2 * (vr_f*vr_t + vim_f*vim_t) -
        b.c1 * (vim_t*vr_f - vr_t*vim_f)
end

function c_phase_angle_diff_rect(b, vr_f, vr_t, vim_f, vim_t)
    return atan(vim_f/vr_f) - atan(vim_t/vr_t) 
end

function c_active_power_balance_demand_rect(b, vr, vim)
    return b.pd + b.gs * (vr^2 + vim^2)
end

function c_reactive_power_balance_demand_rect(b, vr, vim)
    return b.qd - b.bs * (vr^2 + vim^2)
end

function c_voltage_magnitude_rect(vr, vim)
    return vr^2 + vim^2
end
