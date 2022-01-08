;+
; Solve the relativistic Lorentz equation.
; 
; times, an array in [n] of UT sec. times[0] is the initial time,
;   times[1]-times[0] is the time step.
;   
;   If times is omitted, then use rgsms/3 to set nstep. And set ut0 and dt.
;   
; rgsms, an array in [n,3] of position in GSM in Re. The first element is the 
;   initial position, the others are the outputs from the simulation.
; vgsms, an array in [n,3] of velocity in GSM in km/s. The first element is
;   the initial velocity, the others are the outpus from the simulation.
; bmags, an array in [n] of |B| in nT, optional. Set to take |B| as output.
; model, a string sets 'igrf', or 'dipole'. 'igrf' is the default option.
; 
; charge, in Coulumb, required.
; mass, in kg, required.
;-


pro sr_lorentz_eq_solver, rgsms, vgsms, bmags=bmags, times=times, $
    ut0=ut0, dt=step_dt, model=model, charge=charge, mass=mass
    
    re = 6378d      ; in km.
    re_1 = 1d/re
    c = 299792.458d         ; in km/s.
    c_1 = 1d/c



    ; re-configure model.
    nstep = n_elements(times)
    if nstep ne 0 then begin
        ut0 = times[0]
        step_dt = times[1]-times[0]
    endif else nstep = n_elements(rgsms)/3
    if n_elements(ut0) ne 0 then ps = geopack_recalc(ut0)
    if n_elements(step_dt) eq 0 then step_dt = 1d-5 ; sec.
    if n_elements(bmags) eq 0 then bmags = fltarr(nstep)
    if n_elements(model) eq 0 then model = 'igrf'
    if n_elements(charge) eq 0 then message, 'no charge ...'
    if n_elements(mass) eq 0 then message, 'no mass ...'
    
    qdt_2m = (charge*step_dt)/(2*mass)*1d-9  ; q*dt/(2*m) for B in nT and v in km/s.

    case model of
        'igrf': geopack_igrf_gsm, rgsms[0,0], rgsms[0,1], rgsms[0,2], bx,by,bz
        'dipole': geopack_dip, rgsms[0,0], rgsms[0,1], rgsms[0,2], bx,by,bz
    endcase
    bmags[0] = snorm([bx,by,bz])
    
    
    for i=1, nstep-1 do begin
        rgsm0 = rgsms[i-1,*]
        vgsm0 = vgsms[i-1,*]
        
    ;---update position, E and B field.
        rgsm1 = rgsm0 + vgsm0*step_dt*re_1  ; in Re.
        
        case model of
            'igrf': geopack_igrf_gsm, rgsm1[0], rgsm1[1], rgsm1[2], bx,by,bz
            'dipole': geopack_dip, rgsm1[0], rgsm1[1], rgsm1[2], bx,by,bz
        endcase
        bgsm1 = [bx,by,bz]  ; in nT.
        
        egsm1 = [0d,0,0]    ; in muV/m.
        
    ;---update velocity.
        gamma = sr_v2r(vgsm0)
        ugsm = gamma*vgsm0
        ugsm1 = ugsm + qdt_2m*(egsm1+vec_cross(vgsm0, bgsm1))  ; in km/s, for E in muV/m.
        uprime = ugsm1 + qdt_2m*egsm1   ; for E in muV/m.
        tau = qdt_2m*bgsm1  ; in km/s.
        ustar = vec_dot(uprime,tau)*c_1 ; in km/s.
        tmp = uprime*c_1
        gammap = sqrt(1d +vec_dot(tmp,tmp))
        tau2 = vec_dot(tau,tau)
        sigma = gammap^2 - tau2
        gamma1 = sqrt(0.5*(sigma+sqrt(sigma^2+4*(tau2+ustar^2))))
        t = tau/gamma1
        s = 1d/(1d + vec_dot(t,t))
        ugsm1 = s*(uprime + vec_dot(uprime,t)*t + vec_cross(uprime,t))
        vgsm1 = ugsm1/gamma1
        
        rgsms[i,*] = rgsm1
        vgsms[i,*] = vgsm1
        bmags[i] = snorm(bgsm1)
    endfor

end