;+
; Calculate the drift period in sec for a given particle at certain energy, pitch angle, and location,
; using a relativistic Lorentz equation solver.
; To substitute calc_drift_velocity b/c the drift period is the standard quantity.
;
; rgsm, the position of measurements in GSM in Re.
; time, the UT second.
; energy, the energy in keV.
; pitch_angle, the pitch angle in deg.
; species, a string of 'e', 'h+', 'o+', etc.
; model, a string sets 'igrf', or 'dipole'. 'igrf' is the default option.
; duration, the duration of the simulation. Default is 2 sec.
; bounce_period, in sec, output.
; lshell, in #, output.
;-


function calc_drift_period, rgsm=rgsm0, time=ut0, $
    energy=energy0, pitch_angle=pitch_angle, species=species, $
    model=model, duration=duration, $
    bounce_period=bounce_period, lshell=lshell


;---Constants.
    c = 299792.458d         ; in km/s.
    me = 9.10938356d-31     ; kg.
    mp = 1.6726219d-27      ; kg.
    q0 = 1.60217662d-19     ; C.

    deg = 180d/!dpi
    rad = !dpi/180d

    re = 6378d      ; in km.

;---Settings.
    if n_elements(species) eq 0 then message, 'no species ...'
    if n_elements(rgsm0) ne 3 then message, 'no position ...'
    if n_elements(ut0) eq 0 then message, 'no time ...'
    if n_elements(energy0) eq 0 then message, 'no energy ...'
    if n_elements(pitch_angle) eq 0 then pitch_angle = 90.
    pitch_angle0 = pitch_angle*rad
    if n_elements(model) eq 0 then model = 'igrf'

    step_per_gyration = 32d ; # of steps per gyro-motion.

    ; dispatch from species to mass and charge.
    case strlowcase(species) of
        'e': begin
            mass = me
            charge = -q0
            duration_factor = 2
            end
        'p': begin
            mass = mp
            charge = q0
            duration_factor = 40
            end
        'h': begin
            mass = mp
            charge = q0
            duration_factor = 40
            end        
        'h+': begin
            mass = mp
            charge = q0
            duration_factor = 40
            end
        'o+': begin
            mass = mp*16
            charge = q0
            duration_factor = 40
            end
        'he+': begin
            mass = mp*4
            charge = q0
            duration_factor = 40
            end
        'o++': begin
            mass = mp*16
            charge = q0*2
            duration_factor = 40
            end
        else: message, 'do not support '+species+' yet ...'
    endcase
    if n_elements(duration) eq 0 then duration = duration_factor*(400d/energy0) ; sec, inverse scaled with energy.


;---Derived constants.
    re_1 = 1d/re
    c_1 = 1d/c
    twopi = 2*!dpi
    q_mxtwopi = charge/(mass*twopi)     ; C/kg.
    e0 = mass*c^2*1e3/q0                ; the rest mass energy in keV, for c in km/s.

    ; get the initial B field using model.
    ps = geopack_recalc(ut0)
    case model of
        'igrf': geopack_igrf_gsm, rgsm0[0], rgsm0[1], rgsm0[2], bx,by,bz
        'dipole': geopack_dip, rgsm0[0], rgsm0[1], rgsm0[2], bx,by,bz
    endcase
    bgsm0 = [bx,by,bz]  ; in nT.
    bmag0 = snorm(bgsm0)

    gamma0 = energy0/e0+1               ; gamma = 1/sqrt(1-(v/c)^2).
    v0 = sr_r2v(gamma0)                 ; velocity in km/s.
    tauc0 = abs((twopi*mass*gamma0)/(charge*bmag0*1d-9))    ; gyro-period in sec.
    step_dt = tauc0/step_per_gyration   ; time step used in simulation.
    nstep = floor(duration/step_dt)

    ; calculate the velocity in GSM, choose a convenient phase.
    bhat0 = sunitvec(bgsm0)
    what0 = sunitvec(vec_cross(rgsm0,bgsm0))
    ohat0 = vec_cross(bhat0,what0)
    mgsm2fac = [[bhat0],[what0],[ohat0]]
    vfac0 = v0*[cos(pitch_angle0),sin(pitch_angle0),0]
    ; xhat0 = reform(mgsm2fac[0,*])
    ; yhat0 = reform(mgsm2fac[1,*])
    ; zhat0 = reform(mgsm2fac[2,*])
    vgsm0 = mgsm2fac # vfac0


;---Run the relativistic Lorentz equation solver.
    times = ut0+step_dt*findgen(nstep)
    rgsms = dblarr(nstep,3)
    vgsms = dblarr(nstep,3)
    rgsms[0,*] = rgsm0
    vgsms[0,*] = vgsm0

    sr_lorentz_eq_solver, times=times, rgsms, vgsms, bmags=bmags, $
        charge=charge, mass=mass


;---Use |B| to find bounce period.
    ; find equator.
    beqs = smooth(bmags[*],step_per_gyration,/edge_truncate)
    beqs = beqs-mean(beqs)
    nodes = find_node(beqs)
    nnode = n_elements(nodes)
    node_times = times[nodes]
    periods = node_times[1:nnode-1]-node_times[0:nnode-2]

    bounce_period = mean(periods)*2
    lshell = mean(snorm(rgsms[nodes,0:1]))


;---Use SM coord to fit azimuthal velocity.
    rsms = cotran(rgsms, times, 'gsm2sm')
    ; remove modulations due to gyro-motion.
    xsms = smooth(rsms[*,0],step_per_gyration,/edge_truncate)
    ysms = smooth(rsms[*,1],step_per_gyration,/edge_truncate)
    ; remove modulations due to bounce-motion.
    step_per_bounce = 2*bounce_period/step_dt
    xsms = smooth(xsms,step_per_bounce,/edge_truncate)
    ysms = smooth(ysms,step_per_bounce,/edge_truncate)
    xsms = xsms[step_per_bounce:nstep-1-step_per_bounce]
    ysms = ysms[step_per_bounce:nstep-1-step_per_bounce]
    txs = times[step_per_bounce:nstep-1-step_per_bounce]-ut0
    tys = atan(ysms, xsms)
    dys = tys[1:-1]-tys[0:-1]
    index = where(abs(dys) ge 1.5*!dpi, count)
    for ii=0, count-1 do begin
        tmp = total(tys[index[ii]:index[ii]+1]*[-1,1])
        if tmp ge 0 then begin
            tys[index[ii]+1:*] -= 2*!dpi
        endif else begin
            tys[index[ii]+1:*] += 2*!dpi
        endelse
    endfor
    res = linfit(txs, tys)


    ;v_drift = res[1]*deg    ; in deg/sec.
    drift_period = 2*!dpi/res[1]    ; in sec.
    return, drift_period

end


ut0 = time_double('2014-08-28/10:18')
rgsm0 = [-1d,6.2,1.2]           ; in Re.
energys = [40,75,150,275,475]
energy0 = 300d              ; in keV.
pitch_angle0 = 75d          ; in deg.
species = 'e'
model = 'igrf'

foreach energy0, energys do begin
    drift_period = calc_drift_period(rgsm=rgsm0, time=ut0, $
        energy=energy0, pitch_angle=pitch_angle0, species=species, $
        model=model, bounce_period=bounce_period)
    print, energy0, drift_period
endforeach

end
