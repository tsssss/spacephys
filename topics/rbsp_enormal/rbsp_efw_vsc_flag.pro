;+
;-
pro rbsp_efw_vsc_flag, utr, probes = probes, minvsc = minvsc, pad = dut

    if n_elements(utr) eq 0 then message, 'no time ...'
    
    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    
    if n_elements(minvsc) eq 0 then minvsc = -30    ; 30 V.
    
    ; use s/c potential to check. rbsp[ab]_vsc.
    for i = 0, nprobe-1 do begin
        tprobe = probes[i]
        pre0 = 'rbsp'+tprobe+'_'
        tvar = pre0+'vsc'
        load = stplot_load_flag(tvar, trange = utr, mode = 'bigger')
        if load then begin
            ; if vsvy exists, calc vsc from it.
            tmp = stplot_load_flag(pre0+'vsvy', trange = utr, mode = 'bigger')
            if tmp eq 0 then begin
                get_data, pre0+'vsvy', uts, tmp
                store_data, tvar, uts, 0.5*(tmp[*,2]+tmp[*,3])
            endif else begin
                tmp = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
                uts = sfmepoch(tmp.epoch, 'unix', /epoch16)
                tmp = tmp.vsvy
                store_data, tvar, uts, 0.5*(tmp[*,2]+tmp[*,3])
            endelse
        endif
        get_data, tvar, uts, vsc
        nrec = n_elements(uts)
        flag = bytarr(nrec)
        idx = where(vsc le minvsc, cnt)
        if cnt ne 0 then flag[idx] = 1
        tvar = pre0+'vsc_flag'
        store_data, tvar, uts, flag
        
        options, tvar, 'yrange', [-0.5,1.5]
        options, tvar, 'ytitle', 'Quality Flag!C1 for Vsc>'+sgnum2str(minvsc)
        
        ; add pad time.
        if n_elements(dut) ne 0 then $
            stplot_enlarge_flag, tvar, dut = dut, flag0 = 1
    endfor

end

tprobe = 'b'
utr = time_string(['2013-05-01','2013-05-02'])

rgb = [6,4,2]
pre0 = 'rbsp'+tprobe+'_'
tplot_options, 'ystyle', 1

efw = sread_rbsp_efw_l3(utr, probes = tprobe)
store_data, pre0+'desvy_mgse', sfmepoch(efw.epoch, 'unix', /epoch16), $
    efw.efield_inertial_frame_mgse, $
    limits = {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
efw = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'euvw')
store_data, pre0+'deuvw', sfmepoch(efw.epoch, 'unix', /epoch16), $
    efw.efield_uvw, $
    limits = {ytitle:'dE!C(mV/m)', colors:rgb, labels:['u','v','w']}
efw = 0

padtime = 5*60  ; sec.
rbsp_efw_perigee_flag, utr, probe = tprobe, pad = padtime
rbsp_efw_vsc_flag, utr, probe = tprobe, pad = padtime
; use euvw? can be +/-200 in magnitude (in 2013-05-01 -B).
; u and v amplitude should be comparable.


device, decomposed = 0
loadct2, 43

tplot_options, 'constant', 0

vars = pre0+['desvy_mgse','dedespun_mgse','deuvw','perigee_flag','vsc_flag']
tplot, vars, trange = utr
stop

flags = pre0+['perigee_flag','vsc_flag']
vars = pre0+['desvy_mgse','dedespun_mgse']
for i = 0, n_elements(vars)-1 do begin
    get_data, vars[i], uts, dat
    for j = 0, n_elements(flags)-1 do begin
        get_data, flags[j], tmp, flag
        flag = interpol(flag, tmp, uts)
        idx = where(flag eq 1, cnt)
        if cnt ne 0 then dat[idx,*] = !values.d_nan
    endfor
    store_data, vars[i], uts, dat
endfor

vars = pre0+['desvy_mgse','dedespun_mgse','deuvw','perigee_flag','vsc_flag']
tplot, vars, trange = utr

end