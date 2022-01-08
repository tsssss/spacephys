
pro rbsp_efw_euvw_flag, utr, probes = probes, pad = dut, maxdamp = maxdamp, evar = evar

    if n_elements(utr) eq 0 then message, 'no time ...'
    if n_elements(maxdamp) eq 0 then maxdamp = 50   ; 50 mV/m.
    if n_elements(evar) eq 0 then evar = 'euvw'
    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    
    ; use euvw to check, rbsp[ab]_euvw.
    for i = 0, nprobe-1 do begin
        tprobe = probes[i]
        pre0 = 'rbsp'+tprobe+'_'
        tvar = pre0+evar
        load = stplot_load_flag(tvar, trange = utr, mode = 'bigger')
        if load then begin
            tmp = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'euvw')
            uts = sfmepoch(tmp.epoch, 'unix', /epoch16)
            store_data, tvar, uts, tmp.efield_uvw
        endif
        get_data, tvar, uts, euvw
        nrec = n_elements(uts)
        amps = fltarr(nrec,2)
        dr = sdatarate(uts)
        for j = 0,2-1 do begin
            tfs = smooth(euvw[*,j], 3/dr, /nan)
            dfs = [0,tfs[1:nrec-1]-tfs[0:nrec-2]]
            dfs = smooth(dfs, 3/dr, /nan)
            idx = where(dfs[1:nrec-1]*dfs[0:nrec-2] le 0, cnt)
            amps[*,j] = interpol(abs(tfs[idx]),uts[idx],uts)
        endfor
              
        damp = abs(amps[*,0]-amps[*,1])
        flag = bytarr(nrec)
        idx = where(damp ge maxdamp, cnt)
        if cnt ne 0 then flag[idx] = 1
        tvar = pre0+'euvw_flag'
        store_data, tvar, uts, flag
        
        options, tvar, 'yrange', [-0.5,1.5]
        options, tvar, 'ytitle', 'Quality Flag!C1 for |Eu-Ev|<'+sgnum2str(maxdamp)
        
        ; add pad time.
        if n_elements(dut) ne 0 then $
            stplot_enlarge_flag, tvar, dut = dut, flag0 = 1
    endfor

end

tprobe = 'a'
utr = time_string(['2013-05-01','2013-05-02'])

;tprobe = 'b'
;utr = time_string(['2013-06-01','2013-06-02'])

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
efw = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'esvy')
store_data, pre0+'dedespun_mgse', sfmepoch(efw.epoch, 'unix', /epoch16), $
    efw.efield_mgse, $
    limits = {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
efw = 0

padtime = 5*60  ; sec.
rbsp_efw_perigee_flag, utr, probe = tprobe, pad = padtime
rbsp_efw_vsc_flag, utr, probe = tprobe, pad = padtime
rbsp_efw_euvw_flag, utr, probe = tprobe, pad = padtime


device, decomposed = 0
loadct2, 43

tplot_options, 'constant', 0
tplot_options, 'labflag', -1

vars = pre0+['desvy_mgse','dedespun_mgse','deuvw', $
    'perigee_flag','vsc_flag','euvw_flag']
tplot, vars, trange = utr
stop

flags = pre0+['perigee_flag','vsc_flag','euvw_flag']
vars = pre0+['desvy_mgse','dedespun_mgse','deuvw']
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

vars = pre0+['desvy_mgse','dedespun_mgse','deuvw','perigee_flag','vsc_flag','euvw_flag']
tplot, vars, trange = utr
end