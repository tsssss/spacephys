;+
; Quality flag for E field, 1 for bad data.
; Exclude E field around perigee, L<3.
; re: set to indicate the posvar is in Re.
;-
pro rbsp_efw_perigee_flag, utr, probes = probes, mindis = mindis, pad = dut, $
    posvar = posvar, re = re

    if n_elements(utr) eq 0 then message, 'no time ...'
    
    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    
    if n_elements(mindis) eq 0 then mindis = 3
    if n_elements(posvar) eq 0 then posvar = 'pos_gse'
    
    re = (n_elements(re) eq 0)? 6378d: 1d
    
    ; use s/c distance to check. rbsp[ab]_dis.
    for i = 0, nprobe-1 do begin
        tprobe = probes[i]
        pre0 = 'rbsp'+tprobe+'_'
        tvar = pre0+'dis'
        load = stplot_load_flag(tvar, trange = utr, mode = 'bigger')
        if load then begin
            ; if pos exists, calc dis from it.
            tmp = stplot_load_flag(pre0+posvar, trange = utr, mode = 'bigger')
            if tmp eq 0 then begin
                get_data, pre0+posvar, uts, tmp
                store_data, tvar, uts, sqrt(total(tmp^2,2))/re
            endif else begin
                ; load position and calc dis.
                tmp = sread_rbsp_efw_l3(utr, probes = tprobe, $
                    vars = ['epoch','pos_gse'])
                uts = sfmepoch(tmp.epoch, 'unix', /epoch16)
                dis = sqrt(total(tmp.pos_gse^2,2))/re
                store_data, tvar, uts, dis
            endelse
        endif
        get_data, tvar, uts, dis
        nrec = n_elements(uts)
        flag = bytarr(nrec)
        idx = where(dis le mindis, cnt)
        if cnt ne 0 then flag[idx] = 1
        tvar = pre0+'perigee_flag'
        store_data, tvar, uts, flag
        
        options, tvar, 'yrange', [-0.5,1.5]
        options, tvar, 'ytitle', 'Quality Flag!C1 for L>'+sgnum2str(mindis)
        
        ; add pad time.
        if n_elements(dut) ne 0 then $
            stplot_enlarge_flag, tvar, dut = dut, flag0 = 1
    endfor
    
    ; use B magnitude to check. rbsp[ab]_b.

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

rbsp_efw_perigee_flag, utr, probe = tprobe


device, decomposed = 0
loadct2, 43

vars = pre0+['desvy_mgse','perigee_flag']
tplot, vars, trange = utr
end