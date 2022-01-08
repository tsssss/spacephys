;+
; on daily bases, survey plot includes both -A and -B
; vars: vsc, |E| high res, B1 and B2 availability (bars).
;-

; data buffer 2 day long.
; down sample to 1 sec.

utr0 = time_double(['2012-09-25','2016-12-31'])
utr0 = time_double(['2016-03-21','2016-12-31'])
;utr0 = time_double(['2012-09-30','2012-09-31'])
probes = ['a','b']
evals = [5,10,20,50] & evals = evals[sort(evals)]
rootdir = shomedir()+'/psbl_de_32hz/epacket'


dt0 = 86400d
spinrate = 11d  ; sec, 10.9378 sec.
dr0 = 1d/16     ; 16 Hz.
dr1 = 10        ; sec.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
spc = '    '


ut0s = smkarthm(utr0[0], utr0[1], dt0, 'dx')
nday = n_elements(ut0s)

neval = n_elements(evals)
ofns = strarr(neval)
for i = 0, neval-1 do begin
    ofns[i] = rootdir+'/rbsp_e_ge_'+sgnum2str(evals[i])+'.log'
    if file_test(ofns[i]) eq 0 then stouch, ofns[i]
endfor

rbsp_load_spice_kernels, trange = utr0


if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
foreach ut0, ut0s do begin

    ; current day.
    utr = ut0+[0,dt0]


    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'

    ; **** load vsvy and calc esvy (instead of load esvy).
    ; calc E uvw, then |E|.
    ; rbsp_emag0.
    
        tvar = pre0+'emag0'
        dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
        if size(dat,/type) ne 8 then begin
            store_data, tvar, /delete
            message, 'no V data ...', /continue
            continue
        endif
        
        uts = sfmepoch(dat.epoch, 'unix')
        nrec = n_elements(uts)
        vsvy = dat.vsvy

        ; the following is based on linear regression on v1,v2 and eu, v3,v4 and ev.
        eu = (vsvy[*,0]-vsvy[*,1])*10   ; V -> V/m -> mV/m.
        ev = (vsvy[*,2]-vsvy[*,3])*10
        vsvy = 0    ; free memory.
        ; remove offset, at spin rate.
        eu = eu-smooth(eu, spinrate/dr0, /edge_truncate, /nan)
        ev = ev-smooth(ev, spinrate/dr0, /edge_truncate, /nan)
        emag = sqrt(eu^2+ev^2)
        store_data, tvar, uts, emag, limits = $
            {ytitle:'|E|!C(mV/m)', labels:['|E| orig'], yrange:[0,200]}
        

    ; coefficient to remove E field near earth (< 3 Re).
    ; rbspx_dis.
    ; rbspx_emag.

        tvar = pre0+'dis'
        tuts = smkarthm(utr[0],utr[1],60,'dx')
        rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = tuts, /no_spice_load
        get_data, pre0+'state_pos_gse', tmp, posgse
        store_data, pre0+'state_*', /delete
        dis = snorm(posgse*re1)
        store_data, tvar, tuts, dis


        ; make a decaying window.
        mindis0 = 3 ; Re.
        get_data, pre0+'dis', tuts, dis
        dis = interpol(dis, tuts, uts)
        win0 = double(dis le mindis0)
        winwd = 9600/dr0
        winedge = exp(-findgen(winwd)/winwd*5)
        idx = where(win0-shift(win0,1) eq 1, cnt)   ; left edge.
        for i = 0, cnt-1 do begin
            j0 = idx[i]
            dj = (winwd-1)<(j0-0)
            win0[j0-dj:j0] = reverse(winedge[0:dj])
        endfor
        
        idx = where(win0-shift(win0,-1) eq 1, cnt)   ; right edge.
        for i = 0, cnt-1 do begin
            j0 = idx[i]
            dj = (winwd-1)<(nrec-1-j0)
            win0[j0:j0+dj] = (winedge[0:dj])
        endfor
        
        emagbg = dblarr(nrec)
        drec = (dis-mindis0)/(max(dis)-mindis0)>0
        drec = (drec^0.4*100+1)*spinrate/dr0*0.5
        for i = 0, nrec-1 do begin
            ti = i
            i1 = (ti-drec[i])>0
            i2 = (ti+drec[i])<nrec-1
            emagbg[ti] = min(emag[i1:i2])
        endfor
        emag = emag-emagbg
        store_data, pre0+'emag', uts, emag, limits = {ytitle:'|E|!C(mV/m)'}



    ; load the flags.
    ; flag for bad E field.
    ; rbspx_flag_[euv_diff,euv_asym,eclipse,sdt]
    ; rbspx_[euv_mag,euv_mag_diff_abs,euv_avg]
        
        
        ; compare magnitude of Eu and Ev, within 2 spins.
        ; if difference in magnitude is larger than 15 mV/m, then bad field.
        maxdiff = 15
        padt = 60 ; sec.
        drec = padt/dr0
        
        wd0 = 11
        euenv = scalcenv(eu, width = wd0/dr0)
        evenv = scalcenv(ev, width = wd0/dr0)
        euvdiff = abs(euenv-evenv)

        nrec = n_elements(uts)
        tdat = euvdiff ge maxdiff or finite(euvdiff,/nan)
        idx = where(tdat eq 1, cnt)
        for j = 0, cnt-1 do begin   ; expand by the amount of pad time.
            j0 = idx[j]-drec>0
            j1 = idx[j]+drec<nrec-1
            tdat[j0:j1] = 1
        endfor
        store_data, pre0+'flag_euv_diff', uts, tdat, limits = $
            {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
            labels:['1:Diff.Euv>'+sgnum2str(maxdiff)+'mV/m']}
            
            
        ; flag for bad E, when Eu and Ev are not symmetric.
        maxasym = 5 ; mV/m.
        padt = 60 ; sec.
        drec = padt/dr0
        
        tdat = [[eu],[ev]]
        nrec = n_elements(uts)
        
        for i = 0, 1 do $
            tdat[*,i] = smooth(tdat[*,i], 600/dr0, /nan, /edge_truncate)  ; 10 min average.
        tdat = (abs(tdat) ge maxasym)
        for i = 0, 1 do begin
            idx = where(tdat[*,i] eq 1, cnt)
            for j = 0, cnt-1 do begin
                j0 = idx[j]-drec>0
                j1 = idx[j]+drec<nrec-1
                tdat[j0:j1,i] = 1
            endfor
        endfor
        store_data, pre0+'flag_euv_asym', uts, tdat, limits = $
            {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
            colors:[6,4], $
            labels:['1:Asym.Euv>'+sgnum2str(maxasym)+'mV/m']}
            

        ; check eclipse time.
        tmp = sread_rbsp_eclipse_time(utr, probes = tprobe)
        flags = interpol(double(tmp.flags),tmp.uts,uts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
        store_data, pre0+'flag_eclipse', uts, flags, limits = $
            {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
            labels:['1:eclipse']}
            
        ; check sdt time.
        tmp = sread_rbsp_efw_sdt_time(utr[0], probes = tprobe)
        flags = interpol(double(tmp.flags),tmp.uts,uts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
        store_data, pre0+'flag_sdt', uts, flags, limits = $
            {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
            labels:['1:sdt']}
            
            
        flagvars = pre0+'flag_'+['euv_diff','euv_asym','eclipse','sdt']

        flags = bytarr(nrec)
        foreach tvar, flagvars do begin
            get_data, tvar, uts, tdat
            flags = flags or tdat
        endforeach
        store_data, pre0+'flag', uts, flags

        idx = where(flags eq 1, cnt)
        if cnt ne 0 then emag[idx] = !values.d_nan
        store_data, pre0+'emag', uts, emag
        
        
        
        ; down sample to 1 Hz, and only consider the envolope.
        get_data, pre0+'emag', tuts, emag
        uts = smkarthm(utr[0], utr[1], dr1, 'dx')
        nrec = n_elements(uts)-1
        eenv = dblarr(nrec)
        for i = 0, nrec-1 do eenv[i] = max(emag[where(tuts ge uts[i] and tuts lt uts[i+1])],/nan)
        store_data, pre0+'eenv', uts[0:nrec-1], eenv
        

        for i = 0, neval-1 do begin
            idx = where(eenv ge evals[i], cnt)
            if cnt eq 0 then break      ; evals are sorted, no need to check higher values.

            ; add 0 on both sides treat the case when large E event is
            ; at the beginning or end of a day.
            flags = [0,eenv ge evals[i],0]  ; nrec+2 elements.
            flags = flags[1:nrec+1]-flags[0:nrec]   ; nrec+1 elements, same as uts.
            

            ; should be equal in size.
            t1 = uts[where(flags eq 1, cnt1)]
            t2 = uts[where(flags eq -1, cnt2)]
            if cnt1 ne cnt2 then message, 'something wrong here ...'

            lines = strarr(cnt1)
            for j = 0, cnt1-1 do $
                lines[j] = time_string(t1[j], tformat='YYYY_MMDD'+spc+'hh:mm:ss')+spc+$
                    time_string(t2[j], tformat='hh:mm:ss')+spc+$
                    string(max(eenv[where(uts ge t1[j] and uts le t2[j])],/nan),format='(I3)')+spc+$
                    strupcase(tprobe)+spc
            openw, lun, ofns[i], /get_lun, /append
            printf, lun, lines
            free_lun, lun
        endfor
    endforeach
endforeach

rbsp_load_spice_kernels, trange = utr0, /unload
    
end
