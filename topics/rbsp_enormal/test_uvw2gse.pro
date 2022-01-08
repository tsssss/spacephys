;+
; test uvw2gse, use cspice_pxform, mtoq, qslerp, qtom
; conclusion: use interpolation through quaternion is reliable and fast.
;-

utr = time_double(['2012-11-14/04:00','2012-11-14/04:10'])
tprobe = 'b'
pre0 = 'rbsp'+tprobe+'_'
spinrate = 10.937774d
dr0 = 1d/16

tplot_options, 'labflag', -1
tplot_options, 'ystyle', 1



; **** load Vsvy.
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
    
    uts = sfmepoch(dat.epoch, 'unix')
    nrec = n_elements(uts)

; **** calc E UVW.
    tvar = pre0+'euvw'    
    vsvy = dat.vsvy
    eu = (vsvy[*,0]-vsvy[*,1])*10   ; V -> V/m -> mV/m.
    ev = (vsvy[*,2]-vsvy[*,3])*10
    ew = (vsvy[*,4]-vsvy[*,5])*2.4
    ew[*] = 0
    ; remove offset.
    nsec = 48   ; 30 min. tested to be very comparable to Euvw.
    eu = eu-scalcbg(eu, nsection = nsec)
    ev = ev-scalcbg(ev, nsection = nsec)
    
    store_data, tvar, uts, [[eu],[ev],[ew]], limits = $
        {ytitle:'E UVW!C(mV/m)', labels:['Eu','Ev','Ew'], colors:[6,4,2], yrange:[-200,200]}


; **** convert E UVW to GSE.
    tvar = pre0+'egse'
    rbsp_load_spice_kernels, trange = utr, probe = tprobe
    scid = strupcase(pre0+'science')


    ; method 1: use full time resolution.
        tmp = time_string(uts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, et0     ; convert ISO time string to SPICE ET
        ets = et0+uts-uts[0]
	    
        ; test if the last time records agree.
        tmp = time_string(uts[nrec-1],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, et1
        print, et1-ets[nrec-1]  ; diff is 5e-7 sec.
    
        cspice_pxform, scid, 'GSE', ets, pxform1    ; in [3,3,9600].


    ; method 2: use 1 sec time resolution.
        dt0 = 1
        tutr = utr-(utr mod dt0) & if (utr[1] mod dt0) ne 0 then tutr[1]+= dt0
        tuts = smkarthm(tutr[0], tutr[1], dt0, 'dx')
        tnrec = n_elements(tuts)

        tmp = time_string(tuts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+tuts-tuts[0]
    
        ; test if the last time records agree.
        tmp = time_string(tuts[tnrec-1],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet1
        print, tet1-tets[tnrec-1]  ; diff is 6e-8 sec.
        
        cspice_pxform, scid, 'GSE', tets, tpxform
        tqs = mtoq(transpose(tpxform))
        qs = qslerp(tqs, tuts, uts)
        pxform2 = transpose(qtom(qs))
        
        print, minmax(pxform1-pxform2)  ; minmax diff is -3e-7 and 3e-7.
        
        
    ; method 3: use 1/4 sec time resolution, interpolate (quadratic better than spline).
        tuts = smkarthm(tutr[0], tutr[1], dt0/4d, 'dx')
        tnrec = n_elements(tuts)
        
        tmp = time_string(tuts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+tuts-tuts[0]
        
        cspice_pxform, scid, 'GSE', tets, tpxform
        pxform3 = dblarr(3,3,nrec)
        for i = 0, 2 do pxform3[i,*,*] = sinterpol(reform(tpxform[i,*,*]), tuts, uts, /quadratic)
        print, minmax(pxform1-pxform3)

        
    ; use the matrix from method 2.
        eu = eu-smooth(eu, spinrate/dr0, /edge_truncate)
        ev = ev-smooth(ev, spinrate/dr0, /edge_truncate)
        pxform = pxform2    ; u:[0,*], v:[1,*], w:[2,*].
        ex = eu*pxform[0,0,*] + ev*pxform[1,0,*] + ew*pxform[2,0,*]
        ey = eu*pxform[0,1,*] + ev*pxform[1,1,*] + ew*pxform[2,1,*]
        ez = eu*pxform[0,2,*] + ev*pxform[1,2,*] + ew*pxform[2,2,*]

    store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
        {ytitle:'E GSE!C(mV/m)', labels:['Ex','Ey','Ez'], colors:[6,4,2], yrange:[-200,200]}

    ; convert egse to emgse.
    tvar = pre0+'emgse_sheng'
    wscgse = transpose(pxform[2,*,*])
    emgse = sgse2mgse([[ex],[ey],[ez]], ets, wsc = wscgse)
    store_data, tvar, uts, emgse, limits = $
        {ytitle:'E MGSE!C(mV/m)', labels:['Ex','Ey','Ez'], colors:[6,4,2], yrange:[-200,200]}



; **** use rbsp_uvw_to_mgse.
    tvar = pre0+'euvw'
    rbsp_uvw_to_mgse, tprobe, tvar, /no_spice_load
    tvar = pre0+'emgse_tplot'
    stplot_renew, pre0+'euvw_mgse', newname = tvar, /delete
    
    
; **** load E mgse from cdf file.
    tvar = pre0+'emgse'
    tmp = sread_rbsp_efw_l2(utr, probe = tprobe, type = 'esvy')
    dat = tmp.efield_mgse
    tuts = sfmepoch(tmp.epoch, 'unix')
    dat[*,0] = 0.   ; set x to 0.
    dat = sinterpol(dat, tuts, uts)
    store_data, tvar, uts, dat, limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:[6,4,2], labels:'MGSE '+['x','y','z']}


    vars = pre0+['emgse','emgse_sheng','emgse_tplot']
    options, vars, 'yrange', [-60,60]
    tplot, vars

end
