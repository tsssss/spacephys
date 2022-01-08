;+
; download rbsp 'raw' data, save to disk.
; save necessary data, e.g., the ones uses spice kernel.
; data that are easy to load are not saved, e.g., HOPE spec.
; data that are derived are left for a later routine, e.g., FAC, GSE fields.
;-
pro psbl_de_save_32hz_data, id, trange = utr, rootdir = rootdir

    if n_elements(id) eq 0 then message, 'no id ...'
    if n_elements(tprobe) eq 0 then tprobe = strmid(id,strlen(id)-1)
    if n_elements(utr) eq 0 then $
        utr = time_double(strmid(id,0,14),tformat='YYYY_MMDD_hhmm')+[-2,8]*60
    if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/psbl_de_32hz'

    datadir = rootdir+'/data'
    datavar = []

    
    ; expand utr to at least 10 min.
    theutr = utr
    utr = mean(utr) & utr = utr-(utr mod 60) & utr+= [-1,1]*0.5*600

    
    pre0 = 'rbsp'+tprobe+'_'
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    rgb = [6,4,2]
    dt0 = 86400d
    utr0 = utr-(utr mod dt0)+[0,dt0]
    timespan, utr0[0], dt0, /second
    
    
    
; **** load E field. 32Hz dE MGSE.
; rbspx_de_mgse, rbspx_de_survey_mgse.

    tvar = pre0+'de_mgse'
    dat = sread_rbsp_efw_l2(utr0, probes = tprobe, type = 'esvy')
    if size(dat,/type) ne 8 then message, 'no E data ...'
    tmp = sfmepoch(dat.epoch, 'unix')
    dr = sdatarate(tmp[sort(tmp)])
    uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
    dat = sinterpol(dat.efield_mgse, tmp, uts)
    dat[*,0] = 0.   ; set x to 0.
    store_data, tvar, uts, dat[*,[1,2]], limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
        labels:'MGSE '+['y','z'], constant:[-25,0,25]}
    datavar = [datavar,tvar]
    
    
    tvar = pre0+'de_survey_mgse'
    dat = sread_rbsp_efw_l3(utr0, probes = tprobe)
    if size(dat,/type) ne 8 then message, 'no E data ...'
    tmp = sfmepoch(dat.epoch, 'unix')
;    dr = sdatarate(tmp[sort(tmp)])
    dr = 11d
    uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
    tdat = sinterpol(dat.efield_inertial_frame_mgse, tmp, uts)
    tdat[*,0] = 0.   ; set x to 0.
    store_data, tvar, uts, tdat[*,[1,2]], limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
        labels:'MGSE '+['y','z'], constant:[-25,0,25]}
    datavar = [datavar,tvar]
    

; **** load position.
; rbspx_[mlt,lshell,ilat,pos_gse,mlat].

    tvar = pre0+'pos_gse'
    tdat = sinterpol(dat.pos_gse, tmp, uts)*re1
    store_data, tvar, uts, tdat, $
        {ytitle:'R GSE (Re)', colors:rgb, labels:'GSE '+['x','y','z']}
    datavar = [datavar,tvar]
    
    tvar = pre0+'mlt'
    tdat = sinterpol(dat.mlt_lshell_mlat, tmp, uts)
    store_data, tvar, uts, tdat[*,0], limits = {ytitle:'MLT (hr)'}
    datavar = [datavar,tvar]
    
    tvar = pre0+'lshell'
    store_data, tvar, uts, tdat[*,1], limits = {ytitle:'L-shell'}
    datavar = [datavar,tvar]
    
    tvar = pre0+'mlat'
    store_data, tvar, uts, tdat[*,2], limits = {ytitle:'MLat (deg)'}
    datavar = [datavar,tvar]


; **** load B field. l3, gse, emfisis, hires 171 MB for 1 day.
; rbspx_b_gse, rbspx_b.

    emfisis = sread_rbsp_emfisis_l3(utr0, probes = tprobe, $
        type = 'hires', coord = 'gse')
    if size(emfisis,/type) ne 8 then message, 'no B data ...'
    uts = sfmepoch(emfisis.epoch,'unix')
    store_data, pre0+'b_gse', uts, emfisis.mag, $
        limits = {ytitle:'B!C(nT)',colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B magnitude!C(nT)'}    
    datavar = [datavar,pre0+'b_gse']


; **** interpolate the high resolution fields to uniform time.
    get_data, pre0+'de_mgse', uts
    dr = sdatarate(uts)
    dr = 1d/32  ; 32 Hz.
    uts = smkarthm(min(uts), max(uts), dr, 'dx')
    nrec = n_elements(uts)
    vars = pre0+['de_mgse','b_gse','b']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], tuts, tdat
        store_data, vars[i], uts, sinterpol(tdat, tuts, uts)
    endfor


; **** omni Dst and AE.
; ae,dst.
    omni = sread_omni(utr0)
    uts = sfmepoch(omni.epoch,'unix')
    store_data, 'dst', uts, omni.sym_h, limits = {ytitle:'Dst (nT)'}
    store_data, 'ae', uts, omni.ae_index, limits = {ytitle:'AE (nT)'}
    datavar = [datavar,['dst','ae']]


; **** mapping coefficient, footpoint.
; rbspx_fpt_[mlon,mlat,mlt], rbspx_map_coef, rbspx_bmod_gsm.
    scalc_map_coef, pre0+'pos_gse', pre0+'b', model = 't89', $
        prefix = pre0, altitude = 110, coord = 'gse'
    
    tvar = pre0+'fpt_mlat'
    get_data, tvar, tmp, dat

;; this is a test, the plots show that the original mlat is ok.
;; no need to do this complicated interpolation.
;    options, tvar, 'labels', 'orig'
;    tutr = mean(utr)+1800*[-1,1]
;    idx = where(tmp ge tutr[0] and tmp le tutr[1])
;    hem = (dat[idx[0]] gt 0)? 1: -1
;    dat = sinterpmax(abs(dat[idx]))*hem
;    store_data, tvar+'_interp', tmp[idx], dat, $
;        limits = {labels:'interp'}
;    tvar = pre0+'fpt_mlat_compare'
;    store_data, tvar, data = pre0+'fpt_mlat'+['','_interp'], $
;        limits = {psym:-1,colors:[1,2]}

    tvar = pre0+'fpt_mlon'
    get_data, tvar, tmp, dat
    store_data, pre0+'fpt_mlt', tmp, $
        slon2lt(dat/15 ,stoepoch(tmp,'unix'),/mag), $  ; in [0,24].
        limits = {labels:'mapped'}
    options, pre0+'mlt', 'labels', 'in situ'
    tvar = pre0+'mlt_compare'
    store_data, tvar, data = pre0+['mlt','fpt_mlt'], $
        limits = {psym:-1,colors:[1,2]}

    datavar = [datavar,pre0+['fpt_mlon','fpt_mlat','fpt_mlt', $
        'map_coef','bmod_gsm']]
    
;    ofn = rootdir+'/psbl_de_32hz_mlt_compare_'+id+'.pdf'
;;    ofn = 0
;    tplot_options, 'xmargin', [10,10]
;    sgopen, ofn, xsize = 5, ysize = 3, /inch
;    device, decomposed = 0
;    loadct2, 43
;    vars = pre0+['fpt_mlat_compare','mlt_compare']
;;    vars = pre0+['mlt_compare']
;    tplot, vars, trange = mean(utr)+600*[-1,1]
;    sgclose
    
    
    
; **** trim data to tutr, 20 min long.
    tutr = mean(utr)+600*[-1,1]
    for j = 0, n_elements(datavar)-1 do begin
        get_data, datavar[j], tuts, tdat
        idx = where(tuts ge tutr[0] and tuts le tutr[1])
        store_data, datavar[j], tuts[idx], tdat[idx,*,*,*]
    endfor


; **** get w_gse, will be used to convert b/w GSE and mGSE.
; use the survey resolution is enough.
    get_data, pre0+'de_survey_mgse', uts
    nrec = n_elements(uts)
    wgse = dblarr(nrec,3)
;    rbsp_load_spice_kernels, probe = tprobe, trange = utr
    for i = 0, nrec-1 do begin
        tstr = time_string(uts[i], tformat = 'YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tstr, et
        cspice_pxform, 'RBSP'+tprobe+'_SCIENCE', 'GSE', et, pxform
        wgse[i,*] = pxform[2,*]
    endfor
    store_data, pre0+'w_gse', uts, wgse, limits = $
        {ytitle:'W GSE!Cunit vector', labels:'GSE '+['x','y','z'], colors:[6,4,2]}
    datavar = [datavar,pre0+'w_gse']

    
; **** save data to disk.
    ofn = datadir+'/psbl_de_32hz_data_'+id+'.tplot'
    tplot_save, datavar, filename = ofn
    
end

rootdir = shomedir()+'/psbl_de_32hz/'

idinfos = psbl_de_id('detect_step')
nidinfo = n_elements(idinfos)


;rbsp_load_spice_kernels, probe = 'a', /all
;rbsp_load_spice_kernels, probe = 'b', /all


for i = 0, nidinfo-1 do begin
    id = idinfos[i].id
    utr = idinfos[i].utr
    
;    if id ne '2013_0501_0736_b' then continue
        
    store_data, '*', /delete
    psbl_de_save_32hz_data, id, trange = utr, rootdir = rootdir
endfor

end
