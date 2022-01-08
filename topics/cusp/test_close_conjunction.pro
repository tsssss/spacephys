pro map2ilat, var, ilatvar, ilat0s

    get_data, var, t0, dat, val
    get_data, ilatvar, tmp, ilat1s
    
    ; original ilat map to data's ut.
    ilat1s = interpol(ilat1s, tmp, t0)
    
    store_data, var, ilat1s, dat, val
end




ids = cusp_id('close_conjunction')

id = ids[4]

; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = id)



; load polar field.
    fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
    if file_test(fn) eq 0 then message, 'file does not exist ...'
    sdt = ssdtread(fn)
    pre = 'po_'

    rgb = sgcolor(['red','green','blue'])
    labfac = ['v','p','b']   ; for north sunward meridian cross
    labspc  = ['xy','z','56']       ; v~xy, p(vxb)~56, b~z.
    defactor = 1.3d ; correction to E due to shielding.
    
    ;**** get uniform time.
    maxnrec = 100000ul   ; prevent memory overflow.
    t0 = sdt.var.polar_b_spc_z.depend_0
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    tmp = minmax(sdt.var.polar_e_spc_z.depend_0)
    t0 = smkarthm(max([tmp[0],t0[0]]),min([tmp[1],t0[nrec-1]]), dr, 'dx')
    nrec = n_elements(t0)
    tstr = time_string(t0[0], tformat='YYYY_MMDD')
    if n_elements(eventid) ne 0 then tstr = eventid
    print, 'data rate: ', dr
    
    ;**** original b field and spike removal.
    ft  = sdt.var.polar_b_spc_z.depend_0
    fxy = sdt.var.polar_b_spc_x_y.value
    f56 = sdt.var.polar_b_spc_56.value
    fz  = sdt.var.polar_b_spc_z.value
    b_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    b_spc = sinterpol(b_spc, ft, t0)
    ; raw B field data, only interpolated to uniform time.
    store_data, 'bxy', t0, b_spc[*,0], limits = {labels:'B0xy'}
    store_data, 'bz' , t0, b_spc[*,1], limits = {labels:'B0z'}
    store_data, 'b56', t0, b_spc[*,2], limits = {labels:'B056'}
    ; Bxy and Bz have discontinuities and spikes, B56 has short spikes.
    sdespike, t0, b_spc, _extra = extra
    store_data, 'b0xy', t0, b_spc[*,0], limits = {labels:'Bxy'}
    store_data, 'b0z' , t0, b_spc[*,1], limits = {labels:'Bz'}
    store_data, 'b056', t0, b_spc[*,2], limits = {labels:'B56'}

    ;**** total b field. 'po_b'
    btotal = sqrt(total(b_spc^2,2))
    store_data, pre+'b', t0, btotal, limits = {ytitle:'B mag!C(nT)', ynozero:1}
    
    ;**** t96 model.
    ft  = sdt.var.polar_model_b_t96_spc_z.depend_0
    fxy = sdt.var.polar_model_b_t96_spc_x_y.value
    f56 = sdt.var.polar_model_b_t96_spc_56.value
    fz  = sdt.var.polar_model_b_t96_spc_z.value
    bt96_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    bt96_spc = sinterpol(bt96_spc, ft, t0)
    
    ;**** model b field. 'po_b0_spc'.
    db_spc = b_spc-bt96_spc
    bmod_spc = scalcbg(db_spc)+bt96_spc
    store_data, pre+'b0_spc', t0, bmod_spc, $
        limits = {ytitle:'B model SPC!C(nT)', labels:labspc, colors:rgb}
        
    ;**** db field. 'po_db_spc'.
    db_spc = b_spc-bmod_spc
;    idx = where(abs(db_spc) gt 500, cnt)
;    if cnt ne 0 then db_spc[idx] = 0
    store_data, pre+'db_spc', t0, db_spc, $
        limits = {ytitle:'dB SPC!C(nT)', labels:labspc, colors:rgb}
    store_data, pre+'b0_spc', t0, bmod_spc
    
    ;**** de field. 'po_de_spc'.
    ft  = sdt.var.polar_e_spc_z.depend_0
    fxy = sdt.var.polar_e_spc_x_y.value
    f56 = sdt.var.polar_e_spc_56.value
    fz  = sdt.var.polar_e_spc_z.value
    de_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    de_spc = sinterpol(de_spc, ft, t0)*defactor
    de560 = de_spc[*,0]
    ; e56_dot0.
    de56_dot0 = (b_spc[*,0]*de_spc[*,0]+b_spc[*,1]*de_spc[*,1])/b_spc[*,2]
    idx = where(abs(b_spc[*,2]/btotal) le 0.2, cnt)
    if cnt ne 0 then de56_dot0[idx] = !values.f_nan
    ; throw e56 by default, set e56 to keep e56, set edot0 to calc e56.
    if ~keyword_set(e56) then de_spc[*,2] = 0
    if keyword_set(edot0) then de_spc[*,2] = de56_dot0
    store_data, pre+'de_spc', t0, de_spc, $
        limits = {ytitle:'dE SPC!C(mV/m)', labels:labspc, colors:rgb}
    store_data, 'de56', t0, de_spc[*,2]

    ;**** ilat, mlt, dis.
    t1  = sdt.var.polarinvariantlatitude.depend_0
    mlat = sdt.var.polarmaglat.value
    store_data, pre+'mlat', data = {x:t1, y:mlat}, limits = {ytitle:'MLat'}
    ilat = sdt.var.polarinvariantlatitude.value
    tmp = where(mlat lt 0)
    if tmp[0] ne -1 then ilat[tmp] *= -1
    store_data, pre+'ilat', data = {x:t1, y:ilat}, limits = {ytitle:'ILat'}
    mlt  = sdt.var.polarmlt.value
    store_data, pre+'mlt', data = {x:t1, y:mlt}, limits = {ytitle:'MLT'}
    dis  = sdt.var.polarspcraftdist.value
    store_data, pre+'dis', data = {x:t1, y:dis}, limits = {ytitle:'Dist (Re)'}
    
    ;**** rotate from SPC to FAC.
    polar_sdt_spc2fac, bmod_spc, db_spc, de_spc, db_fac, de_fac, ilat, lon, lat
    store_data, pre+'db_fac', t0, db_fac, $
        limits = {ytitle:'dB FAC!C(nT)', labels:labfac, colors:rgb}
    store_data, pre+'de_fac', t0, de_fac, $
        limits = {ytitle:'dE FAC!C(mV/m)', labels:labfac, colors:rgb}
    store_data, pre+'spc2fac', t0, [[lon],[lat]]*(180/!dpi), limits = $
        {ytitle:'SPC2FAC!C(deg)', labels:['lon','lat'], colors:rgb[0:1]}


; load fast field.
    fainfo = loginfo.fast
    fatr = loginfo.fast.plot_time
    trplot = fatr+(fatr[1]-fatr[0])*[-1,1]
    fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_'+ $
        string(fainfo.orbit,format='(I05)')+'.tplot' 
    if file_search(fn) eq '' then $
        fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_'+$
            string(fainfo.orbit,format='(I05)')+'.tplot'
    tplot_restore, filename = fn
    pre = 'fa_'
    suf = ''

    ; truncate to given trange.
    tr = minmax(fatr)
    vars = ['dB_fac_v','B_model','E_ALONG_V','E_NEAR_B','ilat','mlt','dis']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, dat
        idx = where(t0 ge tr[0] and t0 le tr[1], cnt)
        if cnt eq 0 then message, 'invalid time range ...'
        store_data, vars[i], t0[idx], dat[idx,*]
    endfor

    ; get uniform time.
    maxnrec = 100000ul   ; prevent memory overflow.
    get_data, 'dB_fac_v', t0, dbfac
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    t0 = smkarthm(t0[0], t0[nrec-1], dr, 'dx')
    nrec = n_elements(t0)
    print, 'data rate: ', dr
    
    ; db field.
    get_data, 'dB_fac_v', data = tmp
    db_facv = sinterpol(tmp.y, tmp.x, t0)
;    ; remove background field as in polar.
;    db_bg = scalcbg(db_facv)
;    db_facv = db_facv-db_bg
    store_data, pre+'db_fac'+suf, data = {x:t0, y:db_facv}, $
        limits = {ytitle: 'dB FAC!C(nT)', labels:labfac, colors:rgb}
            
    ; model b field in GEI.
    get_data, 'B_model', tmp, bmod_gei
    bmod_gei = sinterpol(bmod_gei, tmp, t0)
    store_data, pre+'b0_gei', t0, bmod_gei, $
        limits = {ytitle:'B model GEI!C(nT)', labels:['x','y','z'], colors:rgb}

    ; total b field.
    get_data, pre+'db_fac'+suf, t0, btotal
    btotal[*,2]+= sqrt(total(bmod_gei^2,2)) & btotal = sqrt(total(btotal^2,2))
    store_data, pre+'b', t0, btotal, limits = $
        {ytitle:'B mag!C(nT)', ynozero:1}
    
    ; de field.
    get_data, 'E_ALONG_V', data = tmp & dev = sinterpol(tmp.y, tmp.x, t0)
    get_data, 'E_NEAR_B', data = tmp & deb = sinterpol(tmp.y, tmp.x, t0)
    ; fix sign error in SDT when FAST in south.
    get_data, 'ilat', tmp, dat
    idx = where(dat le 0, cnt)
    if cnt ne 0 then begin
        idx = where(t0 ge min(tmp[idx]) and t0 le max(tmp[idx]))
        dev[idx] = -dev[idx]
        deb[idx] = -deb[idx]
    endif
    de_facv = [[dev],[dblarr(nrec)],[deb]]      ; [v,bxv,b]
    store_data, pre+'de_fac'+suf, data = {x:t0, y:de_facv}, $
        limits = {ytitle: 'dE FAC!C(mV/m)', labels:labfac, colors:rgb}
    
    ; ilat, mlt, dis.
    stplot_renew, 'ilat', newname = pre+'ilat'
    stplot_renew, 'mlt', newname = pre+'mlt'
    stplot_renew, 'dis', newname = pre+'dis'
    
    


    ; get common ilat.
    potr = loginfo.polar.plot_time
    fatr = loginfo.fast.plot_time
    range = []
    get_data, 'po_ilat', tmp, dat
    range = minmax([range,dat])
    get_data, 'fa_ilat', tmp, dat
    range = minmax([range,dat])
    range = [floor(range[0]),ceil(range[1])]
    ilat0s = smkarthm(range[0],range[1],0.1,'dx')
    
    ; map efluxes based on ilat.
    pre = 'po_'
    ; ut.
    get_data, pre+'ilat', t0, dat
    store_data, pre+'ut', dat, t0
    ; convert.
    vars = pre+['db_fac','de_fac','b','mlt','dis']
    ilatvar = pre+'ilat'
    nvar = n_elements(vars)
    for i = 0, nvar-1 do map2ilat, vars[i], ilatvar, ilat0s
    ; ilat.
    get_data, pre+'ilat', tmp, dat
    store_data, pre+'ilat', dat, dat
    
    ; map efluxes based on ilat.
    pre = 'fa_'
    ; ut.
    get_data, pre+'ilat', t0, dat
    store_data, pre+'ut', dat, t0
    ; convert.
    vars = pre+['db_fac','de_fac','b','mlt','dis']
    ilatvar = pre+'ilat'
    nvar = n_elements(vars)
    for i = 0, nvar-1 do map2ilat, vars[i], ilatvar, ilat0s
    ; ilat.
    get_data, pre+'ilat', tmp, dat
    store_data, pre+'ilat', dat, dat

    get_data, 'po_b', pot0, pob
    get_data, 'fa_b', tmp, fab
    fab = interpol(fab, tmp, pot0)
    get_data, 'fa_dis', tmp, far
    far = interpol(far, tmp, pot0)
    get_data, 'po_dis', tmp, por
    por = interpol(por, tmp, pot0)
    fabmap = fab*(far/por)^3
    store_data, 'po_b', pot0, [[pob],[fabmap]]
    options, 'po_b', 'colors', sgcolor(['black','red'])
    options, 'po_b', 'labels', ['Polar B','FAST B!C  mapped']
    options, 'fa_b', 'labels', 'FAST B'
    
    
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 16
    tplot_options, 'version', 2
    
    poss = sgcalcpos(2, ypad = 15, lmarg = 20, bmarg = 15)
    
    ofn = shomedir()+'/cusp_close_conjunction_'+id+'.pdf
;    ofn = 0
    sgopen, ofn, xsize = 7, ysize = 9, /inch
    !p.charsize = 1
    !p.color = sgcolor('black')
    
    tr = minmax([potr, fatr])
    tr = -[66,80]    
    
    vars = ['po_db_fac','po_de_fac','po_b']
    nvar = n_elements(vars)
    labs = 'po_'+['ilat','mlt','dis']
    tpos = poss[*,0]
    tpos = sgcalcpos(nvar, position = tpos)
    
    ; trim polar data to tr.
    for i = 0, nvar-1 do begin
        get_data, vars[i], t0, dat
        idx = where(t0 eq min(t0))
        store_data, vars[i], t0[0:idx], dat[0:idx,*]
    endfor
    
    tplot, vars, trange = tr, /noerase, position = tpos, var_label = labs, title = 'Polar field '+id, uttick = 'po_ut'
    
    vars = ['fa_db_fac','fa_de_fac','fa_b']
    nvar = n_elements(vars)
    nvar = n_elements(vars)
    labs = 'fa_'+['ilat','mlt','dis']
    tpos = poss[*,1]
    tpos = sgcalcpos(nvar, position = tpos)
    tplot, vars, trange = tr, /noerase, position = tpos, var_label = labs, title = 'FAST field '+ id, uttick = 'fa_ut'
    
    sgclose
end
