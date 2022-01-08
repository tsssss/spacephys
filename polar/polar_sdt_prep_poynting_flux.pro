;+
; Type: procedure.
;
; Purpose: For Polar sdt, prepare E/B fields for poynting flux, remove spike,
;   interpolate to uniform data, rotate the field data from SPC (Spin Plane
;   Coordinates) to FAC (Field Aligned Coordinates).
;       {SPC. xy:away sun, z:north, 56:spin axis}{FAC. v:perp B, along v sense,
;   p:bxv, b:parallel B}. For noon-midnight meridian, xy~v, z~b, 56~p.
;       Interface: {out: 'po_b0_spc','po_b','po_db_spc','po_de_spc', $
;   'po_db_fac','po_de_fac','po_spc2fac','po_ilat','po_mlt','po_dis'}
;
; Parameters:
;   fn, in, string, required. Input filename, exported from sdt.
;
; Keywords:
;   tr, in, dblarr[2], optional. Time range for plot, optional but should
;       provide almost all the time.
;   cusptr, in, dblarr[2], optional. Time range for cusp entry and exist.
;   e56, in, boolean, optional. Set to keep e56.
;   edot0, in, boolean, optional. Set to calc e56 using E dot B = 0.
;   orootdir, in, string, optional. The directory to save the plots.
;
; Notes: Polar SDT E/B field don't have NaN. The main problems are: (1) the
;   spikes in B field; (2) how to separate the model B field from the total B 
;   field; (3) whether to believe E56, if not use Edot0 or set to 0.
;       The default behaviors for the problems are: (1) use sdespike to remove
;   spikes; (2) use T96 model data stored in file as 1st order approx, then use
;   scalcbg to remove the background from the 1st order dB; (3) set E56 = 0.
;       Other bahaviors are: (1) set e56 to use E56; (3) set edot0 to use 
;   E dot B = 0 to calc E56. These behaviors are not recommended to set.
;
; Dependence: slib, tdas.
;
; History: 2013-12-01, Sheng Tian, create.
;-
pro polar_sdt_prep_poynting_flux, fn, tr = tr, cusptr = cusptr, $
    e56 = e56, edot0 = edot0, orootdir = orootdir, noplot = noplot, $
    titpre = titpre, eventid = eventid, _extra = extra

    if file_test(fn) eq 0 then message, 'file does not exist ...'
    sdt = ssdtread(fn)
    pre = 'po_'
    
    if n_elements(orootdir) eq 0 then orootdir = shomedir()
    if ~file_test(orootdir,/directory) then file_mkdir, orootdir
    
    ; some quantities and settings.
    device, decomposed = 0 & loadct2, 43
    red = 6
    green = 4
    blue = 2
    black = 0
    rgb = [red,green,blue]   ; r,g,b in colortable 43,45.
    labfac = ['v','p','b']   ; for north sunward meridian cross
    labspc  = ['xy','z','56']       ; v~xy, p(vxb)~56, b~z.
    defactor = 1.3d ; correction to E due to shielding.
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', -1
    time_stamp, /off
    ct = 43
    posd = [0.15,0.10,0.9,0.45]     ; lower half.
    posu = [0.15,0.55,0.9,0.90]     ; upper half.
    posl = [0.10,0.1,0.30,0.9]      ; left.
    posm = [0.40,0.1,0.60,0.9]      ; middle.
    posr = [0.70,0.1,0.90,0.9]      ; right.

    if n_elements(titpre) eq 0 then titpre = ''

    
    ;**** get uniform time.
    maxnrec = 100000ul   ; prevent memory overflow.
    t0 = sdt.var.polar_b_spc_z.depend_0
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    dr = 1d/8   ; force to use 8 samples/sec.
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
    

    ;**** plot1: original field and spike removal.
    vars = ['bxy','bz','b56','b0xy','b0z','b056']
    options, vars, 'ytitle', '(nT)'
    ylim, ['bxy','b0xy'], min(b_spc[*,0])-20, max(b_spc[*,0])+20, 0
    ylim, ['bz','b0z'], min(b_spc[*,1])-20, max(b_spc[*,1])+20, 0
    ylim, ['b56','b056'], min(b_spc[*,2])-20, max(b_spc[*,2])+20, 0
    fn = orootdir+'/'+tstr+'_polar_b_despike.pdf'
    if ~keyword_set(noplot) then begin
        sgopen, fn, xsize = 5, ysize = 7, /inch
        sgindexcolor
        loadct2, ct
        vars = ['bxy','bz','b56']
        nvar = n_elements(vars)
        titl = 'Polar B_SPC, raw data'
        poss = sgcalcpos(nvar, position = posu)
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        vars = ['b0xy','b0z','b056']
        nvar = n_elements(vars)
        titl = 'Polar B_SPC, after despike'
        poss = sgcalcpos(nvar, position = posd)
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        xyouts, 0.5, 0.95, titpre+'Polar original and despiked B!DSPC!N', /normal, alignment = 0.5, charsize = 1.25
        sgclose
    endif
    vars = ['bxy','bz','b56','b0xy','b0z','b056']
    store_data, vars, /delete
    
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
    
    ;**** plot2: B model on B, dB, original dE and
    ; B total, dB 3D before correction, E56_dot0.
    fn = orootdir+'/'+tstr+'_polar_e56dot0_dbcorr.pdf'
    if ~keyword_set(noplot) then begin
        sgopen, fn, xsize = 10, ysize = 7, /inch
        sgindexcolor
        loadct2, ct
        ; plot B model on top of B total.
        store_data, 'b0xy', t0, [[bmod_spc[*,0]],[b_spc[*,0]],[bt96_spc[*,0]]], $
            limits = {labels:'Bxy'}
        store_data, 'b0z', t0, [[bmod_spc[*,1]],[b_spc[*,1]],[bt96_spc[*,1]]], $
            limits = {labels:'Bz'}
        store_data, 'b056', t0, [[bmod_spc[*,2]],[b_spc[*,2]],[bt96_spc[*,2]]], $
            limits = {labels:'B56'}
        vars = ['b0xy','b0z','b056'] & nvar = n_elements(vars)
        titl = 'Polar B model SPC'
        poss = posl & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'colors', [red,black,green]
        options, vars, 'ytitle', '(nT)'
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; plot dB SPC in components.
        store_data, 'dbxy', t0, db_spc[*,0], limits = {labels:'dBxy'}
        store_data, 'dbz', t0, db_spc[*,1], limits = {labels:'dBz'}
        store_data, 'db56', t0, db_spc[*,2], limits = {labels:'dB56'}
        vars = ['dbxy','dbz','db56'] & nvar = n_elements(vars)
        titl = 'Polar dB SPC'
        poss = posm & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(nT)'
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; plot dE SPC in components.
        store_data, 'dexy', t0, de_spc[*,0], limits = {labels:'dE0xy'}
        store_data, 'dez' , t0, de_spc[*,1], limits = {labels:'dE0z'}
        store_data, 'de56', t0, de560, limits = {labels:'dE056'}
        vars = ['dexy','dez','de56'] & nvar = n_elements(vars)
        titl = 'Polar dE SPC'
        poss = posr & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(mV/m)'
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete & de560 = 0
        ; plot dE56 dot0.
        store_data, 'de56_dot0', t0, de56_dot0
        vars = ['de56_dot0'] & nvar = n_elements(vars)
        titl = 'Polar dE56_dot0'
        poss = posr & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(mV/m)'
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; plot dB SPC before correction.
        store_data, 'db', t0, db_spc, limits = {ytitle:'(nT)', $
            labels:labspc, colors:rgb}
        vars = ['db'] & nvar = n_elements(vars)
        titl = 'Polar dB SPC 3D'
        poss = posm & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(nT)'
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; plot B magnitude.
        tmp = sqrt(total(bmod_spc^2,2))
        tmp2 = sqrt(total(bt96_spc^2,2))
        store_data, 'btotal', t0, [[tmp],[btotal],[tmp2]], limits = $
            {colors:[red,black,green], labels:['Bsmth','in situ','T96']}
        vars = ['btotal'] & nvar = n_elements(vars)
        titl = 'Polar B magnitude'
        poss = posl & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(nT)'
        options, vars, 'ystyle', 1
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        xyouts, 0.5, 0.95, titpre+'Polar original B and E, B!Dmodel!N dB separation, dE!Ddot0!N', /normal, alignment = 0.5, charsize = 1.25
        sgclose
    endif
    
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
    
    ft  = sdt.var.polarpositiongsm_x.depend_0
    fx = sdt.var.polarpositiongsm_x.value
    fy = sdt.var.polarpositiongsm_y.value
    fz  = sdt.var.polarpositiongsm_z.value
    pos = [[fx],[fy],[fz]]
    store_data, pre+'pos_gsm', t1, pos, limits={ytitle:'R GSM (Re)'}
    
    ;**** rotate from SPC to FAC.
    polar_sdt_spc2fac, bmod_spc, db_spc, de_spc, db_fac, de_fac, b1 = lon, b2 = lat
    store_data, pre+'db_fac', t0, db_fac, $
        limits = {ytitle:'dB FAC!C(nT)', labels:labfac, colors:rgb}
    store_data, pre+'de_fac', t0, de_fac, $
        limits = {ytitle:'dE FAC!C(mV/m)', labels:labfac, colors:rgb}
    store_data, pre+'spc2fac', t0, [[lon],[lat]]*(180/!dpi), limits = $
        {ytitle:'SPC2FAC!C(deg)', labels:['lon','lat'], colors:rgb[0:1]}
        
;---Event specific treatment.
    utr1 = time_double(['1998-10-29/10:16','1998-10-29/10:18:15'])
    get_data, pre+'db_fac', uts, dat
    idx = where(uts ge utr1[0] and uts le utr1[1], cnt)
    if cnt ne 0 then begin
        dat[idx,*] = 0
        store_data, pre+'db_fac', uts, dat
    endif
    
    utr2 = time_double(['1998-10-30/03:20:13','1998-10-30/03:20:34'])
    get_data, pre+'db_fac', uts, dat
    idx = where(uts ge utr2[0] and uts le utr2[1], cnt)
    if cnt ne 0 then begin
        dat[idx,*] = 0
        store_data, pre+'db_fac', uts, dat
    endif
    
    utr3 = time_double(['1998-10-30/03:32:03','1998-10-30/03:32:12'])
    get_data, pre+'db_fac', uts, dat
    idx = where(uts ge utr3[0] and uts le utr3[1], cnt)
    if cnt ne 0 then begin
        dat[idx,*] = 0
        store_data, pre+'db_fac', uts, dat
    endif
        
    ; plot3: show dB and dE rotated from SPC to FAC.
    fn = orootdir+'/'+tstr+'_polar_spc2fac.pdf'
    if ~keyword_set(noplot) then begin
        sgopen, fn, xsize = 7, ysize = 7, /inch
        sgindexcolor
        loadct2, ct
        ; dE FAC.
        vars = 'de_fac'+labfac
        nvar = n_elements(vars)
        stplot_split, pre+'de_fac', newname = vars
        tminmax = [0d,0d]
        for i = 0, nvar-1 do begin
            twd = 60d/dr         ; 60 sec smooth.
            get_data, vars[i], t0, tmp
            dat = smooth(tmp, 60d/dr, /edge_mirror)
            store_data, vars[i], t0, [[tmp],[dat]], $
                limits = {colors:[0,red], labels:'dE'+labfac[i]}
            tminmax = minmax([tminmax,dat])
        endfor
        titl = 'Polar dE FAC'
        poss = posu & poss[2] = 0.45 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(mV/m)'
        options, vars, 'yrange', tminmax
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; dB FAC.
        vars = 'db_fac'+labfac
        nvar = n_elements(vars)
        options, pre+'db_fac', 'labels', 'dB'+labfac
        stplot_split, pre+'db_fac', newname = vars
        get_data, pre+'db_fac', t0, tmp
        tminmax = minmax(tmp)
        titl = 'Polar dB FAC'
        poss = posu & poss[0] = 0.60 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(nT)'
        options, vars, 'yrange', tminmax
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; dE SPC.
        vars = 'de_spc'+labspc
        nvar = n_elements(vars)
        stplot_split, pre+'de_spc', newname = vars
        tminmax = [0d,0d]
        for i = 0, nvar-1 do begin
            twd = 60d/dr         ; 60 sec smooth.
            get_data, vars[i], t0, tmp
            dat = smooth(tmp, 60d/dr, /edge_mirror)
            store_data, vars[i], t0, [[tmp],[dat]], $
                limits = {colors:[0,red], labels:'dE'+labspc[i]}
            tminmax = minmax([tminmax,dat])
        endfor
        titl = 'Polar dE SPC'
        poss = posd & poss[2] = 0.45 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(mV/m)'
        options, vars, 'yrange', tminmax
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        ; dB SPC.
        vars = 'db_spc'+labspc
        nvar = n_elements(vars)
        options, pre+'db_spc', 'labels', 'dB'+labspc
        stplot_split, pre+'db_spc', newname = vars
        get_data, pre+'db_spc', t0, tmp
        tminmax = minmax(tmp)
        titl = 'Polar dB SPC'
        poss = posd & poss[0] = 0.60 & poss = sgcalcpos(nvar, position = poss)
        options, vars, 'ytitle', '(nT)'
        options, vars, 'yrange', tminmax
        tplot, vars, title = titl, trange = tr, position = poss, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        store_data, vars, /delete
        xyouts, 0.5, 0.95, titpre+'Polar SPC to FAC rotation of dE and dB', /normal, alignment = 0.5, charsize = 1.25
        sgclose
    endif

end

;fn = sdiskdir('Works')+'/data/cusp/po_sdt_fld_1998_0925_05.sdt'
fn = sdiskdir('Works')+'/data/cusp/po_sdt_fld_1998_1001_02.sdt'
;fn = sdiskdir('Works')+'/data/cusp/po_sdt_fld_1998_1002_16.sdt'
;fn = sdiskdir('Works')+'/data/cusp/po_sdt_fld_1999_1022_14.sdt'
polar_sdt_prep_poynting_flux, fn
end
