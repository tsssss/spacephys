
pro calc_b_spin_plane_angle, ids

    deg = 180d/!dpi
    
    nid = n_elements(ids)
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'

    for i = 0, nid-1 do begin
        ; load data.
        id = ids[i]
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'

        loginfo = cusp_read_conjun_list(logfile, event = id)
        cusptr = loginfo.polar.cusp_time
        if cusptr[1] lt cusptr[0] then cusptr[1]+= 86400d

        sdt = ssdtread(fn)
        pre = 'po_'

        if n_elements(orootdir) eq 0 then orootdir = shomedir()+'/polar_cap_deb'
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
        tplot_options, 'labflag', 1
        tplot_options, 'constant', 0
        time_stamp, /off
        ct = 43
        posd = [0.15,0.10,0.9,0.45]     ; lower half.
        posu = [0.15,0.55,0.9,0.90]     ; upper half.
        posl = [0.10,0.1,0.30,0.9]      ; left.
        posm = [0.40,0.1,0.60,0.9]      ; middle.
        posr = [0.70,0.1,0.90,0.9]      ; right.
            

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

        ;**** angle between background B and spin plane.
        ang = atan(b_spc[*,2],sqrt(b_spc[*,0]^2+b_spc[*,1]^2))*deg
        store_data, pre+'bspin_ang', t0, ang, limits = {ytitle:'B/Spin Angle (Deg)'}

        tvar = pre+'de_fac'
        options, tvar, 'yrange', [-25,25]
        options, tvar, 'ystyle', 1
        
        tvar = pre+'bspin_ang'
        options, tvar, 'yrange', [-5,5]
        options, tvar, 'ystyle', 1
        

        ofn = orootdir+'/polar_cap_deb_'+id+'.pdf'
        sgopen, ofn, xsize = 8.5, ysize = 11, /inch
        sgindexcolor
        loadct2, ct
        vars = pre+['b0_spc','db_fac','de_fac','bspin_ang']
        labs = pre+['ilat','mlt','dis']
        titl = 'Polar Field'
        tplot, vars, var_label = labs, title = titl, trange = tr, /noerase
        if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
        sgclose
    endfor

end

ids = cusp_id('polar_cap_deb')
calc_b_spin_plane_angle, ids
end
