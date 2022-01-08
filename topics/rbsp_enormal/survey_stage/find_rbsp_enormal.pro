; rbspx_de_fac, rbspx_db_fac.
pro prep_rbsp_ebfield, probes

    nprobe = n_elements(probes)
    re = 6374d & re1 = 1d/re

    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        ; ** load e in mgse.
        rbsp_load_efw_waveform_l2, probe = probes[i]
        get_data, pre0+'efw_e-spinfit-mgse_efield_spinfit_mgse', data = tmp
        tmp.y[*,0] = 0      ; throw e_spin.
        store_data, pre0+'e_mgse', data = tmp
        ; convert from mgse to gse.
        rbsp_mgse2gse, pre0+'e_mgse', newname = pre0+'e_gse', $
            probe = probes[i], /no_spice_load
        ; delete vars.
        store_data, pre0+'efw_*', /delete
        
        ; ** load b in gse.
        rbsp_load_emfisis, probe = probes[i], coord = 'gse'
        get_data, pre0+'emfisis_l3_4sec_gse_Mag', data = tmp
        store_data, pre0+'b_gse', data = tmp
        store_data, pre0+'b_mag', data = {x:tmp.x, y:sqrt(total(tmp.y^2,2))}
        
        ; delete vars.
        store_data, pre0+'emfisis_*', /delete
        store_data, pre0+'efw_*', /delete
    endfor

    ; **** get model b field.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'b_gse', data = tmp
        bmod = tmp.y
        len = 1200D/sdatarate(tmp.x)    ; 20 min.
        for j = 0, 2 do bmod[*,j] = smooth(bmod[*,j], len)
        store_data, pre0+'bmod_gse', data = {x:tmp.x, y:bmod}
        store_data, pre0+'db_gse', data = {x:tmp.x, y:tmp.y-bmod}
    endfor

    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        rbsp_load_spice_state, probe = probes[i], $
            coord = 'gsm', /no_spice_load
        get_data, pre0+'state_pos_gsm', data = tmp
        tmp.y *= re1        ; from km to re.
        store_data, pre0+'pos_gsm', data = tmp
        ; delete vars.
        store_data, pre0+'state_*', /delete
    endfor
    
    ; get uniform time.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'e_gse', tmp, dat
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')
        ut0 = uts[0]-(uts[0] mod 86400)
        idx = where(uts ge ut0 and uts lt ut0+86400)    ; only 1 day.
        uts = uts[idx]
        dat = sinterpol(dat, tmp, uts)
        store_data, pre0+'e_gse', uts, dat
    endfor

    ; *** decompose e field, ex is along b, ey is east-west, ez is up.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'e_gse', data = de
        get_data, pre0+'db_gse', data = db
        get_data, pre0+'bmod_gse', data = tmp
        bmod = sinterpol(tmp.y, tmp.x, de.x)
        db = sinterpol(db.y, db.x, de.x)
        bhat = sunitvec(bmod)
        
        p = atan(bhat[*,1],bhat[*,0])
        cosp = cos(p) & sint = bhat[*,2]
        cost = bhat[*,0]/cosp & sinp = bhat[*,1]/cost
        
        vec = de.y
        x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
        y =       -sinp*vec[*,0] + cosp*vec[*,1]
        z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
        store_data, pre0+'de_fac', data = {x:de.x, y:[[x],[y],[z]]}
        ; de_dot0.
        if keyword_set(e_dot0) then begin
            de.y[*,0] = (de.y[*,1]*bmod[*,1]+de.y[*,2]*bmod[*,2])*$
                (-1/bmod[*,0])
            store_data, pre0+'de_gse_dot0', data = {x:de.x, y:de}
            vec = de.y
            x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
            y =       -sinp*vec[*,0] + cosp*vec[*,1]
            z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
            store_data, pre0+'de_fac_dot0', data = {x:de.x, y:[[x],[y],[z]]}
        endif
        vec = db
        x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
        y =       -sinp*vec[*,0] + cosp*vec[*,1]
        z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
        store_data, pre0+'db_fac', data = {x:de.x, y:[[x],[y],[z]]}

        vars = pre0+['db_fac','de_fac']
        options, vars, 'colors', [6,4,2]
        options, vars, 'labels', ['b','e','n']
    endfor

end


pro find_rbsp_enormal, tprobe, utr, bmax, emax

    ; initial settings.
    rbsp_efw_init
    cdf_leap_second_init
    rbsp_emfisis_init
    spc = '    '
    dt = 86400d
    rootdir = shomedir()+'/rbsp_de/rbsp'+tprobe

    if n_elements(bmax) eq 0 then bmax = 4e3    ; nT. threshold for finding eclipse.
    if n_elements(emax) eq 0 then emax = 10     ; mV/m. threshold for e normal.
    if n_elements(tmax) eq 0 then tmax = 300    ; sec. threshold for large field duration.
    
    
    logfn = rootdir+'/list_rbsp'+tprobe+'_large_de.log'
    if file_test(logfn) eq 0 then begin
        tmp = file_dirname(logfn)
        if file_test(tmp,/directory) eq 0 then file_mkdir, tmp
        openw, lun, logfn, /get_lun
        tit = '   event id   '+spc+'  start & end time  '+spc+ $
            'dE max'+spc+' dE normal '+spc+'  dE perp  '+spc+'  dE para  '
        printf, lun, tit
        tit = 'YYYY_MMDD_hhmm'+spc+'hh:mm:ss'+spc+'hh:mm:ss'+spc+ $
            '(mV/m)'+spc+'  min,max  '+spc+'  min,max  '+spc+'  min,max  '
        printf, lun, tit
        tit = '--------------'+spc+'--------'+spc+'--------'+spc+ $
            '------'+spc+'-----------'+spc+'-----------'+spc+'-----------'
        printf, lun, tit
        free_lun, lun
    endif

;    utr[0]-= dt
    
    ; time range to check.
    etr = stoepoch(utr, 'unix')
    utr0 = utr

    rbx = 'rbsp'+tprobe
    pre0 = rbx+'_'

    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10

    ; set up the data buffer.
    uts = []
    defac = []
    dbfac = []
    bmag = []
    emag = []
    posgsm = []

    et1 = sepochfloor(etr[0], 'dy')
    ut1 = sfmepoch(et1, 'unix')



    ; loop through each orbit. 
    go = 1
    load = 1    ; 1 for load a new day into buffer.
    while go do begin
        ; check if need new buffer.
        if load eq 1 then begin
            ut1 = ut1-(ut1 mod dt)
            print, 'loading data '+time_string(ut1)+' ...'
            ; load data for a whole day.
            timespan, ut1, 1, /day
            rbsp_load_spice_kernels
            prep_rbsp_ebfield, tprobe

            get_data, pre0+'de_fac', t0, tmp & defac = [defac,tmp]
            emag = sqrt(total(defac^2,2))
            get_data, pre0+'db_fac', t0, tmp & dbfac = [dbfac,tmp]
            uts = [uts,t0]
            get_data, pre0+'b_mag', t1, tmp
            tmp = interpol(tmp, t1, t0)
            bmag = [bmag,tmp]
            get_data, pre0+'pos_gsm', t1, tmp
            tmp = sinterpol(tmp, t1, t0)
            posgsm = [posgsm,tmp]
        endif

        ; prepare next round. idx0, idx1.
        idx = where(bmag ge bmax, cnt)
        if cnt eq 0 then begin
            ut1 = ut1+dt
            load = 1
            continue
        endif
        idx1 = idx[0]
        idx = where(bmag lt bmax, cnt)
        if cnt eq 0 then begin
            ut1 = ut1+dt
            load = 1
            continue
        endif
        idx0 = idx[0]
        if idx0 ge idx1 then begin
            idx = where(bmag[idx0:*] ge bmax, cnt)
            if cnt eq 0 then begin
                ut1 = ut1+dt
                load = 1
                continue
            endif
            idx1 = idx[0]+idx0
        endif
            

        load = 0
        dr = sdatarate(uts)
        
        tuts = uts[idx0:idx1]
        tdefac = defac[idx0:idx1,*]
        tdbfac = dbfac[idx0:idx1,*]
        temag = emag[idx0:idx1]
        tbmag = bmag[idx0:idx1]
        tposgsm = posgsm[idx0:idx1,*]
        
        uts = uts[idx1+1:*]
        defac = defac[idx1+1:*,*]
        dbfac = dbfac[idx1+1:*,*]
        bmag = bmag[idx1+1:*]
        emag = emag[idx1+1:*]
        posgsm = posgsm[idx1+1:*,*]

        ; check large e field.
        idx = where(temag ge emax, cnt)
        eidx0 = 0
        eidx1 = 0
        while eidx0 lt cnt-1 do begin
            tmp = eidx0
            for i = eidx0, cnt-2 do if idx[i+1]-idx[i] gt tmax/dr then break  ; about 5 min.
            eidx1 = i
            if eidx1 eq eidx0 then begin    ; 1 point, interpreted as spike.
                eidx0 = eidx1+1
                continue
            endif
            
            tidx0 = idx[eidx0]
            tidx1 = idx[eidx1]
            
            ; generate plots to file.
            ofn = rootdir+time_string(tuts[tidx0], tformat='/YYYY/MM/')+ $
                pre0+'large_de_'+ $
                time_string(tuts[tidx0], tformat='YYYY_MMDD_hhmm')+'.pdf'
            tutr = 0.5*(tuts[tidx0]+tuts[tidx1])+[-1,1]*tmax*3
            
            tvar = pre0+'de_fac'
            store_data, tvar, tuts, tdefac
            options, tvar, 'ytitle', 'dE FAC!C(mV/m)'
            
            tvar = pre0+'db_fac'
            store_data, tvar, tuts, tdbfac
            options, tvar, 'ytitle', 'dB FAC!C(nT)'
            
            tvar = pre0+'b_mag'
            store_data, tvar, tuts, tbmag
            options, tvar, 'ynozero', 1
            options, tvar, 'ytitle', 'B mag (nT)'
            tvar = pre0+['de_fac','db_fac','pos_gsm']
            options, tvar, 'colors', sgcolor(['red','green','blue'])
            
            tvar = pre0+'pos_gsm'
            store_data, tvar, tuts, tposgsm
            stplot_split, tvar, newname = tvar+'_'+['x','y','z']
            tvar = pre0+'pos_gsm_x'
            options, tvar, 'ytitle', 'X GSM (Re)'
            tvar = pre0+'pos_gsm_y'
            options, tvar, 'ytitle', 'Y GSM (Re)'
            tvar = pre0+'pos_gsm_z'
            options, tvar, 'ytitle', 'Z GSM (Re)'

            vars = pre0+['de_fac','db_fac','b_mag']
            nvar = n_elements(vars)
            labs = pre0+['pos_gsm_'+['z','y','x']]
            titl = 'RBSP-'+strupcase(tprobe)+' dE,dB survey '+$
                time_string(tuts[tidx0])
            poss = sgcalcpos(nvar, lmarg = 12)
            sgopen, ofn, xsize = 8, ysize = 11.5, /inch
            tplot, vars, var_label = labs, trange = tutr, title = titl, /noerase
            sgclose
            
            
            cmd = time_string(tuts[tidx0], tformat='YYYY_MMDD_hhmm')+spc
            cmd+= time_string(tuts[tidx0], tformat='hh:mm:ss')+spc
            cmd+= time_string(tuts[tidx1], tformat='hh:mm:ss')+spc
            tmp = sgnum2str(max(temag[tidx0:tidx1]),nsgn=3)
            for i = 1, 6-strlen(tmp) do tmp = ' '+tmp
            cmd = cmd+tmp+spc
            tmp = sgnum2str(min(tdefac[tidx0:tidx1,2]),nsgn=2)+','+$
                sgnum2str(max(tdefac[tidx0:tidx1,2]),nsgn=2)
            for i = 1, 11-strlen(tmp) do tmp = ' '+tmp
            cmd = cmd+tmp+spc
            tmp = sgnum2str(min(tdefac[tidx0:tidx1,1]),nsgn=2)+','+$
                sgnum2str(max(tdefac[tidx0:tidx1,1]),nsgn=2)
            for i = 1, 11-strlen(tmp) do tmp = ' '+tmp
            cmd = cmd+tmp+spc
            tmp = sgnum2str(min(tdefac[tidx0:tidx1,0]),nsgn=2)+','+$
                sgnum2str(max(tdefac[tidx0:tidx1,0]),nsgn=2)
            for i = 1, 11-strlen(tmp) do tmp = ' '+tmp
            cmd = cmd+tmp+spc
            
            openw, lun, logfn, /get_lun, /append
            printf, lun, cmd
            free_lun, lun

            eidx0 = eidx1+1
        endwhile
        
        ; trim the near earth field.
        idx = where(bmag lt bmax, cnt)
        if cnt eq 0 then begin
            ut1 = ut1+dt
            load = 1
            continue
        endif else begin
            uts = uts[idx[0]:*]
            defac = defac[idx[0]:*,*]
            dbfac = dbfac[idx[0]:*,*]
            bmag = bmag[idx[0]:*]
            emag = emag[idx[0]:*]
            posgsm = posgsm[idx[0]:*,*]
        endelse
        
        ; add a new day when less than 9 hour.
        if (uts[0] mod dt) gt 15d*3600 then begin
            ut1 = ut1+dt
            load = 1
            continue
        endif
        
        if ut1 gt max(utr0) then go = 0
    endwhile

end

tr = ['2013-07-01','2013-08-01']
;tr = ['2013-05-01','2013-05-02']
utr = time_double(tr)
find_rbsp_enormal, 'b', utr
find_rbsp_enormal, 'a', utr
end
