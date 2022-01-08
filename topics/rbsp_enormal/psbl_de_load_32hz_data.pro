;+
; load data for given id.
;-

pro psbl_de_load_32hz_data, id, trange = utr, rootdir = rootdir

    if n_elements(rootdir) eq 0 then $
        rootdir = shomedir()+'/Google Drive/works/data/psbl_de_32hz'
    
    
    datfn = rootdir+'/psbl_de_32hz_'+id+'.tplot'
    if file_test(datfn) eq 0 then begin
        if n_elements(utr) eq 2 then $
            psbl_de_save_32hz_data, id, trange = utr, rootdir = file_dirname(datfn)
        return
    endif
    
    
    ; should have the following vars.
    ; rbspx_de_[mgse,gse,fac,svy]
    ; rbspx_de_dot0_[mgse,gse,fac]
    ; rbspx_b_[mgse,gse,fac]
    ; rbspx_[b0,db]_gse
    ; rbspx_fpt_[mlat,mlon,mlt]
    ; rbspx_[ilat,mlt,lshell,n,t,beta,bx_ratio,ebangle]
    tplot_restore, filename = datfn
    
    
    tprobe = strmid(id,strlen(id)-1)
    pre0 = 'rbsp'+tprobe+'_'
    
    
    reload = 0
    if n_elements(utr) eq 2 then begin
        get_data, pre0+'de_mgse', t0
        if max(utr)-min(utr)-max(t0)+min(t0) gt 60 then reload = 1  ; required more data.
    endif
    
    if keyword_set(reload) then begin
        if n_elements(utr) eq 2 then psbl_de_save_data, id, trange = utr, rootdir = file_dirname(datfn)
        return
    endif
    
    tplot_restore, filename = datfn
    
    
    ; patches.
    
    updatefile = 0
    
    ; patch 1: dst and ae.
    get_data, 'dst', t0, dat
    if n_elements(t0) eq 1 and t0[0] eq 0 then begin
        updatefile = 1
        get_data, pre0+'de_mgse', t0
        utr = minmax(t0)
        dt0 = 86400d
        utr0 = utr-(utr mod dt0)+[0,dt0]
        omni = sread_omni(utr0)
        uts = sfmepoch(omni.epoch,'unix')
        store_data, 'dst', uts, omni.symh, limits = {ytitle:'Dst (nT)'}
        store_data, 'ae', uts, omni.ae, limits = {ytitle:'AE (nT)'}
    endif
    
    ; patch 2: rbspx_db_fac.
    get_data, pre0+'db_fac', t0, dat
    if n_elements(t0) eq 1 and t0[0] eq 0 then begin
        updatefile = 1
        get_data, pre0+'de_gse', data = de
        get_data, pre0+'db_gse', data = db
        get_data, pre0+'b0_gse', data = tmp
        bmod = sinterpol(tmp.y, tmp.x, de.x)
        db = sinterpol(db.y, db.x, de.x)
        bhat = sunitvec(bmod)
        
        p = atan(bhat[*,1],bhat[*,0])
        cosp = cos(p) & sint = bhat[*,2]
        cost = bhat[*,0]/cosp & sinp = bhat[*,1]/cost
        
        vec = db
        x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
        y =       -sinp*vec[*,0] + cosp*vec[*,1]
        z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
        
        store_data, pre0+'db_fac', data = {x:dedot0.x, y:[[x],[y],[z]]}, $
            limits = {colors:rgb, labels:'dB FAC '+['b','p','v'], ytitle:'dB!C(nT)'}
    endif
    
    ; patch 3: rbspx_[de,db,de_dot0]_fac_mat['',1,2]
    ; rbspx_pf['',_dot0]_fac_mat['',1,2]_map
    ; rbspx_map_coef
    get_data, pre0+'pf_fac_mat_map', t0, dat
    if n_elements(t0) eq 1 and t0[0] eq 0 then begin
        updatefile = 1
        
        tfilters = [20,60,250]
        filterids = ['1','2']
        nfilter = 2
        coordlabs = ['b','p','v']
        rgb = [6,4,2]
        dt = 20*60
        sclinfo = [10,dt/2,40]
        
        denames = pre0+['de_fac','de_dot0_fac']
        dbnames = pre0+['db_fac','db_fac']
        pfnames = pre0+['pf_fac','pf_dot0_fac']
        suffixs = '_mat'+filterids
        
        ; calc poynting flux.
        for i = 0, n_elements(denames)-1 do begin
            dename = denames[i]
            dbname = dbnames[i]
            pfname = pfnames[i]
            
            stplot_calc_pflux_mat, dename, dbname, pfname, $
                filter = tfilters, scaleinfo = sclinfo
                
            vars = [dename,dbname,pfname]+'_mat'
            options, tvar, 'constant', 0
            
            ; map poynting flux.
            vars = pfname+[suffixs,'_mat']
            for j = 0, n_elements(vars)-1 do begin
                tvar = vars[j]
                get_data, pre0+'b0_gse', t0, dat
                store_data, pre0+'b0', t0, sqrt(total(dat^2,2))
                smap2iono, tvar, pre0+'fpt_mlat', b = pre0+'b0', $
                    newname = tvar+'_map', coef = coef
                tvar = vars[j]+'_map'
                options, tvar, 'constant', 0
                options, tvar, 'colors', rgb
            endfor
            store_data, pre0+'_map_coef', t0, coef
        endfor
    endif
    
    
    
    if updatefile then begin
        vars = pre0+['de_mgse','de_gse','de_fac','de_svy', $
            'de_dot0_mgse','de_dot0_gse','de_dot0_fac', $
            'b_mgse','b_gse','db_gse','db_fac','b0_gse', $
            'bx_ratio','ebangle','ilat','mlt','lshell', $
            'fpt_mlat','fpt_mlon','fpt_mlt', $
            'n','t','beta', $
            'de_fac_mat'+['','1','2'], $
            'de_dot0_fac_mat'+['','1','2'], $
            'db_fac_mat'+['','1','2'], $
            'pf_fac_mat'+['','1','2']+'_map', $
            'pf_dot0_fac_mat'+['','1','2']+'_map', 'map_coef']
        vars = [vars,'dst','ae']
        tplot_save, vars, filename = datfn
    endif
    

end


psbl_de_load_data, '2013_0126_2114_b'
end
