;+
; loads RBSP and ASI data.
; RBSP E/B at 16 Hz, in situ and mapped position, HOPE and MAGEIS spectrogram
; ASI keogram.
; RBSP-A DF 04:54:32UT, MLT 22.6+/-0.1, MLat 64.3+/-1.8.
; RBSP-B DF 04:54:59UT, MLT 22.2+/-0.1, MLat 64.4+/-1.8.
;-
pro _2013_0501_load_data

;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    secofday = 86400d

    fac = ['b','w','n'] ; para, west, north.
    xyz = ['x','y','z']
    rgb = [6,4,2]

;---settings.
    utr0 = time_double(['2013-05-01/07:30','2013-05-01/07:45'])  ; long time range.
    utr1 = utr0+[-1,1]*300   ; 5 min pad time.
    utr2 = time_double(['2013-05-01/07:30','2013-05-01/07:45'])  ; short time range.
    utr3 = time_double(['2013-05-01/07:30','2013-05-01/07:45'])  ; ewo time range.
    utr4 = time_double(['2013-05-01/07:30','2013-05-01/07:45'])  ; keo time range.
    dfuts = time_double(['2013-05-01/07:38'])  ; DF time.


    rootdir = shomedir()+'/Google Drive/works/works/rbsp_psbl1'
    datdir = rootdir+'/data'
    if file_test(datdir,/directory) eq 0 then file_mkdir, datdir
    datfn = datdir+'/2013_0501_all_data.tplot'
    if file_test(datfn) eq 1 then tplot_restore, filename=datfn

    ; control loading behavior.
    reload_pos = 0
    reload_map = 0
    reload_info = 1
    reload_eb = 0
    reload_pf = 0
    reload_asi = 0
    reload_esa = 0


    ; for mapping
    models = ['t89','t96','t01','t04s']  ; the B0 model ~ measured. tried t04s, t96, t89, B0 are different.
    modelcolors = [6,4,3,2]
    nmodel = n_elements(models)
    model0 = 't01'

    ; for rbsp.
    probes = ['b']
    nprobe = n_elements(probes)
    pres = 'rb'+probes+'_'
    dr0 = 1d/16     ; the common datarate for E/B field.
    spinrate = 12d  ; sec.
    ; for e calc.
    blen1s = [100d,100,12]  ; m, twice of boom lengths.
    ; for edot0.
    bmin = 5d   ; nT.
    rmin = 0.25 ; percent.
    ; for hope. 
    ionminen = 30
    eleminen = 200
    ; for pflux calc.
    scaleinfo = {s0:0.25d, s1:1200d, dj:1d/8, ns:0d}    ; pflux calc.

    ; for asi.
    sites = ['atha']
    hgtidx = 1      ; for 110 km.
    minelev = 5    ; deg.
    ; mlats for ewogram.
    ewo_dlat = 0.5d     ; width in mlat.
    ewo_lats = string(smkarthm(63.5,64.5,ewo_dlat,'dx'),format='(F4.1)')
    ewo_res = 0.05d     ; resolution in mlt.
    mltrng = [-1.7,0]
    mltrng = mltrng-(mltrng mod ewo_res)+[0,ewo_res]
    ; mlts for keogram.
    keo_dmlt = 0.2d     ; width in mlt.
    ;keo_mlts = string([22.1,22.5]-24,format='(F4.1)') use einfo to calc.
    keo_res = 0.1d      ; resolution in mlat.
    latrng = [60d,68]
    keovals = smkarthm(latrng[0], latrng[1], keo_res)

    ; for mapping.
    dir = -1    ; -1 for northern hemisphere.



;---One-time settings.
    allvars = []

    ; event info vars.
    infovars = ['event_info']

    ; position vars.
    tvar = ['r_gsm','mlt','mlat','lshell','dis','q_uvw2gsm']
    posvars = []
    foreach pre0, pres do posvars = [posvars,pre0+tvar]

    ; mapping vars.
    tvar = ['fpt_'+['mlt','mlat','mlon'],'bmod_gsm_'+models,'map_coef0']
    mapvars = models+'_par'
    foreach pre0, pres do mapvars = [mapvars,pre0+tvar]

    ; E/B field vars.
    tvar = ['b_gsm','b0_gsm','db_gsm','bmag','e_gsm','edot0_gsm']
    ebvars = []
    foreach pre0, pres do ebvars = [ebvars,pre0+tvar]

    ; pflux vars.
    tvar = [['de','db']+'_fac',['de','db']+'_mor_spec','pf_fac','pf_para_mor_spec','map_coef1']
    pfvars = []
    foreach pre0, pres do pfvars = [pfvars,pre0+tvar]

    ; particle vars.
    tvar = ['e_en_mageis','e_en','e_pa','h_en','h_pa','o_en','o_pa', $
        'ni','n','t','n_efw','n_combine','vsc']
    esavars = []
    foreach pre0, pres do esavars = [esavars,pre0+tvar]

    ; asi vars.
    tvar = ['asi','ewo_'+ewo_lats]
    asivars = []
    foreach tsite, sites do asivars = [asivars, tsite+'_'+tvar]
    foreach pre0, pres do foreach tsite, sites do asivars = [asivars, pre0+'keo_'+tsite]


    allvars = [infovars,posvars,mapvars,ebvars,pfvars,esavars,asivars]



    ; plot settings.
    device, decomposed = 0
    loadct2, 43
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'xticklen', -0.02
    tplot_options, 'yticklen', -0.01
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 6
    tplot_options, 'ynozero', 1



    
;---Load position.
; load in-situ position and spice kernel uvw pointing dir.
    vars = posvars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_pos) then load = 1

    if load then begin
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            spice = sread_rbsp_spice_product(utr0, probe=tprobe)

            uts = spice.ut_pos
            store_data, pre0+'r_gsm', uts, spice.pos_gsm, limits={ytitle:'(Re)',colors:rgb, labels:xyz}
            store_data, pre0+'mlt', uts, spice.mlt
            store_data, pre0+'mlat', uts, spice.mlat*deg
            store_data, pre0+'lshell', uts, spice.lshell
            store_data, pre0+'dis', uts, snorm(spice.pos_gsm)

            uts = spice.ut_cotran
            quvw2gsm = spice.q_uvw2gsm
            store_data, pre0+'q_uvw2gsm', uts, quvw2gsm
        endforeach
        tplot_save, allvars, filename=datfn
    endif


;---Load mapped position.
; load geopack pars for models, trace fpt, calc mapping coef based on model B.
    vars = mapvars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_map) then load = 1

    if load then begin
        ; load geopack par's. [models_par].
        foreach tmodel, models do sgeopack_par, utr0, tmodel

        foreach pre0, pres do begin
            get_data, pre0+'r_gsm', tuts, posgsm
            tets = stoepoch(tuts, 'unix')
            tnrec = n_elements(tuts)
            hfmlats = dblarr(tnrec,nmodel)
            hfmlons = dblarr(tnrec,nmodel)
            hfmlts = dblarr(tnrec,nmodel)
            bmods = dblarr(tnrec,3)
            mapcoefs = dblarr(tnrec,nmodel)
            
            
            for j=0, nmodel-1 do begin
                tmodel = models[j]
                suf0 = '_'+tmodel
                t89 = (tmodel eq 't89')? 1: 0
                t96 = (tmodel eq 't96')? 1: 0
                t01 = (tmodel eq 't01')? 1: 0
                t04s = (tmodel eq 't04s')? 1: 0
                storm = (tmodel eq 't04s')? 1: 0
                
                get_data, tmodel+'_par', tmp, pars
                pars = sinterpol(pars, tmp, tuts)
                tnrec = n_elements(tuts)
                
                for i=0, tnrec-1 do begin
                    tet = tets[i]
                    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
                    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
                    
                    xp = posgsm[i,0]
                    yp = posgsm[i,1]
                    zp = posgsm[i,2]
                    ; map to the ionosphere.
                    ; [xyz]f is the footpoint position.
                    par = reform(pars[i,*])
                    geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, $
                        /refine, /ionosphere, $
                        t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                    ; convert footpoint from gsm to mag.
                    geopack_conv_coord, xf,yf,zf, /from_gsm, $
                        txf,tyf,tzf, /to_mag
                    hfmlats[i,j] = asin(tzf/r0)*deg
                    hfmlons[i,j] = atan(tyf,txf)*deg
                    hfmlts[i,j] = slon2lt(hfmlons[i,j], tet, /mag, /deg)/15

                    ; get the model B field at footprint and s/c pos.
                    geopack_igrf_gsm, xf,yf,zf, bxf,byf,bzf

                    ; internal B field.
                    geopack_igrf_gsm, xp,yp,zp, bxp,byp,bzp
                    ; external B field.
                    case tmodel of
                        't89': geopack_t89, par, xp,yp,zp, tbx,tby,tbz
                        't96': geopack_t96, par, xp,yp,zp, tbx,tby,tbz
                        't01': geopack_t01, par, xp,yp,zp, tbx,tby,tbz
                        't01s': geopack_t01, par, xp,yp,zp, tbx,tby,tbz /storm
                        't04s': geopack_ts04, par, xp,yp,zp, tbx,tby,tbz
                        else: begin tbx = 0 & tby = 0 & tbz = 0 & end
                    endcase
                    bxp+= tbx & byp+= tby & bzp+= tbz
                    bmods[i,*] = [bxp,byp,bzp]
                    mapcoefs[i,j] = snorm([bxf,byf,bzf])/snorm([bxp,byp,bzp])
                endfor
                store_data, pre0+'bmod_gsm'+suf0, tuts, bmods, limits=$
                    {ytitle:'(nT)', labels:'Bmod GSM'+xyz, colors:rgb}
            endfor
            store_data, pre0+'fpt_mlat', tuts, hfmlats, limits={ytitle:'(deg)', colors:modelcolors, labels:models}
            store_data, pre0+'fpt_mlon', tuts, hfmlons, limits={ytitle:'(deg)', colors:modelcolors, labels:models}
            ; fix mlt, use the same conversion for aurora.
            midns = -15.5179*(tuts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
            hfmlts = dblarr(tnrec,nmodel)
            for i=0, nmodel-1 do hfmlts[*,i] = (hfmlons[*,i]-midns)/15
            store_data, pre0+'fpt_mlt', tuts, hfmlts, limits={ytitle:'(hr)', colors:modelcolors, labels:models}
            store_data, pre0+'map_coef0', tuts, mapcoefs, limits=$
                {ytitle:'', labels:models, colors:modelcolors}
        endforeach
        tplot_save, allvars, filename=datfn
    endif


;---event info.
    load = 0
    foreach tvar, infovars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_info) then load = 1
    
    if load then begin
        info0 = {dfut:0d, fmlat:dblarr(nmodel), fmlt:dblarr(nmodel)}
        tid = time_string(utr0[0],tformat='YYYY_MMDD')
        einfo = {id: tid, model: models}

        foreach tprobe, probes, i do begin
            pre0 = pres[i]
            tinfo = info0
            
            tdfut = dfuts[i]
            get_data, pre0+'fpt_mlt', uts, fmlt
            tfmlts = sinterpol(fmlt, uts, tdfut)
            get_data, pre0+'fpt_mlat', uts, fmlat
            tfmlats = sinterpol(fmlat, uts, tdfut)
            
            tinfo.fmlat = tfmlats
            tinfo.fmlt = tfmlts
            tinfo.dfut = tdfut
            einfo = create_struct(pre0,tinfo, einfo)
            
            print, 'RB-'+strupcase(tprobe)
            print, 'Dipolarization time:  ', time_string(dfuts[i])
            print, 'Footpoint MLT (hr):   ', sgnum2str((mean(tfmlts)+24) mod 24,ndec=1), '+/-', sgnum2str(stddev(tfmlts),ndec=1)
            print, 'Footpoint MLat (deg): ', sgnum2str(mean(tfmlats),ndec=1), '+/-', sgnum2str(stddev(tfmlats),ndec=1)
        endforeach
        
        store_data, 'event_info', 0, einfo
        
        tplot_save, allvars, filename=datfn
    endif




;---Load E/B field.
; load E/B field, separate B0 and dB from B.
    vars = ebvars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_eb) then load = 1

    if load then begin
        ; uniform time for E/B.
        uts = smkarthm(utr0[0], utr0[1], dr0, 'dx')

        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            
            ; load Q.
            get_data, pre0+'q_uvw2gsm', tuts, quvw2gsm
            quvw2gsm = qslerp(quvw2gsm, tuts, uts)
            muvw2gsm = transpose(qtom(quvw2gsm))

            ; load B.
            emfisis = sread_rbsp_emfisis_l3(utr0, type='hires', probe=tprobe, coord='gsm')
            tuts = sfmepoch(emfisis.epoch,'unix')
            bgsm = sinterpol(emfisis.mag,tuts,uts)
            store_data, pre0+'b_gsm', uts, bgsm, $
                limits = {ytitle:'(nT)', labels:'GSM B'+xyz, colors:rgb}
            b0gsm = bgsm
            for i=0,2 do b0gsm[*,i] = scalcbg(bgsm[*,i])
            store_data, pre0+'b0_gsm', uts, b0gsm, $
                limits = {ytitle:'(nT)', labels:'GSM B0'+xyz, colors:rgb}
            store_data, pre0+'db_gsm', uts, bgsm-b0gsm, $
                limits = {ytitle:'(nT)', labels:'GSM dB'+xyz, colors:rgb}
            bmag = snorm(bgsm)
            bmag = smooth(bmag, spinrate/dr0, /edge_truncate)
            store_data, pre0+'bmag', uts, bmag, limits = {ytitle:'(nT)', labels:'|B|'}

            ; efield.
            rbsp_preprocess_efield, utr0, probe=tprobe
            get_data, 'rbsp'+tprobe+'_euvw', tuts, dat
            dat = sinterpol(dat, tuts, uts)
            store_data, pre0+'e_uvw', uts, dat
            store_data, 'rbsp'+tprobe+'*', /delete
           
            ; calc e gsm.
            get_data, pre0+'e_uvw', uts, euvw
            euvw[*,2] = 0
            ex = euvw[*,0]*muvw2gsm[0,0,*] + euvw[*,1]*muvw2gsm[1,0,*] + euvw[*,2]*muvw2gsm[2,0,*]
            ey = euvw[*,0]*muvw2gsm[0,1,*] + euvw[*,1]*muvw2gsm[1,1,*] + euvw[*,2]*muvw2gsm[2,1,*]
            ez = euvw[*,0]*muvw2gsm[0,2,*] + euvw[*,1]*muvw2gsm[1,2,*] + euvw[*,2]*muvw2gsm[2,2,*]
            tvar = pre0+'e_gsm'
            store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
                {ytitle:'(mV/m)', colors:rgb, labels:'GSM E'+xyz}
                
            ; calc e dot0 gsm.
            get_data, pre0+'b0_gsm', uts, b0gsm
            bx = b0gsm[*,0] & by = b0gsm[*,1] & bz = b0gsm[*,2]
            bu = bx*muvw2gsm[0,0,*] + by*muvw2gsm[0,1,*] + bz*muvw2gsm[0,2,*]
            bv = bx*muvw2gsm[1,0,*] + by*muvw2gsm[1,1,*] + bz*muvw2gsm[1,2,*]
            bw = bx*muvw2gsm[2,0,*] + by*muvw2gsm[2,1,*] + bz*muvw2gsm[2,2,*]
                        
            buvw = [[bu],[bv],[bw]]
            bmag = snorm(buvw)
            euvw[*,2] = -(euvw[*,0]*buvw[*,0]+euvw[*,1]*buvw[*,1])/buvw[*,2]
            idx = where(bmag le bmin or abs(buvw[*,2])/bmag le rmin, cnt)
            if cnt ne 0 then euvw[idx,2] = !values.d_nan
            
            ex = euvw[*,0]*muvw2gsm[0,0,*] + euvw[*,1]*muvw2gsm[1,0,*] + euvw[*,2]*muvw2gsm[2,0,*]
            ey = euvw[*,0]*muvw2gsm[0,1,*] + euvw[*,1]*muvw2gsm[1,1,*] + euvw[*,2]*muvw2gsm[2,1,*]
            ez = euvw[*,0]*muvw2gsm[0,2,*] + euvw[*,1]*muvw2gsm[1,2,*] + euvw[*,2]*muvw2gsm[2,2,*]

            tvar = pre0+'edot0_gsm'
            store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
                {ytitle:'(mV/m)', colors:rgb, labels:'GSM Edot0'+xyz}
        endforeach
        tplot_save, allvars, filename=datfn
    endif



;---Load pflux.
; convert E/B in GSM to FAC, calc pflux in FAC, modify mapping coef.
    vars = pfvars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_pf) then load = 1

    if load then begin
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            
            ; convert E/B field to fac.
            get_data, pre0+'b0_gsm', uts, b0gsm
            get_data, pre0+'r_gsm', tuts, rgsm
            rgsm = sinterpol(rgsm, tuts, uts)
            rhat = sunitvec(rgsm)
            bhat = sunitvec(b0gsm)
            phat = sunitvec(scross(rhat,bhat))
            vhat = scross(bhat,phat)

            tvar = keyword_set(use_edot0)? pre0+'edot0_gsm': pre0+'e_gsm'
            get_data, tvar, uts, egsm
            defac = [[sdot(egsm,bhat)],[sdot(egsm,phat)],[sdot(egsm,vhat)]]
            store_data, pre0+'de_fac', uts, defac, limits = {ytitle:'(mV/m)', $
                labels:'FAC dE'+fac, colors:rgb}
            tvar = pre0+'db_gsm'
            get_data, tvar, uts, dbgsm
            dbfac = [[sdot(dbgsm,bhat)],[sdot(dbgsm,phat)],[sdot(dbgsm,vhat)]]
            store_data, pre0+'db_fac', uts, dbfac, limits = {ytitle:'(nT)', $
                labels:'FAC dB'+fac, colors:rgb}


            tvar = pre0+'de_mor_spec'
            get_data, pre0+'de_fac', uts, vec
            store_data, pre0+'demag', uts, snorm(vec)
            stplot_mor, pre0+'demag', scaleinfo=scaleinfo, newname = tvar
            options, tvar, 'zrange', [0,1]
            options, tvar, 'ztitle', '|dE|!U2!N (mV/m)'
            options, tvar, 'ytitle', 'Period!C(sec)'
            
            tvar = pre0+'db_mor_spec'
            get_data, pre0+'db_fac', uts, vec
            store_data, pre0+'dbmag', uts, snorm(vec)
            stplot_mor, pre0+'dbmag', scaleinfo = scaleinfo, newname = tvar
            options, tvar, 'zrange', [0,1]
            options, tvar, 'ztitle', '|dB|!U2!N (nT)'
            options, tvar, 'ytitle', 'Period!C(sec)'
            
            ; instataneous pflux.
            stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', $
                scaleinfo=scaleinfo
            
            ; save the parallel pflux spectrogram.
            stplot_renew, pre0+'pf_fac_mor_spec_1', newname=pre0+'pf_para_mor_spec'

            ; map coef updated with real B field.
            get_data, pre0+'map_coef0', tuts, mcoef, limits=lim
            get_data, pre0+'bmag', uts, bmag
            bmag = sinterpol(bmag, uts, tuts)
            foreach tmodel, models, j do begin
                get_data, pre0+'bmod_gsm_'+tmodel, tuts, bmodgsm
                bmodmag = snorm(bmodgsm)
                mcoef[*,j] *= (bmodmag/bmag)
            endforeach
            store_data, pre0+'map_coef1', tuts, mcoef, limits=lim
        endforeach
        
        ; update data file.
        tplot_save, allvars, filename = datfn
    endif



;---Load asi.
    vars = asivars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_asi) then load = 1

    if load then begin
        foreach tsite, sites do begin
            pre0 = tsite+'_'
            
            ; load site info.
            asc = sread_thg_asc(0, tsite, type='asf', vars=['mlat','mlon','elev'])
            mlats = reform(asc.(0).mlat[hgtidx,*,*])
            mlons = reform(asc.(0).mlon[hgtidx,*,*])
            elevs = asc.(0).elev
            edge = where(elevs lt minelev or ~finite(elevs))
            imgsz = size(elevs,/dimensions)
            midnrng = -15.5179*(utr3*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.

            ewovals = smkarthm(mltrng[0], mltrng[1], ewo_res, 'dx')
            newoval = n_elements(ewovals)

            tmp = minmax(mlats)
            nkeoval = n_elements(keovals)

            ; converts corner position to center position.
            mlats = (mlats[1:*,1:*]+mlats[1:*,0:imgsz[1]-1]+mlats[0:imgsz[0]-1,1:*]+mlats[0:imgsz[0]-1,0:imgsz[1]-1])*0.25
            mlons = (mlons[1:*,1:*]+mlons[1:*,0:imgsz[1]-1]+mlons[0:imgsz[0]-1,1:*]+mlons[0:imgsz[0]-1,0:imgsz[1]-1])*0.25

            ; load auroral images.
            asf = sread_thg_asi(utr3, tsite, type='asf')
            uts = asf.utsec
            nrec = n_elements(uts)
            midns = -15.5179*(uts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
            
            ; image preprocess.
            imgs = fltarr([imgsz,nrec])
            for i=0, nrec-1 do begin
                timg = float(reform(asf.img[i,*,*]))
                ; remove edge.
                timg[edge] = mean(timg[0:10,0:10], /nan)
                ; crude corner average subtraction. use thg_asf_site_offset?
                timg = (timg-timg[0]) > 0
                ; scale luminosity, adopted from thm_asi_merge_mosaic.
                timg *= 64d/(median(timg) > 1)
                imgs[*,*,i] = timg
            endfor
            store_data, pre0+'asi', uts, imgs, {mlats:mlats, mlons:mlons, elevs:elevs, midns:midns}
            
            ; ewogram at given mlats.
            foreach tewo_lat, ewo_lats do begin
                tlatrng = float(tewo_lat)+[0,1]*ewo_dlat
                suf0 = '_'+tewo_lat
                ewos = fltarr(nrec,newoval)
                for i=0, nrec-1 do begin
                    timg = imgs[*,*,i]
                    tmlts = (mlons-midns[i])/15
                    ; map to ewogram.
                    for j=0, newoval-1 do begin
                        idx = where(mlats ge tlatrng[0] and mlats le tlatrng[1] $
                            and tmlts ge ewovals[j]-ewo_res*0.5 and tmlts lt ewovals[j]+ewo_res*0.5, cnt, complement=idx2)
                        if cnt ne 0 then ewos[i,j] = mean(timg[idx],/nan)
                    endfor
                endfor
                tvar = pre0+'ewo'+suf0
                store_data, tvar, uts, ewos, ewovals, limits= $
                    {ytitle:'MLT (hr)', spec:1, yrange:mltrng, no_interp:1, ystyle:1}
            endforeach


            ; load auroral images for long time range.
            asf = sread_thg_asi(utr4, tsite, type='asf')
            uts = asf.utsec
            nrec = n_elements(uts)
            midns = -15.5179*(uts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
            
            ; image preprocess.
            imgs = fltarr([imgsz,nrec])
            for i=0, nrec-1 do begin
                timg = float(reform(asf.img[i,*,*]))
                ; remove edge.
                timg[edge] = mean(timg[0:10,0:10], /nan)
                ; crude corner average subtraction. use thg_asf_site_offset?
                timg = (timg-timg[0]) > 0
                ; scale luminosity, adopted from thm_asi_merge_mosaic.
                timg *= 64d/(median(timg) > 1)
                imgs[*,*,i] = timg
            endfor

            get_data, 'event_info', 0, einfo
            vnames = strlowcase(tag_names(einfo))

            ; keogram at given mlts.
            foreach pre0, pres, k do begin
                tidx = where(vnames eq pre0)
                tkeo_mlt = einfo.(tidx).fmlt
                tmltrng = float(tkeo_mlt)+[-1,1]*0.5*keo_dmlt
                suf0 = '_'+tsite
                keos = fltarr(nrec,nkeoval)
                for i=0, nrec-1 do begin
                    timg = imgs[*,*,i]
                    tmlts = (mlons-midns[i])/15
                    ; map to keogram.
                    for j=0, nkeoval-1 do begin
                        idx = where(tmlts ge tmltrng[0] and tmlts le tmltrng[1] $
                            and mlats ge keovals[j]-keo_res*0.5 and mlats le keovals[j]+keo_res*0.5, cnt, complement=idx2)
                        if cnt ne 0 then keos[i,j] = mean(timg[idx],/nan)
                    endfor
                endfor
                tvar = pre0+'keo'+suf0
                store_data, tvar, uts, keos, keovals, limits= $
                    {ytitle:'MLat (deg)', spec:1, yrange:latrng, no_interp:1, ystyle:1}
            endforeach
        endforeach
        tplot_save, allvars, filename=datfn
    endif


;---load particle data.
    vars = esavars
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_esa) then load = 1

    if load then begin
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            
            mageis = sread_rbsp_mageis_l3(utr0, probes = tprobe)
            uts = sfmepoch(mageis.epoch,'unix')
            fedu = total(mageis.fedu,3,/nan)
            store_data, pre0+'e_en_mageis', uts, fedu, mageis.fedu_energy, limits = $
                {spec:1, no_interp:1, ylog:1, zlog:1, yrange:[4e1,4e3], ystyle:1, $
                ytitle:'MAGEIS!CElectron!CEenergy!C(keV)', $
                zticks:3, zrange:[1e-2,1e5], $
                ztitle:'(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'}
            
            vars = ['PITCH_ANGLE',$
                'Epoch_Ele','HOPE_ENERGY_Ele',$
                'Epoch_Ion','HOPE_ENERGY_Ion',$
                'FEDU','FPDU','FODU']
                
            hopel3 = sread_rbsp_hope_l3(utr0, probes = tprobe, vars = vars)
            pas = hopel3.pitch_angle

            ; ion spec.
            uts = sfmepoch(hopel3.epoch_ion,'unix')
            nrec = n_elements(uts)
            ens = hopel3.hope_energy_ion

            dat = hopel3.fodu
            type = 'o'
            zr = [1e4,1e7]
            idx = where(finite(dat,/nan),cnt)
            if cnt ne 0 then dat[idx] = 0
            for i = 0, nrec-1 do begin
                idx = where(ens[i,*] le ionminen, cnt)
                if cnt ne 0 then dat[i,idx,*] = 0
            endfor
            
            paspec = total(dat, 2)
            idx = where(paspec eq 0)
            paspec[idx] = 7
            
            tvar = pre0+type+'_pa'
            store_data, tvar, uts, paspec, pas, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
                yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
                ytitle:'Pitch Angle (deg)'}
                
            enspec = total(dat, 3)
            idx = where(enspec eq 0)
            enspec[idx] = 7
            
            tvar = pre0+type+'_en'
            store_data, tvar, uts, enspec, ens, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr, $
                yrange:[3,3e4], ystyle:1, ylog:1, $
                ytitle:'Energy (eV)'}
            
            
            dat = hopel3.fpdu
            type = 'h'
            zr = [1e4,1e7]
            idx = where(finite(dat,/nan),cnt)
            if cnt ne 0 then dat[idx] = 0
            for i = 0, nrec-1 do begin
                idx = where(ens[i,*] le ionminen, cnt)
                if cnt ne 0 then dat[i,idx,*] = 0
            endfor
            
            paspec = total(dat, 2)
            idx = where(paspec eq 0)
            paspec[idx] = 7
            
            tvar = pre0+type+'_pa'
            store_data, tvar, uts, paspec, pas, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
                yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
                ytitle:'Pitch Angle (deg)'}
                
            enspec = total(dat, 3)
            idx = where(enspec eq 0)
            enspec[idx] = 7
            
            tvar = pre0+type+'_en'
            store_data, tvar, uts, enspec, ens, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr, $
                yrange:[30,3e4], ystyle:1, ylog:1, $
                ytitle:'Energy (eV)'}
            
            
            ; electron.
            uts = sfmepoch(hopel3.epoch_ele,'unix')
            nrec = n_elements(uts)
            ens = hopel3.hope_energy_ele
            
            dat = hopel3.fedu
            type = 'e'
            zr = [1e7,1e10]
            idx = where(finite(dat,/nan),cnt)
            if cnt ne 0 then dat[idx] = 0
            for i = 0, nrec-1 do begin
                idx = where(ens[i,*] le ionminen, cnt)
                if cnt ne 0 then dat[i,idx,*] = 0
            endfor
            
            paspec = total(dat, 2)
            idx = where(paspec eq 0)
            paspec[idx] = 7
            
            tvar = pre0+type+'_pa'
            store_data, tvar, uts, paspec, pas, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
                yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
                ytitle:'Pitch Angle (deg)'}
                
            enspec = total(dat, 3)
            idx = where(enspec eq 0)
            enspec[idx] = 7
            
            tvar = pre0+type+'_en'
            store_data, tvar, uts, enspec, ens, limits = $
                {spec:1, no_interp:1, zlog:1, zrange:zr, $
                yrange:[200,2e4], ylog:1, ystyle:1, $
                ytitle:'Energy (eV)'}
            
            vars = pre0+['h_en','h_pa','e_en','e_pa','o_en','o_pa']
            options, vars, 'spec', 1
            options, vars, 'no_interp', 1
            options, vars, 'zlog', 1
            options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'
            
            vars = pre0+['h_pa','e_pa','o_pa']
            options, vars, 'yrange', [0,180]
            options, vars, 'ystyle', 1
            
            tvar = pre0+'e_en'
            options, tvar, 'ytitle', 'e- energy!C(eV)'
            tvar = pre0+'h_en'
            options, tvar, 'ytitle', 'H+ energy!C(eV)'
            
            tvar = pre0+'e_pa'
            options, tvar, 'ytitle', 'e- pitch!C(deg)'
            tvar = pre0+'h_pa'
            options, tvar, 'ytitle', 'H+ pitch!C(deg)'
            tvar = pre0+'o_pa'
            options, tvar, 'ytitle', 'O+ pitch!C(deg)'
            
            
            hopemom = sread_rbsp_hope_l3(utr0, probes=tprobe, type='mom', $
                vars = ['Epoch_Ion','Dens_o_30','Dens_p_30'])
            if size(hopemom,/type) ne 8 then return
            uts = sfmepoch(hopemom.epoch_ion,'unix')
            store_data, pre0+'ni', uts, [[hopemom.dens_p_30],[hopemom.dens_o_30]], $
                limits = {ytitle:'n!Di!N!C(cm!U-3!N)', ylog:1, constant:1, labels:$
            ['n!DH!N>30eV','n!DO!N>30eV'], colors:[0,6]}
            
            
            hopemom = sread_rbsp_hope_l3(utr0, probes=tprobe, type='mom')
            if size(hopemom,/type) ne 8 then return
            uts = sfmepoch(hopemom.epoch_ele,'unix')
            store_data, pre0+'n', uts, hopemom.dens_e_200, $
                limits = {ytitle:'n!De!N!C(cm!U-3!N)', ylog:1, constant:1, labels:'n!De!N>200eV'}
            store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
                limits = {ytitle:'T!De!N!C(eV)', ylog:1, colors:[6,0], labels:['Tperp','Tpara'], constant:1000}
            tvar = pre0+'t'
            get_data, tvar, uts, dat
            idx = where(dat eq 1e20, cnt)
            if cnt ne 0 then dat[idx] = !values.d_nan
            store_data, tvar, uts, dat
            

            ; fit Vsc to density.
            get_data, pre0+'e_gsm', uts
            dat = sread_rbsp_efw_l2(utr0, probes = tprobe, type = 'vsvy')
            tuts = sfmepoch(dat.epoch, 'unix')
            vsvy = dat.vsvy
            vsc = mean(vsvy[*,0:1], dimension=2)
            vsc = sinterpol(vsc, tuts, uts)
            store_data, pre0+'vsc', uts, vsc

            get_data, pre0+'n', tuts, hopen
            idx = where(finite(hopen))
            tuts = tuts[idx]
            hopen = hopen[idx]
            hopen = sinterpol(hopen, tuts, uts, /nan)
            
            ; hopen = a*exp(vsc*b)
            ; alog(hopen) = alog(a)+vsc*b
            tx0 = vsc
            ty0 = alog(hopen)
            
            a = regress(tx0, ty0, const = b)
            ty1 = a[0]*vsc+b[0]
            
            ybar = mean(ty0)
            sstot = total((ty0-ybar)^2)
            ssres = total((ty0-ty1)^2)
            
            r2 = 1-ssres/sstot
            print, a, b, r2
            
            efwn = exp(ty1)
            store_data, pre0+'n_efw', uts, efwn, limits = $
                {ytitle:'Density!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
                labels:'n!De,EFW Vsc'}
            store_data, pre0+'n_combine', uts, [[efwn],[hopen]], limits = $
                {ytitle:'RBSP-'+strupcase(tprobe)+'!CDensity!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
                labels:['n!De,EFW Vsc!N', 'n!De,HOPE > 200 eV'], colors:[0,6]}
        endforeach
        tplot_save, allvars, filename = datfn
    endif



;---further settings for above vars.
    if load then begin
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            get_data, pre0+'map_coef', tuts, tdat
            tdat = reform(tdat[*,where(models eq model0)])
            get_data, pre0+'pf_fac', uts, pf
            tdat = sinterpol(tdat, tuts, uts)
            for i = 0, 2 do pf[*,i] *= tdat
            store_data, pre0+'pf_fac_map', uts, pf, limits = $
                {ytitle:'S@100 km!C(mW/m!U2!N)', colors:rgb, labels:['b','p','v']}
        

            tvar = pre0+'e_pa'
            options, tvar, 'ytitle', 'HOPE!CElectron!CPitch Angle!C(deg)'
            
            tvar = pre0+'e_en'
            options, tvar, 'ytitle', 'HOPE!CElectron!CEnergy!C(eV)'
            options, tvar, 'yticks', 2
            options, tvar, 'ytickv', [1e2,1e3,1e4]
            options, tvar, 'yminor', 10
            options, tvar, 'ystyle', 1
            options, tvar, 'yrange', [200,2e4]
            options, tvar, 'zrange', [1e4,1e10]
            options, tvar, 'zticks', 3
            
            tvar = pre0+'h_en'
            options, tvar, 'ytitle', 'HOPE!CProton!CEnergy!C(eV)'
            options, tvar, 'yticks', 2
            options, tvar, 'ytickv', [1e2,1e3,1e4]
            options, tvar, 'yminor', 10
            options, tvar, 'ystyle', 1
            options, tvar, 'yrange', [40,4e4]
            options, tvar, 'zrange', [1e4,1e7]
            options, tvar, 'zticks', 3
            
            tvar = pre0+'o_en'
            options, tvar, 'ytitle', 'HOPE!COxygen!CEnergy!C(eV)'
            options, tvar, 'yticks', 2
            options, tvar, 'ytickv', [1e2,1e3,1e4]
            options, tvar, 'yminor', 10
            options, tvar, 'ystyle', 1
            options, tvar, 'yrange', [40,4e4]
            options, tvar, 'zrange', [1e4,1e7]
            options, tvar, 'zticks', 3
            
            tvar = pre0+'vsc'
            options, tvar, 'yticks', 2
            
            tvar = pre0+'b_gsm'
            options, tvar, 'ystyle', 1
            options, tvar, 'yrange', [-150,300]
            options, tvar, 'yticks', 3
            options, tvar, 'yminor', 5
            
            tvar = pre0+'pf_fac_map'
            perp = '!9'+string(94b)+'!X'
            options, tvar, 'ystyle', 1
            options, tvar, 'yrange', [-90,180]
            options, tvar, 'yticks', 3
            options, tvar, 'yminor', 5
            options, tvar, 'ytitle', 'S @100 km!C(mW/m!U2!N)'
            options, tvar, 'labels', 'S!D'+['||',perp+',West',perp+',North']
        endforeach
    endif



end

