;+
; on daily bases, preprocess efield for certain probe.
; subtract background E field near earth (<3 Re)
; remove bad E field when Eu and Ev do not match and
; when they are asymetric around 0.
;-


pro themis_preprocess_efield, ut0, probes = probes, write_log=write_log


;---Settings.
    if n_elements(probes) eq 0 then probes = ['a','b','c','d','e']
    rootdir = shomedir()+'/thm_efield_preprocess'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    logfn = rootdir+'/thm_efield_preprocess.log'
    
    dt0 = 86400d
    utr0 = ut0-(ut0 mod dt0)+[0,dt0]
    dr0 = 1d/8      ; data rate, Hz.
    spinrate = 3    ; sec.
    mintlen = 3600  ; at least 1 hr of data.
    dis0 = 7d   ; Re, show B field beyond dis0.
    
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re

    noplots = 1
    
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    uvw = ['u','v','w']
    
    if ~noplots then begin
        device, decomposed = 0
        loadct2, 43
        tplot_options, 'version', 2
        tplot_options, 'num_lab_min', 10
        tplot_options, 'labflag', -1
        tplot_options, 'xticklen', -0.03
        tplot_options, 'yticklen', -0.005
        tplot_options, 'constant', 0
    endif
    
    
    date = time_string(utr0[0],tformat='YYYY_MMDD')

    
    foreach tprobe, probes do begin
        
        pre0 = 'th'+tprobe+'_'

    ;---Load data.
        ; thx_[pos_gsm,mlt,dis].
        posvar = ['Epoch','XYZ_GSM','DNEUTS','RADIUS','SM_LCT_T']
        dat = sread_thm_orbit(utr0, probes=tprobe, vars=posvar)
        if size(dat,/type) ne 8 then begin
            if keyword_set(write_log) then begin
                openw, lun, logfn, /get_lun, /append
                printf, lun, date+': th-'+tprobe+' no orbit data ...'
                free_lun, lun
            endif
            return
        endif
        
        tuts = sfmepoch(dat.epoch,'unix')
        store_data, pre0+'pos_gsm', tuts, [[dat.xyz_gsm],[dat.dneuts]], limits=$
            {ytitle:'(Re)', colors:[6,4,2,0], labels:['GSM '+['X','Y','Z'],'D NeutS!N'], constant:[0,[-1,1]*dis0],panel_size:0.6}
        store_data, pre0+'mlt', tuts, dat.sm_lct_t, limits={labels:'MLT',ytitle:'(hr)',ystyle:1,yrange:[0,24],yticks:2,yminor:6,constant:[6,12,18],panel_size:0.4}
        store_data, pre0+'dis', tuts, dat.radius


        ; thx_[edsl1,emag].
        efsvar0s = pre0+['efs_dot0_dsl','efs_dot0_time']
        efsvar1s = ['efs_dsl','efs_utsec']
        efs = sread_thm_efi_l2(utr0, vars=efsvar0s, newname=efsvar1s, probes=tprobe)
        if size(efs,/type) ne 8 then begin
            if keyword_set(write_log) then begin
                openw, lun, logfn, /get_lun, /append
                printf, lun, date+': th-'+tprobe+' no EFS L2 data ...'
                free_lun, lun
            endif
            return
        endif
        tvar = pre0+'e_dsl1'
        uts = efs.efs_utsec
        edsl = efs.efs_dsl
        nrec = n_elements(edsl)/3
        if nrec le mintlen/dr0 then begin
            if keyword_set(write_log) then begin
                openw, lun, logfn, /get_lun, /append
                printf, lun, date+': th-'+tprobe+' EFS L2 less than '+sgnum2str(mintlen/3600d,nsgn=2)+' hr ...'
                free_lun, lun
            endif
            return
        endif
        edsl[*,2] = !values.d_nan
        ;for i=0,1 do edsl[*,i] = edsl[*,i]-smooth(edsl[*,i],twin/spinrate, /edge_truncate)
        store_data, tvar, uts, edsl[*,0:1], limits={ytitle:'(mV/m)', colors:rgb[0:1], labels:'DSL E0'+xyz[0:1]}
        emag = sqrt(edsl[*,0]^2+edsl[*,1]^2)
        store_data, pre0+'emag', uts, emag, limits={ytitle:'(mV/m)', labels:['|E|']}


        
        ; thx_[b_dsl].
        fgmvar0s = pre0+['fgs_dsl','fgs_time']
        fgmvar1s = ['fgs_dsl','fgs_utsec']
        fgm = sread_thm_fgm_l2(utr0, vars=fgmvar0s, newname=fgmvar1s, probes=tprobe)
        if size(fgm,/type) ne 8 then begin
            if keyword_set(write_log) then begin
                openw, lun, logfn, /get_lun, /append
                printf, lun, date+': th-'+tprobe+' no FGM L2 data ...'
                free_lun, lun
            endif
            return
        endif
        tvar = pre0+'b_dsl'
        tuts = fgm.fgs_utsec
        bdsl = fgm.fgs_dsl
        tnrec = n_elements(bdsl)/3
        if tnrec le mintlen/dr0 then begin
            if keyword_set(write_log) then begin
                openw, lun, logfn, /get_lun, /append
                printf, lun, date+': th-'+tprobe+' FGM L2 less than '+sgnum2str(mintlen/3600d,nsgn=2)+' hr ...'
                free_lun, lun
            endif
            return
        endif
        get_data, pre0+'dis', tmp, dis
        dis = sinterpol(dis, tmp, tuts)
        idx = where(dis gt dis0, cnt)
        if cnt ne 0 then brng = minmax(bdsl[idx,*]) else brng = minmax(bdsl)
        store_data, tvar, tuts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz, ystyle:0, yrange:brng}
        


    ;---Remove co-rotation field.
        ; thx_[e0,coef_e0].
        mindis0 = 3 ; Re.
        get_data, pre0+'dis', tuts, dis
        dis = interpol(dis, tuts, uts)
        win0 = double(dis le mindis0)
        winwd = 9600/spinrate
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
            i2 = (ti+drec[i])<(nrec-1)
            emagbg[ti] = min(emag[i1:i2])
        endfor
        ec = (1-emagbg/emag)
        store_data, pre0+'coef_e0', uts, ec
        store_data, pre0+'e0', uts, emagbg, limits=$
            {ytitle:'(mV/m)',labels:'|E| BG'}


        ; apply bg removal coefficient and update Euvw.
        get_data, pre0+'e_dsl1', limits=lim
        edsl[*,0] *= ec
        edsl[*,1] *= ec
        store_data, pre0+'e_dsl2', uts, edsl[*,0:1], limits=lim


        


    ;---Calc E dot0.
        ; thx_[b_dsl1,b0,db_dsl,e_dsl].
        bmin = 5d   ; nT.
        rmin = 0.25 ; percent.
        twin = 600d  ; sec. time window to remove offset.

        
        get_data, pre0+'b_dsl', tuts, bdsl
        bdsl = sinterpol(bdsl, tuts, uts)
        for i=0, 2 do bdsl[*,i] = smooth(bdsl[*,i], twin/spinrate, /edge_mirror)
        bmag = snorm(bdsl)
        edsl[*,2] = -(edsl[*,0]*bdsl[*,0]+edsl[*,1]*bdsl[*,1])/bdsl[*,2]
        ; remove those when E.B is dangerous.
        idx = where(bmag le bmin or abs(bdsl[*,2])/bmag le rmin, cnt)
        if cnt ne 0 then begin
            tdat = uts
            tdat[idx] = !values.d_nan
            
            ; ignore small gap, and small chunk of data.
            maxdrec = round(twin/dr0)
            mindrec = round(maxdrec/20)>2
            

            ; remove small chunk of data.
            tfs = where(finite(tdat), cnt)
            if cnt gt nrec*0.1 then begin
                dfs = tfs[1:cnt-1]-tfs[0:cnt-2]
                tmp = where(dfs ge mindrec, cnt)
                if cnt ne 0 then begin
                    idx2 = [tfs[tmp],nrec-1]
                    idx1 = [0,tfs[tmp+1]]
                    cnt = n_elements(idx1)
                    for i=0, cnt-1 do begin
                        if (idx2[i]-idx1[i]) le maxdrec then begin
                            tdat[idx1[i]:idx2[i]] = !values.d_nan
                        endif
                    endfor
                endif
            endif
            
            ; ignore small gaps.
            tfs = where(finite(tdat), cnt)
            if cnt gt nrec*0.1 then begin
                dfs = tfs[1:cnt-1]-tfs[0:cnt-2]
                tmp = where(dfs ge mindrec, cnt)
                if cnt ne 0 then begin
                    idx1 = tfs[tmp]+1
                    idx2 = tfs[tmp+1]-1
                    cnt = n_elements(idx1)
                    for i=0, cnt-1 do begin
                        if (idx2[i]-idx1[i]) le maxdrec then begin
                            tdat[idx1[i]:idx2[i]] = 0
                        endif
                    endfor
                endif
            endif
            
            idx = where(finite(tdat,/nan), cnt)
            if cnt ne 0 then edsl[idx,*] = !values.d_nan
        endif
        
        store_data, pre0+'e_dsl', uts, edsl, limits={ytitle:'(mV/m)', colors:rgb, labels:'DSL E'+xyz}
        store_data, pre0+'b0_dsl', uts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz}
        
        
<<<<<<< HEAD
        vmin = fmin
        vmax = fmax
        store_data, pre0+'ev_env', uts, [fmin+fmax], limits=$
            {yrange:[-1,1]*1}
        
        dmax = smooth(umax-vmax,padt/dr0)
        dmin = smooth(umin-vmin,padt/dr0)
        store_data, pre0+'dmax', uts, dmax, limits={yrange:[-1,1]*5}
        store_data, pre0+'dmin', uts, dmin, limits={yrange:[-1,1]*5}
        store_data, pre0+'euvw_env', uts, [[umin],[umax],[vmin],[vmax]], limits=$
            {ytitle:'(mV/m)',labels:['Eu','','Ev',''],colors:[6,6,2,2],panel_size:0.5}

        flags = (abs(dmax) ge 5) or (abs(dmin) ge 5)
        idx = where(flags eq 1, cnt)
        drec = padt/dr0
        if cnt ne 0 then begin
            for i=0, cnt-1 do flags[(idx[i]-drec)>0:(idx[i]+drec)<(nrec-1)] = 1
            idx = where(flags eq 1)
            get_data, pre0+'eff', uts, edsl
            edsl[idx,*] = !values.d_nan
            store_data, pre0+'eff', uts, edsl
        endif
        store_data, pre0+'bade_flag', uts, flags, limits=$
            {labels:'1: bad E', yrange:[-0.5,1.5], ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.2, ytitle:''}



    ;---E dot B = 0.
    ; thx_[bgsm,e0gsm].
        thm_load_fgm, probe=tprobe, trange=utr0, coord='dsl', datatype='fgl'
        get_data, pre0+'eff', uts, edsl, limits=lim
        get_data, pre0+'fgl', tuts, bdsl & bdsl=sinterpol(bdsl,tuts, uts)
        options, pre0+'fgl', 'colors', [6,4,2]
        edsl[*,2] = -(edsl[*,0]*bdsl[*,0]+edsl[*,1]*bdsl[*,1])/bdsl[*,2]
        bmag = snorm(bdsl)
        idx = where(bmag le 5 or abs(bdsl[*,2])/bmag le 0.25, cnt)
        if cnt ne 0 then edsl[idx,2] = !values.d_nan
        store_data, pre0+'e0gsm', uts, edsl, limits=$
            {ytitle:'(mV/m)', labels:'GSM E'+['x','y','z'], colors:[6,4,2]}
        thm_cotrans, pre0+'e0gsm', in_coord='dsl', out_coord='gsm'
        options, pre0+'e0gsm', 'ytitle', '(mV/m)'
=======

    ;---Rotate E/B fields to GSM, degap.
        thm_load_state, /get_support, probe=tprobe, trange=utr0
        thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'e_dsl', pre0+'e_gsm'
        thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'b_dsl', pre0+'b_gsm'
        store_data, pre0+'state*', /delete
>>>>>>> c08a1c821974d1b664af9993718ed6f9a52980f9
        
        tvar = pre0+'e_gsm'
        options, tvar, 'ytitle', '(mV/m)'
        options, tvar, 'colors', rgb
        options, tvar, 'labels', 'GSM E'+xyz
        
        tvar = pre0+'b_gsm'
        options, tvar, 'ytitle', '(nT)'
        options, tvar, 'colors', rgb
        options, tvar, 'labels', 'GSM B'+xyz
        
<<<<<<< HEAD
        options, pre0+'pos_gsm', 'panel_size', 0.6
        options, pre0+'mlt', 'panel_size', 0.4
        tdegap, pre0+'bgsm', /overwrite, /nowarning
        tdegap, pre0+'e0gsm', /overwrite, /nowarning
        
        get_data, pre0+'bgsm', tuts, bgsm
        get_data, pre0+'dis', tmp, dis & dis = interpol(dis,tmp, tuts)
        idx = where(dis ge 8, cnt)
        if cnt ne 0 then begin
            brng = minmax(bgsm[idx,*])
            options, pre0+'bgsm', 'yrange', brng
            options, pre0+'bgsm', 'ystyle', 0
        endif
        

=======
        get_data, pre0+'e_gsm', uts, egsm
        erng = max(abs(egsm),/nan)<80
        get_data, pre0+'dis', tuts, dis
        dis = interpol(dis, tuts, uts)
        idx = where(dis ge dis0, cnt)
        if cnt ne 0 then erng = max(abs(egsm[idx,*]),/nan)<80
        erng = [-1,1]*erng[0]
        options, pre0+'e_gsm', 'yrange', erng
        options, pre0+'e_gsm', 'ystyle', 0
        
        
        ;tplot, pre0+['e_dsl','e_gsm','b_dsl','b_gsm'], trange=utr0
    
>>>>>>> c08a1c821974d1b664af9993718ed6f9a52980f9

    ;---Plot.
        if ~keyword_set(noplots) then begin
            tdir = rootdir+'/th'+tprobe+'/'+time_string(ut0,tformat='YYYY')
            if file_test(tdir,/directory) eq 0 then file_mkdir, tdir
<<<<<<< HEAD
            ofn = tdir+'/'+date+'_th'+tprobe+'_bade_flag.pdf'
;            ofn = 0
;            stop
=======
            ofn = tdir+'/'+date+'_th'+tprobe+'_preprocess_field.pdf'
            ;ofn = 0
            
>>>>>>> c08a1c821974d1b664af9993718ed6f9a52980f9
            sgopen, ofn, xsize=11, ysize=8.5, /inch
            
            device, decomposed=0
            loadct2, 43
            xchsz = double(!d.x_ch_size)/!d.x_size
            ychsz = double(!d.y_ch_size)/!d.y_size
            
            tplot_options, 'ymargin', [5,5]
            tplot_options, 'xmargin', [15,15]
            
<<<<<<< HEAD
            figlabs = ['a','b','c','d','e','f','g']+'.'
            vars = pre0+['vsc','euvw_env','bade_flag','e0gsm','bgsm','pos_gsm','mlt']
=======
            figlabs = ['a','b','c','d','e','f']+'.'
            vars = pre0+['b_gsm','e_gsm','e_dsl1','pos_gsm','mlt']
>>>>>>> c08a1c821974d1b664af9993718ed6f9a52980f9
            nvar = n_elements(vars)
            
            titl = 'TH-'+strupcase(tprobe)+' '+date+' E field preprocess: remove co-rotation E, and E.B=0'
            tplot, vars, trange=utr0, get_plot_position=poss
            xyouts, (poss[0,0]+poss[2,0])*0.5, poss[3,0]+ychsz*0.8, /normal, alignment=0.5, titl, charsize=1.2
            for i=0, nvar-1 do xyouts, poss[0,i]-xchsz*8, poss[3,i]-ychsz*0.5, figlabs[i], /normal

            sgclose
        endif
        
    endforeach
end


secofday = 86400d
probes = ['b']
utr0 = time_double(['2012-03-25','2015-12-31'])

; 2013-03-27. SDT
; need to recover the SDT and eclipse flag.
; and manuver flag??

nday = (utr0[1]-utr0[0])/secofday
uts = smkarthm(utr0[0], secofday, nday+1, 'x0')

foreach tprobe, probes do foreach tut, uts do themis_preprocess_efield, tut, probes=tprobe, /write_log

end
