

; load data.
_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 5
tplot_options, 'yticklen', -0.01
tplot_options, 'xticklen', -0.02
tplot_options, 'zcharsize', 0.8
tplot_options, 'ycharsize', 1



; constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1

; settings.
    ionminen = 30   ; eV.
    eleminen = 200  ; eV.

    models = ['t89','t04s','t01','t96']
    modelidx = 2
    probes = ['a','b']
    nprobe = n_elements(probes)

    utr1 = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    utr2 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
    midut = time_double('2013-06-07/04:55')

    diputa = time_double('2013-06-07/04:54:30')
    diputb = time_double('2013-06-07/04:54:59')


; **** calc keogram and photon count.
    yr = [60,68]
    mlat0s = smkarthm(yr[0],yr[1],0.1,'dx')
    nmlat0 = n_elements(mlat0s)

    get_data, 'asf_info', tmp, info
    mlats = info.mlats
    mlts = info.mlts
    imgsz = info.imgsz
    get_data, 'asf_mos', uts, mos, pxidx
    nrec = n_elements(uts)

    ; make keogram.
    if tnames('rbspa_keogram_t01') eq '' then begin
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            get_data, pre0+'fpt_mlat', tuts, scmlats
            get_data, pre0+'fpt_mlt', tuts, scmlts
            
            scmlats = reform(scmlats[*,modelidx])
            scmlts = reform(scmlts[*,modelidx])
            
            scmlats = interpol(scmlats,tuts,uts)
            scmlts = interpol(scmlts,tuts,uts)*15-360
            
            keos = fltarr(nrec,nmlat0)
            pcnts = fltarr(nrec)
            for i = 0, nrec-1 do begin
                idx = where(abs(mlts-scmlts[i]) le 0.1, cnt)
                if cnt lt 3 then continue
                timg = fltarr(imgsz)
                timg[pxidx] = mos[i,*]
                
                timg = timg[idx]
                tmlats = mlats[idx]
                idx = sort(tmlats)
                timg = timg[idx]
                tmlats = tmlats[idx]
                
                keos[i,*] = interpol(timg,tmlats,mlat0s)
                
                ; add a line for Fpt.
                ;tmp = min(abs(mlat0s-scmlats[i]), idx)
                ;keos[i,idx] = !values.f_nan

                
                ; photon count around fpt.
                idx = where(abs(mlts-scmlts[i]) le 0.1 and $
                    abs(mlats-scmlats[i]) le 0.1, cnt)
                timg = fltarr(imgsz)
                timg[pxidx] = mos[i,*]
                if cnt ge 0 then pcnts[i] = max(timg[idx])
                
            endfor
            
            tvar = (pre0+'keogram_'+models[modelidx])[0]
            store_data, tvar, uts, keos, mlat0s, $
                limits = {spec:1, no_interp:1, $
                ytitle:'MLat!C(deg)', $
                yrange:yr, ztitle:'Photon count', yticks:2, zrange:[50,300], zticks:1, ystyle:1}
            tvar = (pre0+'count_'+models[modelidx])[0]
            store_data, tvar, uts, pcnts, $
                limits = {ytitle:'Photon count',yrange:[100,600], ylog:1, yminor:5, labels:'Photon!C  @SC footpnt'}
                
        endforeach
    endif


; **** simple model to calculate the source and time lag.
    if tnames('rbspa_o_alpha') eq '' then begin
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            
            get_data, pre0+'o_pa', uts, dat, val
            nrec = n_elements(uts)
            
            alphas = dblarr(nrec)
            for i = 0, nrec-1 do begin
                tmp = max(reform(dat[i,*]), idx)
                alphas[i] = val[idx]
            endfor
            store_data, pre0+'o_alpha', uts, alphas
            
            
            get_data, pre0+'o_en', uts, dat, val
            nrec = n_elements(uts)
            
            ens = dblarr(nrec)
            for i = 0, nrec-1 do begin
                tidx = where(val[i,*] ge 100)
                tmp = max(reform(dat[i,tidx]), idx)
                ens[i] = val[i,tidx[idx]]
            endfor
            v0s = sqrt(2*ens*1.6e-19/(16*1.67e-27))*1e-3    ; in km/s.
            store_data, pre0+'o_energy', uts, ens
            
            
            
            ; trace use the middle time and pos.
            tet = stoepoch(midut,'unix')
            ; get sc position.
            get_data, pre0+'pos_gsm', tuts, dat
            tpos = sinterpol(dat, tuts, midut)
            get_data, 't01_par', tuts, pars
            par = sinterpol(pars, tuts, midut)
            
            r0 = 2d
            geopack_epoch, tet, year, mo, dy, hr, mi, sc, /breakdown_epoch
            geopack_recalc, year, mo, dy, hr, mi, sc, /date, tilt = tilt
            geopack_trace, tpos[0], tpos[1], tpos[2], -1, par, xf, yf, zf, $
                fline = fline, /igrf, r0 = r0, /t01
            geopack_igrf_gsm, fline[*,0], fline[*,1], fline[*,2], bx0, by0, bz0
            geopack_t01, par, fline[*,0], fline[*,1], fline[*,2], dbx, dby, dbz
            
            b0 = 240d 
            bgsms = [[bx0+dbx],[by0+dby],[bz0+dbz]]
            bmags = snorm(bgsms)
            bcoef = 1+(b0-min(bmags))/bmags
            for i = 0, 2 do bgsms[*,i] *= bcoef
            bmags = snorm(bgsms)
            
            bs = b0/sin(alphas*rad)^2
            
            r0s = dblarr(nrec,3)
            t0s = dblarr(nrec)
            for i = 0, nrec-1 do begin
                for j = 0, 2 do r0s[i,j] = sinterpol(fline[*,j], bmags, bs[i])
                idx = where(fline[*,0] lt r0s[i,0])
                tfline = [fline[idx,*],r0s[i,*]]
                tbmags = bmags[idx]
                tnrec = n_elements(tfline[*,0])
                t0s[i] = 0
                for j = 1, tnrec-2 do begin
                    tds = snorm(tfline[j,*]-tfline[j-1,*])*re
                    tdt = tds/(v0s[i]*sqrt(1-tbmags[j-1]^2/bs[i]^2))
                    t0s[i] = t0s[i]+tdt
                endfor
            endfor
            
            t1s = smooth(t0s,3)
            r1s = snorm(r0s)
            get_data, pre0+'o_en', uts, dat
            for i = 0, nrec-1 do begin
                if max(dat[i,*]) lt 0 then begin
                    t1s[i] = !values.d_nan
                    r1s[i] = !values.d_nan
                endif
            endfor
            
            store_data, pre0+'dt', uts, t1s, limits = {psym:-1, ytitle:'Time lag!C(sec)'}
            store_data, pre0+'ds', uts, r1s, limits = {psym:-1, ytitle:'Energized R!C(Re)',yrange:[1,6]}

        endforeach
    endif

; **** integrate spec to get ke flux.
    if tnames('rbspa_keflux_oxygen') eq '' then begin
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            plot_hope_l3_keflux, utr2, 'oxygen', probe = tprobe, min_energy = ionminen
        endforeach
    endif


; **** calc mapped pflux, keflux.
    pflabs = 'S!D'+['||','!9^!X,West','!9^!X,North']

    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        get_data, pre0+'pf_fac', uts, dat, limits = lim
        get_data, pre0+'map_coef', tuts, tmp
        tmp = interpol(tmp[*,modelidx],tuts, uts)
        for i = 0, 2 do dat[*,i] *= tmp
        store_data, pre0+'pf_fac_mor_map', uts, dat, limits = lim
        options, pre0+'pf_fac_mor_map', 'labels', pflabs
        
        store_data, pre0+'pf_fac_mor_map', uts, dat[*,0], limits = $
            {ytitle:'(mW/m!U2!N)', labels:pflabs[0]+'!N!C@100km', colors:0}
        
        get_data, pre0+'keflux_oxygen', uts, dat, limits = lim
        get_data, pre0+'map_coef', tuts, tmp
        tmp = interpol(tmp[*,modelidx],tuts, uts)
        dat = dat*tmp
        store_data, pre0+'keflux_oxygen_map', uts, dat, limits = $
            {ytitle:'(mW/m!U2!N)', labels:'KEflux O!U+!N!C  @100km', colors:0}
        
    endforeach



; **** options for variables.
    vars = 'rbsp?_keogram_t01'
    options, vars, 'ytitle', 'MLat (deg)'
    options, vars, 'yminor', 4

    vars = tnames('rbsp?_count_t01')
    options, vars, 'ylog', 1
    options, vars, 'yrange', [100,800]
    options, vars, 'ystyle', 1
    options, vars, 'constant', 0
    options, vars, 'ytitle', ''
    options, vars, 'labels', 'Fpt T01'
    options, vars, 'yticks', 1
    options, vars, 'yminor', 2
    options, vars, 'ytickv', [100,400]
    options, vars, 'ytickname', ['100','400']

    tvar = 'rbspa_keflux_oxygen_map'
    options, tvar, 'yrange', [-3.5,0.5]
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-3,-1.5,0]
    options, tvar, 'yminor', 5
    options, tvar, 'ystyle', 1
    options, tvar, 'psym', -1
    options, tvar, 'symsize', 0.5

    tvar = 'rbspb_keflux_oxygen_map'
    options, tvar, 'yrange', [-9.5,0.5]
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-9,-4.5,0]
    options, tvar, 'yminor', 5
    options, tvar, 'ystyle', 1
    options, tvar, 'psym', -1
    options, tvar, 'symsize', 0.5

    tvar = tnames('rbsp?_dt')
    get_data, tvar, uts, dat
    store_data, tvar, uts, dat/60d
    options, tvar, 'ytitle', '(min)'
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [0,4,8]
    options, tvar, 'yminor', 4

    tvar = tnames('rbsp?_ds')
    options, tvar, 'ytitle', '(Re)'
    options, tvar, 'labels', 'R!Dmax'
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [1,4,7]
    options, tvar, 'yminor', 3
    
    tvar = tnames('rbsp?_pf_fac_mor_map')
    options, tvar, 'labels', 'Map to!C100 km'
    options, tvar, 'constant', !values.d_nan

    tvar = tnames('rbsp?_keflux_oxygen_map')
    options, tvar, 'labels', 'Map to!C  100 km'
    
    
    vars = tnames('rbsp?_o_en')
    options, vars, 'ytitle', '(eV)'
    options, vars, 'yrange', [40,4e4]
    options, vars, 'ytickv', [100,1000,10000]
    options, vars, 'ystyle', 1
    options, vars, 'yminor', 5
    options, vars, 'yticks', 2
    options, vars, 'ytickname', '1e'+['2','3','4']
    options, vars, 'zrange', [1e4,1e7]
    options, vars, 'zticks', 3
    options, vars, 'ztickname', '1e'+['4','5','6','7']
    
    vars = tnames('rbsp?_o_pa')
    options, vars, 'zrange', [1e5,1e8]
    options, vars, 'zticks', 3
    options, vars, 'ztickname', '1e'+['5','6','7','8']
    
    

; **** tplot variable and other settings.
    ;vars = ['keogram_t01','count_t01','pf_fac_mor_map','keflux_oxygen_map','o_en','o_pa']
    vars = ['keogram_t01','pf_fac_mor_map','keflux_oxygen_map','o_en','o_pa']
    ;vars = ['count_t01','pf_fac_mor_map','keflux_oxygen_map','o_en','o_pa']
    logidx = (where(vars eq 'pf_fac_mor_map'))[0]
    nvar = n_elements(vars)

    labs = ['mlt','lshell','mlat']


    tposs = sgcalcpos(nvar, region=pos1)

    ofn = 0
    ofn = shomedir()+'/fig_new_pflux_compare.pdf'
    ;sgopen, ofn, xsize = 8.5, ysize = 11, /inch
    ;sgopen, ofn, xsize = 5, ysize = 9, /inch
    sgopen, ofn, xsize = 10, ysize = 4.5, /inch
    device, decomposed = 0
    loadct2, 43

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    utr = utr1
    lmarg = 15
    cols = [6,4,2]
    pos1 = [0.1,0.5,1,1]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1
    pos2 = [0.1,0,1,0.5]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1
    pos1 = [0.01,0,0.5,1.05]
    pos2 = [0.5,0,0.99,1.05]
    

    pre0 = 'rbspa_'
    diput = diputa
    tposs = sgcalcpos(nvar, region=pos1)
    tplot, pre0+vars[0:logidx-1], trange = utr, /nouttick, $
        /noerase, position = tposs[*,0:logidx-1], vlab_margin = lmarg, title = 'RBSP-A'
    stplot_linlog, pre0+vars[logidx], trange = utr, $
        /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 1, $
        linyr = [-5,5], logyr = [5,200], logytickv = [10,100], linytickv = [-5,0,5]
    tplot, pre0+vars[logidx+1:*], trange = utr, /noerase, $
        position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs
    get_data, pre0+'fpt_mlat', uts, dat
    dat = reform(dat[*,modelidx])
    plot, uts, dat, /noerase, position=tposs[*,0], $
        xstyle=5, xrange=utr, $
        ystyle=5, yrange=[60,68], $
        color=sgcolor('white'), thick=4


    lab2s = ['Keogram!C       @SC MLT', $
        'S parallel!C       1/4-1200s', 'KEflux O!U+','O!U+!N EN','O!U+!N PA']
    ;lab2s = ['Photon Count!C       @SC Fpt', $
    ;    'S parallel!C       1/4-1200s', 'KEflux O!U+','O!U+!N EN','O!U+!N PA']
    lab1s = ['a','b','c','d','e']+'-1.'
    for i = 0, nvar-1 do begin
        xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
            lab1s[i]+' '+lab2s[i]
    endfor

    ; add notations.
    ttut = diputa
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty1 = tposs[3,0]
    ty2 = tposs[1,nvar-1]
    plots, tx+[0,0], [ty1,ty2], /normal, color = 6
    tidx = where(vars eq 'keflux_oxygen_map')
    ty = tposs[1,tidx]+ychsz*1.5
    xyouts, tx+xchsz*0.5, ty, /normal, alignment=0, 'T_Dip', color = 0

    get_data, pre0+'keflux_oxygen_map', tuts, dat
    idx = where(tuts ge utr[0] and tuts le utr[1])
    tuts = tuts[idx]
    dat = dat[idx]
    tmp = min(dat, idx)
    okeut = tuts[idx]
    ttut = okeut
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty1 = tposs[3,0]
    ty2 = tposs[1,nvar-1]
    plots, tx+[0,0], [ty1,ty2], /normal, color = 6, linestyle = 2

    tidx = where(vars eq 'keflux_oxygen_map')
    tpos = tposs[*,tidx]
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, alignment = 0, 'Upward', color = 0
    ttut = okeut
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty = tposs[1,tidx]+ychsz*0.6
    xyouts, tx+xchsz*1, ty, /normal, alignment = 0, 'T_Max|KE|', color = 0
    
    tidx = where(vars eq 'pf_fac_mor_map')
    tpos = tposs[*,tidx]
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, alignment = 0, 'Earthward', color = 0
    plots, tpos[[0,2]], tpos[1]+(tpos[3]-tpos[1])*0.25, /normal, linestyle = 5

    tidx = where(vars eq 'o_pa') 
    ttx = tpos[0,0]+xchsz*1
    tty = tposs[1,tidx]+ychsz*0.5
    get_data, pre0+'dt', tuts, dat
    idx = where(tuts eq okeut)
    tdt = dat[idx]
    get_data, pre0+'ds', tuts, dat
    idx = where(tuts eq okeut)
    tds = dat[idx]
        xyouts, ttx, tty+ychsz*0, /normal, alignment=0, color=255, charsize=0.9, $
        'dT = '+sgnum2str(tdt/60d,ndec=1)+' min, Source below '+sgnum2str(tds,ndec=1)+' Re'
    xyouts, ttx, tty+ychsz*1, /normal, alignment=0, color=255, charsize=0.9, $
        'T_Max|KE| - T_Dip = '+sgnum2str((okeut-diput)/60,ndec=1)+' min'
        
        
        
        
    pre0 = 'rbspb_'
    diput = diputb
    tposs = sgcalcpos(nvar, region=pos2)
    tplot, pre0+vars[0:logidx-1], trange = utr, /nouttick, $
        /noerase, position = tposs[*,0:logidx-1], vlab_margin = lmarg, title = 'RBSP-B'
    stplot_linlog, pre0+vars[logidx], trange = utr, $
        /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 1, $
        linyr = [-5,5], logyr = [5,200], logytickv = [10,100], linytickv = [-5,0,5]
    tplot, pre0+vars[logidx+1:*], trange = utr, /noerase, $
        position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs
    get_data, pre0+'fpt_mlat', uts, dat
    dat = reform(dat[*,modelidx])
    plot, uts, dat, /noerase, position=tposs[*,0], $
        xstyle=5, xrange=utr, $
        ystyle=5, yrange=[60,68], $
        color=sgcolor('white'), thick=4
    
    ;lab2s = ['Keogram!C      @SC MLT', 'Brightness!C      @SC Fpt', $
    ;    'S FAC', 'KEflux O!U+','O!U+!N EN','O!U+!N PA']
    ;lab1s = ['a','b','c','d','e','f']+'-1.'
    ;lab2s = ['Keogram!C      @SC MLT', $
    ;    'S FAC', 'KEflux O!U+','O!U+!N EN','O!U+!N PA']        
    lab1s = ['a','b','c','d','e']+'-2.'
    for i = 0, nvar-1 do begin
        xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
            lab1s[i]+' '+lab2s[i]
    endfor
    
    ; add notations.
    ttut = diputb
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty1 = tposs[3,0]
    ty2 = tposs[1,nvar-1]
;    plots, tx+[0,0], [ty1,ty2], /normal, color = 6
    tidx = where(vars eq 'keflux_oxygen_map')
    ty = tposs[1,tidx]+ychsz*1.5
;    xyouts, tx+xchsz*0.5, ty, /normal, alignment=0, 'T_Dip', color = 0
    
    get_data, pre0+'keflux_oxygen_map', tuts, dat
    idx = where(tuts ge utr[0] and tuts le utr[1])
    tuts = tuts[idx]
    dat = dat[idx]
    tmp = min(dat, idx)
    okeut = tuts[idx]
    ttut = okeut
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty1 = tposs[3,0]
    ty2 = tposs[1,nvar-1]
    plots, tx+[0,0], [ty1,ty2], /normal, color = 6, linestyle = 2
    
    tidx = where(vars eq 'keflux_oxygen_map')
    tpos = tposs[*,tidx]
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, alignment = 0, 'Upward', color = 0
    ttut = okeut
    tx = (ttut-utr[0])/(utr[1]-utr[0])
    tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
    ty = tposs[1,tidx]+ychsz*0.6
    xyouts, tx+xchsz*1, ty, /normal, alignment = 0, 'T_Max|KE|', color = 0

    tidx = where(vars eq 'pf_fac_mor_map')
    tpos = tposs[*,tidx]
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, alignment = 0, 'Earthward', color = 0
    plots, tpos[[0,2]], tpos[1]+(tpos[3]-tpos[1])*0.25, /normal, linestyle = 5
    
    tidx = where(vars eq 'o_pa')
    ttx = tpos[0,0]+xchsz*1
    tty = tposs[1,tidx]+ychsz*0.5
    get_data, pre0+'dt', tuts, dat
    idx = where(tuts eq okeut)
    tdt = dat[idx]
    get_data, pre0+'ds', tuts, dat
    idx = where(tuts eq okeut)
    tds = dat[idx]
;    xyouts, ttx, tty+ychsz*0, /normal, alignment=0, color=255, charsize=0.9, $
;        'dT = '+sgnum2str(tdt/60d,ndec=1)+' min, Source below '+sgnum2str(tds,ndec=1)+' Re'
;    xyouts, ttx, tty+ychsz*1, /normal, alignment=0, color=255, charsize=0.9, $
;        'T_Max|KE| - T_Dip = '+sgnum2str((okeut-diput)/60,ndec=1)+' min'
    xyouts, ttx, tty+ychsz*0, /normal, alignment=0, color=255, charsize=0.9, $
        'Source below '+sgnum2str(tds,ndec=1)+' Re'


    sgclose


end
