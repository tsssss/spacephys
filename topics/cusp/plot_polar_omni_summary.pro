
pro catdata, vname, flag, t0, omni, utr
    vnames = strlowcase(tag_names(omni))
    i = where(vnames eq vname)
    if flag eq 1 then store_data, vname, data = {x:t0, y:omni.(i)} $
    else begin
        get_data, vname, data = tmp
        idx = where(tmp.x ge utr[0])
        store_data, vname, data = {x:[tmp.x[idx[0]:*],t0], $
            y:[tmp.y[idx[0]:*,*],omni.(i)]}
    endelse
end

pro plot_polar_omni_summary
    compile_opt idl2
    
    ; plot options.
    sgwindow, 0, xsize = 720, ysize = 852, xpos = 0, ypos = 25
    sgindexcolor, ct = 43
    !p.font = 1 & !y.charsize = 0.7 & !x.charsize = 0.7 &!z.charsize = 0.6
    tplot_options, 'zcharsize', 0.7

    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    time_stamp, /off
    
    ; time span.
    et1 = stoepoch('2002-09-05')
    et2 = stoepoch('2002-09-06')
    
    ; use polar pos to detect orbital period and plot time range.
    poet = sepochfloor(et1)         ; used to read polar data.
    omniet = sepochfloor(et1)       ; used to read omni data.
    trids = 1 & et = 0 & mlat = 0 & ilat = 0 & mlt = 0 & dis = 0

    ; loop each orbit.
    go = 1
    while go do begin
        ; check if need new buffer.
        if n_elements(trids) eq 1 then begin
            tmp = sread_polar_pos(t = poet, $
                vars = ['Epoch','mlt','ilat','mlat','dis'])
            poet += 86400000D
            if poet gt et2 then go = 0
            trids = [trids, sfindnode(tmp.mlat,2)+n_elements(et)]
            et = [et, tmp.epoch]
            mlat = [mlat, abs(tmp.mlat)]
            ilat = [ilat, abs(tmp.ilat)]    ; abs ilat.
            mlt = [mlt, tmp.mlt]
            dis = [dis, tmp.dis]
        endif
        
        ; get time and pos info.
        tr = et[trids[0:1]]
        utr = sfmepoch(tr, 'unix')
        t0 = sfmepoch(et[trids[0]:trids[1]], 'unix')
        print, 'time range: ', sfmepoch(tr[0]), ', ', sfmepoch(tr[1])
        store_data, 'lat', data = {x:t0, y:[[ilat[trids[0]:trids[1]]],$
            [mlat[trids[0]:trids[1]]]]}
        store_data, 'mlt', data = {x:t0, y:mlt[trids[0]:trids[1]]}
        store_data, 'dis', data = {x:t0, y:dis[trids[0]:trids[1]]}
        options, 'lat', 'ytitle', 'Lat (deg)'
        options, 'lat', 'labels', ['ILat', 'MLat']
        options, 'lat', 'colors', [4,2]
        options, 'mlt', 'ytitle', 'MLT (hr)'
        options, 'dis', 'ytitle', 'R (Re)'
        ylim, 'lat', 0, 90, 0
        options, 'lat', 'labpos', [30,60]
        ylim, 'mlt', 0, 24, 0
        ylim, 'dist', 2, 9, 0
        
        ; prepare next round.
        et = et[trids[1]:*]
        mlat = mlat[trids[1]:*]
        ilat = ilat[trids[1]:*]
        mlt = mlt[trids[1]:*]
        dis = dis[trids[1]:*]
        trids = trids[1:*]-trids[1]
        
        ; read polar hydra.
        tmp = sread_polar_hydra(utr)
        if size(tmp, /type) eq 8 then begin     ; is structure, not -1.
            t0 = sfmepoch(tmp.epoch, 'unix')
            store_data, 'JEi', data = {x:t0, y:tmp.jei, v:tmp.eni}
            store_data, 'JEe', data = {x:t0, y:tmp.jee, v:tmp.ene}
            options, ['JEi','JEe'], 'ztitle', '1/(cm^2-s-sr)'
            options, ['JEi','JEe'], 'spec', 1
            options, 'JEi', 'ytitle', 'JEi (eV)'
            options, 'JEe', 'ytitle', 'JEe (eV)'
            ylim, 'JEi', tmp.eni[0], tmp.eni[-1], 1
            ylim, 'JEe', tmp.ene[0], tmp.ene[-1], 1
            zlim, 'JEi', 1e4, 1e8, 1
            zlim, 'JEe', 1e4, 5e8, 1
        endif
        
        ; read omni.
        tmp = tplot_tr('flow', /full)
        ; check if need to load next month.
        if n_elements(tmp) eq 1 then flag = 1 $     ; no old data.
        else if tmp[1] eq utr[1] then flag = 1 $    ; no old data.
        else if tmp[1] gt utr[1] then flag = -1 $   ; no new data.
        else flag = 0                               ; new data.
        if flag gt -1 then begin
            omniut = sfmepoch(omniet,'unix')
            omni = sread_omni(omniut)
            omniet = sepochadd(omniet, 1, 'mo')
            t0 = sfmepoch(omni.epoch, 'unix')
            catdata, 'ae', flag, t0, omni, utr
            catdata, 'symh', flag, t0, omni, utr
            catdata, 'pdyn', flag, t0, omni, utr
            catdata, 'bgse', flag, t0, omni, utr
            catdata, 'vgse', flag, t0, omni, utr
            catdata, 'flow', flag, t0, omni, utr
            
            options, 'ae', 'ytitle', 'AE (nT)'
            options, 'symh', 'ytitle', 'SymH (nT)'
            options, 'pdyn', 'ytitle', 'P (nPa)'
            options, 'bgse', 'ytitle', 'B GSE (nT)'
            options, 'vgse', 'ytitle', 'V GSE (km/s)'
            options, 'flow', 'ytitle', 'V (km/s)'
            options, 'bgse', 'labels', ['x','y','z']
            options, 'bgse', 'colors', [6,4,2]
            get_data, 'vgse', t0, tmp & store_data, 'vgse1', t0, tmp[*,1:2]
            options, 'vgse1', 'labels', ['y','z']
            options, 'vgse1', 'colors', [4,2]
        endif

        ; prepare plot.
        dirptn = shomedir()+'/omni_polar_halforbit/yyyy/MM/'
        fnptn = 'omnipolar_yyyy_MMdd_hh.eps'
        dir = sptn2fn(dirptn, tr[0])
        if file_test(dir) eq 0 then file_mkdir, dir
        fn = dir+sptn2fn(fnptn, tr[0])
        tl = 'omni+polar. '+time_string(t0[0],tformat='YYYY-MM-DD')
        vars = ['bgse', 'vgse1', 'pdyn', 'ae', 'symh', $
            'JEi', 'JEe', 'lat', 'mlt', 'dis']
        lbls = ['flow']
        
        tplot, vars, var_label = lbls, trange = utr, title = tl
        pstplot, filename = fn
    endwhile
end
