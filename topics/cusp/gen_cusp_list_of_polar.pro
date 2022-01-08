;+
; Type: procedure.
; Purpose: Loop through 1 month of Polar orbit and record all cusp crossings.
; Parameters:
;   etr, in, double/dblarr[2], req. Time range in epoch. If in double, it is
;       the start time, the end time is automatically set to be 1 month later. 
;       If in dblarr[2], it is in [start, end time].
;   logfile, in, string, opt. Log file's name. Default is 
;       '~/cusp_list_polar_yyyy_MM.log'
; Keywords:
;   datadir, in, string, opt. Root directory for Sheng's data file system.
;   plotdir, in, string, opt. Root directory for the plots.
; Return: none.
; Notes: none.
; Dependence: none.
; History:
;   2014-07-15, Sheng Tian, create.
;-
pro gen_cusp_list_of_polar, etr, logfile, datadir = datadir, plotdir = plotdir
    compile_opt idl2
    
    ; const.
    re = 6378d
    
    ; plot options.
    sgwindow, 0, xsize = 720, ysize = 852, xpos = 0, ypos = 25
    !p.font = 1 & !y.charsize = 0.7 & !x.charsize = 0.7 &!z.charsize = 0.6
    tplot_options, 'zcharsize', 0.7

    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', -1
    time_stamp, /off
    
    ; time span.
    if n_elements(etr) eq 1 then etr = [etr,sepochadd(etr,1,'mo')]
    if n_elements(etr) ne 2 then message, 'invalid time range ...'
    
    ; directories.
    if n_elements(plotdir) eq 0 then plotdir = shomedir()+'/'+sfmepoch(etr[0],'yyyy_MM')+'/'
    if file_test(plotdir) eq 0 then file_mkdir, plotdir

    ; log file.
    if n_elements(logfile) eq 0 then $
        logfile = shomedir()+'/cusp_list_polar_'+sfmepoch(etr[0],'yyyy_MM')+'.log'
    if file_test(logfile) eq 0 then begin
        spawn, 'touch '+logfile
        openw, lun, logfile, /get_lun
        printf, lun, ' Polar cusp date/time           Polar position          SymH   AE    FAST'
        printf, lun, '----------------------  ------------------------------  ----  ----  ------'
        printf, lun, 'yyyy-mmdd  hh:mm-hh:mm   R (Re)  ILat (deg)   MLT (hr)  (nT)  (nT)  orbit#'
        printf, lun, '---------  -----------  -------  ----------  ---------  ----  ----  ------'
        free_lun, lun
    endif
    
    ; use polar pos to detect orbital period and plot time range.
    tet = sepochfloor(etr[0])         ; used to read polar data.
    ; read the previous day's orbit to find the first node.
    tmp = sread_polar_orbit(sfmepoch(sepochadd(tet,-1,'dy'),'unix'), $
            vars = ['Epoch','MAG_LATITUDE'])
    nodeids = sfindnode(tmp.mag_latitude,2)
    nodes = [tmp.epoch[nodeids]]
    while tet lt etr[1] do begin
        tmp = sread_polar_orbit(sfmepoch(tet,'unix'), $
            vars = ['Epoch','MAG_LATITUDE'])
        tet += 86400000d
        nodeids = sfindnode(tmp.mag_latitude,2)
        nodes = [nodes,tmp.epoch[nodeids]]
    endwhile
    
    ; loop through adjacent nodes.
    nodes = suniq(nodes)
    nodes = nodes[sort(nodes)]
    for i = 0, n_elements(nodes)-2 do begin
        tetr = nodes[i:i+1]
        ; check starting time of nodes.
        if tetr[1] le etr[0] then continue
        
        utr = sfmepoch(tetr, 'unix')
        
        ; read polar hydra.
        tmp = sread_polar_hydra(utr)
        if size(tmp, /type) eq 8 then begin
            t0 = sfmepoch(tmp.epoch, 'unix')
            store_data, 'poJEi', t0, tmp.jei, tmp.eni
            store_data, 'poJEe', t0, tmp.jee, tmp.ene
            ylim, 'poJEi', tmp.eni[0], tmp.eni[-1], 1
            ylim, 'poJEe', tmp.ene[0], tmp.ene[-1], 1
        endif
        options, 'po'+['JEi','JEe'], 'ztitle', '1/(cm!U2!N-s-sr)'
        options, 'po'+['JEi','JEe'], 'spec', 1
        options, 'po'+['JEi','JEe'], 'no_interp', 1
        options, 'poJEi', 'ytitle', 'Polar ion (eV)'
        options, 'poJEe', 'ytitle', 'Polar ele (eV)'
        zlim, 'poJEi', 1e4, 1e8, 1
        zlim, 'poJEe', 1e4, 5e8, 1
        
        ; read fast esa.
        tmp = sread_fast_esa(utr, type = 'ies')
        if size(tmp, /type) eq 8 then begin
            store_data, 'faJEi', sfmepoch(tmp.epoch, 'unix'), $
                (tmp.ion_0+tmp.ion_90+tmp.ion_180)/3, tmp.ion_en
        endif
        tmp = sread_fast_esa(utr, type = 'ees')
        if size(tmp, /type) eq 8 then begin
            store_data, 'faJEe', sfmepoch(tmp.epoch, 'unix'), $
                (tmp.el_0+tmp.el_90+tmp.el_180)/3, tmp.el_en
        endif
        options, 'fa'+['JEi','JEe'], 'ztitle', '1/(cm!U2!N-s-sr)'
        options, 'fa'+['JEi','JEe'], 'spec', 1
        options, 'fa'+['JEi','JEe'], 'no_interp', 1
        options, 'faJEi', 'ytitle', 'Fast ion (eV)'
        options, 'faJEe', 'ytitle', 'Fast ele (eV)'
        ylim, 'faJEi', 4, 40000, 1
        ylim, 'faJEe', 4, 40000, 1
        zlim, 'faJEi', 1e4, 1e8, 1
        zlim, 'faJEe', 1e6, 1e10, 1
        
        ; read omni.
        tmp = sread_omni(utr)
        t0 = sfmepoch(tmp.epoch, 'unix')
        store_data, 'ae', t0, tmp.ae
        store_data, 'symh', t0, tmp.symh
        store_data, 'pdyn', t0, tmp.pdyn
        store_data, 'bgse', t0, tmp.bgse
        store_data, 'flow', t0, tmp.vgse[*,0]
        store_data, 'vgse', t0, tmp.vgse[*,1:2]
        options, 'ae', 'ytitle', 'AE (nT)'
        options, 'ae', 'format', '(I5)'
        options, 'symh', 'ytitle', 'SymH (nT)'
        options, 'pdyn', 'ytitle', 'P (nPa)'
        options, 'bgse', 'ytitle', 'B GSE (nT)'
        options, 'vgse', 'ytitle', 'V GSE (km/s)'
        options, 'flow', 'ytitle', 'V (km/s)'
        options, 'flow', 'format', '(F6.1)'
        options, 'bgse', 'labels', ['x','y','z']
        options, 'bgse', 'colors', [6,4,2]
        options, 'vgse', 'labels', ['y','z']
        options, 'vgse', 'colors', [4,2]
        
        ; read polar pos.
        vars = ['Epoch','GSE_POS','EDMLT_TIME','L_SHELL','MAG_LATITUDE']
        tmp = sread_polar_orbit(utr, vars = vars)
        t0 = sfmepoch(tmp.epoch, 'unix')
        store_data, 'pomlt', t0, tmp.edmlt_time
        ilat = acos(sqrt(1d/tmp.l_shell))*(180d/!dpi)
        idx = where(tmp.mag_latitude lt 0, cnt)
        if cnt gt 0 then ilat[idx] = -ilat[idx]
        dis = sqrt(total(tmp.gse_pos^2,2))*(1d/re)
        store_data, 'pohem', t0, ((ilat ge 0)*2-1)*90, $
            limits = {thick:4, colors:2}
        store_data, 'poilat', t0, abs(ilat)
        store_data, 'podist', t0, dis
        tvars = 'po'+['mlt','ilat','dist']
        options, tvars, 'colors', 2
        options, tvars, 'labels', 'polar'
        options, tvars, 'labflag', 3
        
        ; read fast pos.
        tmp = sread_fast_pos(utr, vars = ['Epoch','mlt','ilat','dis','orbit'])
        t0 = sfmepoch(tmp.epoch, 'unix')
        store_data, 'famlt', t0, tmp.mlt
        store_data, 'fahem', t0, ((tmp.ilat ge 0)*2-1)*60, $
            limits = {thick:4, colors:4}
        store_data, 'failat', t0, abs(tmp.ilat)
        store_data, 'fadist', t0, tmp.dis
        store_data, 'faorbit', t0, tmp.orbit
        tvars = 'fa'+['mlt','ilat','dist']
        options, tvars, 'colors', 4
        options, tvars, 'labels', 'fast'
        options, tvars, 'labflag', 3
        options, 'faorbit', 'ytitle', 'Fast orbit'
        options, 'faorbit', 'format', '(I05)'
        
        ; combine polar and fast pos.
        store_data, 'mlt', data = ['po','fa']+'mlt'
        store_data, 'ilat', data = [['po','fa']+'ilat',['po','fa']+'hem']
        store_data, 'dist', data = ['po','fa']+'dist'
        ylim, 'mlt', 0,24, 0
        ylim, 'ilat', 60, 90, 0
        ylim, 'dist', 1,10, 0
        options, 'mlt', 'ytitle', 'MLT (hr)'
        options, 'ilat', 'ytitle', 'Ilat (deg)'
        options, 'dist', 'ytitle', 'R (Re)'
        options, 'pomlt', 'labpos', 16
        options, 'famlt', 'labpos', 8
        options, 'poilat', 'labpos', 80
        options, 'failat', 'labpos', 70
        options, 'podist', 'labpos', 7
        options, 'fadist', 'labpos', 4
        
        ; plot.
        vars = ['bgse','pdyn','symh','poJEi','poJEe','faJEi','faJEe',$
            'ilat','mlt','dist']
        labs = ['flow','ae','faorbit']
        erase
        sgindexcolor, 43
        tplot, vars, var_label = labs, trange = utr, $
            title = 'Cusp, omni+polar+fast. '+time_string(t0[0],tformat='YYYY-MM-DD')
        print, 'use tlimit to select time, type .c to continue'
        tlimit

        ; select polar cusp time.
        flag = 'l'
        while strlowcase(flag) eq 'l' do begin
            ctime, cusputr, npoints = 2
            cuspetr = stoepoch(cusputr,'unix') & midutr = mean(cusputr)
            printf, -1, 'Polar cusp time: ', sfmepoch(cuspetr)
            read, flag, prompt = '(L)oop, (C)ancel, (S)ave: '
            if flag ne 'c' and n_elements(cuspetr) ne 2 then flag = 'l'
        endwhile

        if strlowcase(flag) eq 'c' then continue
        
        ; prepare output.
        str = '' & tab = '  '
        str += sfmepoch(cuspetr[0],'yyyy-MMdd')+tab
        str += sfmepoch(cuspetr[0],'hh:mm')+'-'+sfmepoch(cuspetr[1],'hh:mm')+tab
        ; dist.
        get_data, 'podist', t0, dat & tmp = sinterpol(dat, t0, cusputr)
        str += string(tmp[0], format='(F3.1)')+'-'+string(tmp[1], format='(F3.1)')+tab
        ; ilat.
        get_data, 'pohem', t0, dat & tmp = sinterpol(dat, t0, cusputr[0])
        tmp = (tmp ge 0)? 'N':'S' & str += tmp
        get_data, 'poilat', t0, dat & tmp = sinterpol(dat, t0, cusputr)
        str += string(tmp[0], format='(F4.1)')+'-'+string(tmp[1], format='(F4.1)')+tab
        ; mlt.
        get_data, 'pomlt', t0, dat & tmp = sinterpol(dat, t0, cusputr)
        str += string(tmp[0], format='(F4.1)')+'-'+string(tmp[1], format='(F4.1)')+tab
        ; dst.
        get_data, 'symh', t0, dat & idx = where(t0 ge cusputr[0] and t0 le cusputr[1])
        str += string(max(dat[idx]), format='(I4)')+tab
        ; ae.
        get_data, 'ae', t0, dat & idx = where(t0 ge cusputr[0] and t0 le cusputr[1])
        str += string(max(dat[idx]), format='(I4)')+tab
        
        printf, -1, str
        openw, lun, logfile, /append, /get_lun
        printf, lun, str
        free_lun, lun
        ylim, 'mlt', 06, 18, 0
        fn = plotdir+'cusp_'+sfmepoch(cuspetr[0],'yyyy_MMdd_hh')+'.eps'
        sgpsopen, fn
        tplot, vars, var_label = labs, $
            title = 'Cusp, omni+polar+fast. '+time_string(t0[0],tformat='YYYY-MM-DD')
        sgpsclose, /pdf
    endfor
end

gen_cusp_list_of_polar, stoepoch(['1999-08-28','1999-08-29']), datadir = sdiskdir('Research')
end