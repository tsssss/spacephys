
pro plot_polar_fast_summary, tr, no_delete = no_delete, title = titl, $
    nofast = nofast
    
    compile_opt idl2
    
    ; plot options.
    if ~keyword_set(no_delete) then store_data, tnames(), /delete
    dev0 = !d.name
    device, decomposed = 0
    loadct2, 43

    utr = sfmepoch(tr,'unix')
    
    ; read polar hydra.
    tmp = sread_polar_hydra(utr)
    t0 = sfmepoch(tmp.epoch, 'unix')
    store_data, 'poJEi', data = {x:t0, y:tmp.jei, v:tmp.eni}
    store_data, 'poJEe', data = {x:t0, y:tmp.jee, v:tmp.ene}
    options, 'po'+['JEi','JEe'], 'ztitle', '1/(cm!U2!N-s-sr)'
    options, 'po'+['JEi','JEe'], 'spec', 1
    options, 'po'+['JEi','JEe'], 'no_interp', 1
    options, 'poJEi', 'ytitle', 'Polar ion!C(eV)'
    options, 'poJEe', 'ytitle', 'Polar ele!C(eV)'
    ylim, 'poJEi', tmp.eni[0], tmp.eni[-1], 1
    ylim, 'poJEe', tmp.ene[0], tmp.ene[-1], 1
    zlim, 'poJEi', 1e4, 1e8, 1
    zlim, 'poJEe', 1e4, 5e8, 1

    ; read fast esa.
    if ~keyword_set(nofast) then begin
        tmp = sread_fast_esa(utr, type = 'ies')
        t0 = sfmepoch(tmp.epoch, 'unix')
        store_data, 'faJEi', data = {x:t0, $
            y:(tmp.ion_0+tmp.ion_90+tmp.ion_180)/3, v:tmp.ion_en}
        tmp = sread_fast_esa(utr, type = 'ees')
        t0 = sfmepoch(tmp.epoch, 'unix')
        store_data, 'faJEe', data = {x:t0, $
            y:(tmp.el_0+tmp.el_90+tmp.el_180)/3, v:tmp.el_en}
        options, 'fa'+['JEi','JEe'], 'ztitle', '1/(cm!U2!N-s-sr)'
        options, 'fa'+['JEi','JEe'], 'spec', 1
        options, 'fa'+['JEi','JEe'], 'no_interp', 1
        options, 'faJEi', 'ytitle', 'Fast ion!C(eV)'
        options, 'faJEe', 'ytitle', 'FAST ele!C(eV)'
        ylim, 'faJEi', 4, 40000, 1
        ylim, 'faJEe', 4, 40000, 1
        zlim, 'faJEi', 1e4, 1e8, 1
        zlim, 'faJEe', 1e6, 1e10, 1
    endif
    
    ; read polar pos.
    re = 6378d
    vars = ['Epoch','GSE_POS','EDMLT_TIME','L_SHELL','MAG_LATITUDE']
    tmp = sread_polar_orbit(utr, vars = vars)
    t0 = sfmepoch(tmp.epoch, 'unix')
    ilat = acos(sqrt(1d/tmp.l_shell))*(180d/!dpi)
    idx = where(tmp.mag_latitude lt 0, cnt)
    if cnt gt 0 then ilat[idx] = -ilat[idx]
    dis = sqrt(total(tmp.gse_pos^2,2))*(1d/re)
    store_data, 'pomlt', t0, tmp.edmlt_time
    store_data, 'pohem', t0, ((ilat ge 0)*2-1)*(90-1), $
        limits = {thick:4, colors:2}
    store_data, 'poilat', t0, abs(ilat)
    store_data, 'podist', t0, dis
    tvars = 'po'+['mlt','ilat','dist']
    options, tvars, 'colors', 2
    options, tvars, 'labels', 'polar'
    options, tvars, 'labflag', 3
    
    ; read fast pos.
    tmp = sread_fast_pos(utr, $
        vars = ['Epoch','mlt','ilat','dis','orbit'])
    t0 = sfmepoch(tmp.epoch, 'unix')
    store_data, 'famlt', t0, tmp.mlt
    store_data, 'fahem', t0, ((tmp.ilat ge 0)*2-1)*(60+1), $
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
    options, 'mlt', 'ytitle', 'MLT!C(hr)'
    options, 'ilat', 'ytitle', 'Ilat!C(deg)'
    options, 'dist', 'ytitle', 'R!C(Re)'
    options, 'pomlt', 'labpos', 16
    options, 'famlt', 'labpos', 8
    options, 'poilat', 'labpos', 80
    options, 'failat', 'labpos', 70
    options, 'podist', 'labpos', 7
    options, 'fadist', 'labpos', 4    
    
    ; read dst.
    tmp = sread_omni(utr)
    t0 = sfmepoch(tmp.epoch, 'unix')
    store_data, 'dst', t0, tmp.sym_h
    store_data, 'ae', t0, tmp.ae_index
    store_data, 'symh', t0, tmp.sym_h
    store_data, 'pdyn', t0, tmp.pressure
    store_data, 'bgse', t0, [[[tmp.bx_gse]],[[tmp.by_gse]],[[tmp.bz_gse]]]
    store_data, 'flow', t0, tmp.flow_speed
    store_data, 'vgse', t0, [[[tmp.vy]],[[tmp.vz]]]
    options, 'ae', 'ytitle', 'AE (nT)'
    options, 'ae', 'format', '(I5)'
    options, 'symh', 'ytitle', 'SymH!c(nT)'
    options, 'pdyn', 'ytitle', 'P!C(nPa)'
    options, 'bgse', 'ytitle', 'B GSE!C(nT)'
    options, 'vgse', 'ytitle', 'V GSE!C(km/s)'
    options, 'flow', 'ytitle', 'V!DSW!N (km/s)'
    options, 'flow', 'format', '(F6.1)'
    options, 'bgse', 'labels', ['x','y','z']
    options, 'bgse', 'colors', [6,4,2]
    options, 'vgse', 'labels', ['y','z']
    options, 'vgse', 'colors', [4,2]

    ; prepare plot.
    dirptn = shomedir()+'/polar_fast/yyyy/mm/'
    fnptn = 'yyyymmddhh_polarfast.ps'
    dir = sptn2fn(dirptn, tr[0])
    if file_test(dir) eq 0 then file_mkdir, dir
    fn = dir+sptn2fn(fnptn, tr[0])
    fn = shomedir()+'/'+sptn2fn(fnptn, tr[0])
    tl = 'Cusp, omni+polar+fast. '+sfmepoch(tr[0],'YYYY-MM-DD')
    if n_elements(titl) ne 0 then tl = titl
    vars = ['bgse','pdyn','symh', 'poJEi','poJEe', 'faJEi','faJEe', $
        'ilat','mlt','dist']

    utr = sfmepoch(tr,'unix')
    lbls = ['flow','ae','faorbit']
    tplot, vars, var_label = lbls, trange = utr, title = tl
;    pstplot, filename = fn, chsize = 0.8
    set_plot, dev0

;    store_data, tnames(), /delete
end

tr = stoepoch(['2003-11-20/15:30','2003-11-20/16:50'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2001-04-01/16:45','2001-04-01/17:35'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2003-10-29/10:15','2003-10-29/12:20'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2004-11-07/19:00','2004-11-07/21:00'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2004-11-08/14:10','2004-11-08/15:40'], 'yyyy-MM-dd/hh:mm')  ; not plotted.
tr = stoepoch(['2001-11-03/11:00','2001-11-03/12:20'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2001-11-03/09:20','2001-11-03/10:40'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2001-11-04/03:30','2001-11-04/06:00'], 'yyyy-MM-dd/hh:mm')  ; nice n-s conjugate too.
tr = stoepoch(['2001-11-06/12:10','2001-11-06/13:50'], 'yyyy-MM-dd/hh:mm')
tr = stoepoch(['2001-03-22/19:10','2001-03-22/19:35'], 'yyyy-MM-dd/hh:mm')

; storms when FAST had E.
tr = stoepoch(['2001-04-01/12:30','2001-04-01/13:30'], 'yyyy-MM-dd/hh:mm')

; tmp.
tr = stoepoch(['2001-10-28/07:30','2001-10-28/10:30'])
;pstplot
plot_polar_fast_summary, tr

end