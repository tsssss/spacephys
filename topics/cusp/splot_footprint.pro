; Given altitude in km, NOT distance in Re.
; Time resolution is 1 sec. (Too high?, change to 20 sec?)

function fa_sdt_plot_footpoint, fn, trange
    compile_opt idl2
    on_error, 0

    if file_test(fn) eq 0 then $
        message, 'file does not exist ...'

    if n_elements(trange) eq 0 then $
        message, 'time range is undefined ...'  

    if n_elements(nm) eq 0 then nm = 'fa_'

    ; =====================
    ; read mlt, ilat, dist.
    ; =====================
    tplot_restore, filename = fn

    ; time, unix time.
    get_data, 'ALT', data = tmp
    t0 = tmp.x

    ; ilat, mlt, alt, they all share the same epoch.
    alt = tmp.y
    get_data, 'ILAT', data = tmp
    ilat = tmp.y
    get_data, 'MLT', data = tmp
    mlt = tmp.y

    ; trim to given time range.
    idx = where(t0 ge trange[0] and t0 le trange[1])
    if idx[0] eq -1 then $
        message, 'wrong time range ...'

    t0   = t0[idx]
    ilat = ilat[idx]
    mlt  = mlt[idx]
    alt  = alt[idx]

    ; special care.
    ; change altitude in km to distance in Re.
    re1 = 1D/6356.7523      ; polar radius = 6356.7523 km.
    dis = alt*re1+1

    ; ===========
    ; store data.
    ; ===========
    ; need t0, ilat, mlt, dis.
    return, {name:nm, t0:t0, ilat:ilat, mlt:mlt, dis:dis}

end


; time resolution is 1 min.
; ilat doesn't have sign.

function po_sdt_plot_footpoint, fn, trange
    compile_opt idl2
    on_error, 0

    if file_test(fn) eq 0 then $
        message, 'file does not exist ...'

    if n_elements(trange) eq 0 then $
        message, 'time range is undefined ...'  

    if n_elements(nm) eq 0 then nm = 'po_'

    ; =====================
    ; read mlt, ilat, dist.
    ; =====================
    sdt = ssdtread(fn)

    ; time, unix time.
    t0 = sdt.var.polarinvariantlatitude.depend_0

    ; ilat, mlt, dis, they all share the same epoch.
    ilat = sdt.var.polarinvariantlatitude.value
    mlt  = sdt.var.polarmlt.value
    dis  = sdt.var.polarspcraftdist.value
    mlat = sdt.var.polarmaglat.value

    ; trim to given time range.
    idx = where(t0 ge trange[0] and t0 le trange[1])
    if idx[0] eq -1 then $
        message, 'wrong time range ...'

    t0   = t0[idx]
    ilat = ilat[idx]
    mlat = mlat[idx]
    mlt  = mlt[idx]
    dis  = dis[idx]

    ; special care.
    ; give ilat sign: positive is north, negative is south.
    idx = where(mlat lt 0)
    if idx[0] ne -1 then ilat[idx] *= -1

    ; ===========
    ; store data.
    ; ===========
    ; need t0, ilat, mlt, dis.
    return, {name:nm, t0:t0, ilat:ilat, mlt:mlt, dis:dis}

end

pro splot_footprint, fafn, fatr, pofn, potr, ps = ps

    coef = 1D/60

    fa = fa_sdt_plot_footpoint(fafn, fatr)
    po = po_sdt_plot_footpoint(pofn, potr)

    mn = min([min(fa.t0),min(po.t0)])
    mx = max([max(fa.t0),max(po.t0)])
    rg = [mn,mx]
    rgt0 = rg
    fa.t0 -= rgt0[0]
    po.t0 -= rgt0[0]
    rgt0 -= rgt0[0]
    rgt0 *= coef

    mn = min([min(fa.mlt),min(po.mlt)])
    mx = max([max(fa.mlt),max(po.mlt)])
    rg = [mn,mx]        ; expand by 1.2.
    rgmlt = [rg[1]*1.1-rg[0]*0.1,rg[0]*1.1-rg[1]*0.1]
    rgmlt = [9,15]

    mn = min([min(fa.ilat),min(po.ilat)])
    mx = max([max(fa.ilat),max(po.ilat)])
    rg = [mn,mx]        ; expand by 1.2.
    rgilat = [rg[1]*1.1-rg[0]*0.1,rg[0]*1.1-rg[1]*0.1]
    rgilat = [50,90]*((fa.ilat[0] gt 0)*2-1)

    mn = min([min(fa.dis),min(po.dis)])
    mx = max([max(fa.dis),max(po.dis)])
    rg = [mn,mx]        ; expand by 1.2.
    rgdis = [rg[0]*1.1-rg[1]*0.1,rg[1]*1.1-rg[0]*0.1]

    ; plot.
    if keyword_set(ps) then begin
        thedev = !d.name
        set_plot, 'ps'
        fn = time_string(fatr[0], tformat = 'YYYYMMDD')+'_footprint.ps'
        device, filename = fn
        ;      !p.charsize = 1.5
        ;      !p.charthick = 2.0
    endif else begin
        ;      !p.charsize = 2.0
        ;      !p.charthick = 1.0
    endelse
    !p.multi = [0,1,3,0,0]
    title = 'fast and polar footprint, ' + $
        time_string(min([fatr[0], potr[0]]))+', po: --, fa: ..'
    fa.t0 *= coef
    po.t0 *= coef
    plot, fa.t0, fa.mlt, xrange = rgt0, yrange = rgmlt, $
        xstyle = 1, ystyle = 1, ytitle = 'mlt (hr)', /nodata, $
        title = title
    oplot, fa.t0, fa.mlt, linestyle = 1
    oplot, po.t0, po.mlt, linestyle = 2

    plot, fa.t0, fa.ilat, xrange = rgt0, yrange = rgilat, $
        xstyle = 1, ystyle = 1, ytitle = 'ilat (deg)', /nodata
    oplot, fa.t0, fa.ilat, linestyle = 1
    oplot, po.t0, po.ilat, linestyle = 2

    plot, fa.t0, fa.dis, xrange = rgt0, yrange = rgdis, $
        xstyle = 1, ystyle = 1, ytitle = 'dis (Re)', /nodata, $
        xtitle = 'time (min)'
    oplot, fa.t0, fa.dis, linestyle = 1
    oplot, po.t0, po.dis, linestyle = 2

    if keyword_set(ps) then begin
        device, /close
        set_plot, thedev
    endif
end

workroot = shomedir() + '/Works'
fafn = workroot + '/data/event/fa_sdt_fld_19970527_v01.tplot'
fatr = time_double(['1997-05-27/18:15', '1997-05-27/18:21'])
pofn = workroot + '/data/event/po_sdt_fld_19970527_v02.sdt'
potr = time_double(['1997-05-27/18:05', '1997-05-27/18:15'])

fafn = workroot + '/data/event/fa_sdt_fld_19970914_v01.tplot'
fatr = time_double(['1997-09-14/13:31', '1997-09-14/13:35'])
pofn = workroot + '/data/event/po_sdt_fld_19970914_v02.sdt'
potr = time_double(['1997-09-14/13:35', '1997-09-14/13:45'])

fafn = workroot + '/data/event/fa_sdt_fld_19970817_v01.tplot'
fatr = time_double(['1997-08-17/13:25', '1997-08-17/13:28'])
pofn = workroot + '/data/event/po_sdt_fld_19970817_v02.sdt'
potr = time_double(['1997-08-17/13:15', '1997-08-17/13:30'])

fafn = workroot + '/data/event/fa_sdt_fld_19970914_v01.tplot'
fatr = time_double(['1997-09-14/13:31', '1997-09-14/13:34'])
pofn = workroot + '/data/event/po_sdt_fld_19970914_v02.sdt'
potr = time_double(['1997-09-14/13:35', '1997-09-14/13:45'])

fafn = workroot + '/data/event/fa_sdt_fld_19981002_v01.tplot'
fatr = time_double(['1998-10-02/15:35', '1998-10-02/15:41'])
pofn = workroot + '/data/event/po_sdt_fld_19981002_v02.sdt'
potr = time_double(['1998-10-02/15:50', '1998-10-02/16:20'])

splot_footprint, fafn, fatr, pofn, potr, /ps
end
