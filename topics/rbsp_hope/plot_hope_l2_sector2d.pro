;+
; Type: crib.
; Purpose: Plot HOPE L2 2d sector distribution at certain time.
; Parameters:
;   date, in, string, req. 'YYYY-MM-DD'.
;   ut0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   probe, in, string, req. 'a','b'.
;   ion, in, string, req. 'proton','oxygen','helium'.
; Keywords: none.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-12, Sheng Tian, create.
;-

pro plot_hope_l2_sector2d, t0, type0, probe = probe, unit = unit0, log = log, $
    hopel2 = hopel2, datb = datb

    ; date, time, and settings.
    ut0 = time_double(t0)
    date = time_string(ut0,tformat='YYYY-MM-DD')
    utr = time_double(date)+[0,86400d]
    
    if n_elements(probe) eq 0 then probe = ['b']
    if n_elements(type0) eq 0 then type = 'proton' else type = strlowcase(type0)
    if n_elements(unit0) eq 0 then unit = 'velocity' else unit = strlowcase(unit0)

    case type of
        'electron':zrng = [4,8]
        'proton':zrng = [3.5,6]
        'oxygen':zrng = [3,5]
        'helium':zrng = [2,5]
    endcase
    
    print, 'rbsp'+probe+', '+type0+', '+time_string(ut0)

    ; suppress math excepttype.
    !except = 0
    
    ; constants.
    npxl = 5
    
    rad = !dpi/180d
    deg = 180d/!dpi
        
    case type of
        'electron': type0 = 'FEDU'
        'proton': type0 = 'FPDU'
        'oxygen': type0 = 'FODU'
        'helium': type0 = 'FHEDU'
    endcase

    case type of
        'electron': mass0 = 1d/1836
        'proton': mass0 = 1d
        'oxygen': mass0 = 16d
        'helium': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.

    
    ; load hope l2 data.
    if n_elements(hopel2) eq 0 then begin
        timespan, date, 1, /day
        hopel2 = sread_rbsp_hope_l2(utr, probes = probe)
    endif

    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel2) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel2) eq enidx)
    tmp = min(hopel2.(epidx)-stoepoch(ut0,'unix'),rec,/absolute)
    ut0 = sfmepoch(hopel2.(epidx)[rec],'unix')
    i0 = where(tag_names(hopel2) eq type0)
    type_dt = 12

    ; level 2.
    datl2 = reform((hopel2.(i0))[rec,*,*,*])            ; in [72,16,5].
    enl2s = reform(hopel2.(enidx)[rec,*])               ; energy bins.
    enmode = hopel2.energy_collapsed[rec]               ; 0 for collaped bins.
    idx = where(datl2 eq -1e31, cnt)
    if cnt ne 0 then datl2[idx] = !values.d_nan

    ; remove duplicated energy bins.
;    enl2s = enl2s[0:*:2]
;    nenl2 = n_elements(enl2s)
;    datl2 = datl2[0:*:2,*,*]    
    idx = uniq(enl2s,sort(enl2s))
    enl2s = enl2s[idx]
    nenl2 = n_elements(enl2s)
    datl2 = datl2[idx,*,*]


    ; remove duplicated sectors. datl2 > dat2a.
    scmode = reform(hopel2.sector_collapse_cntr[rec,*]) ; sector collapse mode.
    ; scmode = intarr(npxl)+1                           ; no collapsing.
    ; adjust value for collapsing.
    for j = 0, npxl-1 do begin
        if scmode[j] eq 1 then continue
        ;    datl2[*,*,j] = (datl2[*,*,j]-16384)*64+16384
        ;    datl2[*,*,j] = datl2[*,*,j]/scmode[j]
    endfor

    secs = 16/scmode
    nsa2a = total(secs)
    sa2as = dblarr(nsa2a)       ; raw sector angle, have duplicated angles.
    for j = 0, npxl-1 do begin
        nsec = secs[j]
        ang0 = 11.25*scmode[j]  ; offset angle accounts for collapsing.
        sa2as[total(secs[0:j])-secs[j]] = smkarthm(ang0,360d/nsec,nsec,'x0')
    endfor

    ; dat2a.
    dat2a = dblarr(nsa2a,nenl2)
    for i = 0, nenl2-1 do begin
        tmp = []
        for j = 0, npxl-1 do tmp = [tmp,reform(datl2[i,0:*:scmode[j],j])]
        dat2a[*,i] = tmp
    endfor

    ; massage into sector angle, i.e., averaging for same sector angle.
    sa2bs = suniq(sa2as)        ; uniq sector angle.
    nsa2b = n_elements(sa2bs)
    en2bs = enl2s
    nen2b = n_elements(en2bs)
    dat2b = dblarr(nsa2b,nen2b)     ; in [nsa,nen].
    for i = 0, nsa2b-1 do begin
        idx = where(sa2as eq sa2bs[i])
        for j = 0, nen2b-1 do dat2b[i,j] = mean(dat2a[idx,j])
    endfor



; **** the data for polar contour.
    tdat = dat2b                    ; in [nsa,nen].
    tang = sa2bs
    case unit of
        'energy': begin
            tdis = enl2s
            xtitl = 'E (eV)'
            end
        'velocity': begin
            tdis = sqrt(2*enl2s/mass0)*1e-3
            xtitl = 'V (km/s)'
            end
    endcase
    if keyword_set(log) then begin
        tdis = alog10(tdis)
        xtitl = 'Log!D10!N '+xtitl
    endif

    tang = tang # ((bytarr(nen2b)+1)+smkarthm(0,0.001,nen2b,'n'))
    tdis = tdis ## (bytarr(nsa2b)+1)

    tang = tang*rad

    ; remove nan.
    idx = where(finite(tdat,/nan))
    tdat[idx] = 0


    idx = where(tdat ne 0)
    min0 = min(tdat[idx],/nan)
    max0 = max(tdat[idx],/nan)
    nztick = 10
    if n_elements(zrng) eq 0 then zrng = [floor(alog10(min0)),ceil(alog10(max0))-2]
    tpos = [0.15,0.15,0.85,0.85]
    titl = 'RBSP-'+strupcase(probe)+' HOPE L2 Sector Angle!C'+type+' Eflux!C'+$
        time_string(ut0)+' - '+time_string(ut0+type_dt,tformat='hh:mm:ss')
    
    rootdir = shomedir()+'/rbsp_hope'
    ofn = rootdir+'/rbsp'+probe+'/'+type+'/hope_l2_sector2d_'+type+'_'+time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.pdf'
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    sgindexcolor, 43, file = 'ct2'
    sgdistr2d, tdat, tang, tdis, position = tpos, zrange = zrng, title = titl, xtitle = xtitl, ncolor = 10

    ; plot B field.
    ; load b_uvw.
    if n_elements(datb) eq 0 then begin
        emfisis = sread_rbsp_emfisis(utr)
        vars = emfisis.name
        b_uts = sfmepoch(*emfisis[where(vars eq 'Epoch')].value,'unix',/tt2000)
        b_gse = *emfisis[where(vars eq 'Mag')].value
        store_data, 'b_gse', b_uts, b_gse
        rbsp_uvw2gse, 'b_gse', newname = 'b_uvw', probe = prob, /inverse
        get_data, 'b_uvw', b_uts, b_uvw
        datb = {b_uts:b_uts,b_uvw:b_uvw,b_gse:b_gse}
    endif
    b_uts = datb.b_uts
    b_uvw = datb.b_uvw
    tbuvw = reform(sinterpol(b_uvw,b_uts,ut0-type_dt, /spline))
    anga = atan(tbuvw[1],tbuvw[0])
    angb = atan(tbuvw[2],sqrt(tbuvw[0]^2+tbuvw[1]^2))

    plot, [-1,1],[-1,1], /noerase, /nodata, position = tpos, $
        xstyle = 5, ystyle = 5
    arrow, 0,0, cos(angb)*cos(anga),cos(angb)*sin(anga), /data

    sgclose


; ;sgpsopen, shomedir()+'/hope_l2_sector_'+ion+'.eps', xsize = 5, ysize = 4, /inch
; sgopen, 0, xsize = 5, ysize = 4, /inch
; sgindexcolor, 43, file = 'ct2'
; 
; clrs = smkarthm(10,250,nen2b,'n')
; xrng = [0,360]
; yrng = [1e3,1e7]
; zrng = minmax(en2bs)
; zttl = 'Energy (eV)'
; tpos = [0.15,0.15,0.85,0.85]
; 
; pos1 = tpos
; plot, xrng, yrng, position = pos1, /nodata, $
;     color = 0, background = 255, $
;     yrange = yrng, ystyle = 1, ylog = 1, ytitle = 'Flux (s!E-1!Ncm!E-2!Nsr!E-1!NeV!E-1!N)', $
;     xrange = xrng, xstyle = 1, xlog = 0, xtitle = 'Sector Angle (Deg)'
; for i = 0, nen2b-1 do begin
;     idx = sort(sa2bs)
;     oplot, sa2bs[idx], dat2b[idx,i], color = clrs[i], psym = -4, symsize = 0.2
; endfor
; 
; emsz = double(!d.x_ch_size)/!d.x_size
; exsz = double(!d.y_ch_size)/!d.y_size
; pos2 = tpos
; pos2[0] = pos2[2]+emsz & pos2[2] = pos2[0]+emsz
; sgcolorbar, clrs, position = pos2, $
;     zrange = zrng, ztitle = zttl, zcharsize = 0.8
; 
; ;sgclose

end



utr = time_double(['2013-05-01/07:37','2013-05-01/07:42'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 0
unit = 'velocity'
types = ['electron','proton','oxygen']
probes = ['b']

for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l2_sector2d, uts[i], types[j], unit = unit, log = log, probe = probes[k], $
            hopel2 = hopel2, datb = datb
end

