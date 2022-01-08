;+
; Type: crib.
; Purpose: Plot HOPE L3 2d pitch distribution at certain time (polygon).
; Parameters:
;   t0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   type0, in, string, opt. 'proton','oxygen','helium','electron'.
;       Default is 'proton'.
; Keywords:
;   probe, in, string, opt. 'a','b'. Default is 'a'
;   unit, in, string, opt. 'velocity', 'energy'.
;   log, in, boolean, opt. Default is linear.
;   hopel3, in/out, struct, opt. Pass among plotting routines to fast plot.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-12, Sheng Tian, create.
;-

pro plot_hope_l3_pitch2d_polygon, t0, type0, probe = probe, unit = unit0, $
    log = log, hopel3 = hopel3

    ; date and time.
    ut0 = time_double(t0)
    date = time_string(ut0,tformat='YYYY-MM-DD')
    utr = time_double(date)+[0,86400d]
    
    if n_elements(probe) eq 0 then probe = ['a']
    if n_elements(type0) eq 0 then type = 'proton' else $
        type = strlowcase(type0)
    if n_elements(unit0) eq 0 then unit = 'velocity' else $
        unit = strlowcase(unit0)

    case type of
        'electron':zrng = [4,8]
        'proton':zrng = [3.5,6]
        'oxygen':zrng = [3,6]
        'helium':zrng = [2,5]
    endcase
    
    print, 'rbsp'+probe+', '+type0+', '+time_string(ut0)
    
    ; suppress math excepttype.
    !except = 0
    
    ; constants.
    npixel = 5
    
    rad = !dpi/180d
    deg = 180d/!dpi
        
    case type of
        'electron': type0 = 'FEDU'
        'proton': type0 = 'FPDU'
        'oxygen': type0 = 'FODU'
        'helium': type0 = 'FHEDU'
    endcase
    
    ; load hope l3 data.
    load = 0
    if n_elements(hopel3) eq 0 then load = 1 else begin
        get_data, 'hopel3', t0, info
        if info.probe ne probe then load = 1
    endelse

    if load eq 1 then begin
        timespan, date, 1, /day
        hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
        store_data, 'hopel3', t0, {probe:probe}
    endif
    
    case type of
        'electron': mass0 = 1d/1836
        'proton': mass0 = 1d
        'oxygen': mass0 = 16d
        'helium': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
    
    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel3) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel3) eq enidx)
    tmp = min(hopel3.(epidx)-stoepoch(ut0,'unix'),rec,/absolute)
    ut0 = sfmepoch(hopel3.(epidx)[rec],'unix')  ; time in the middle of a spin.
    i0 = where(tag_names(hopel3) eq type0)
    dtidx = (type eq 'electron')? 'EPOCH_ELE_DELTA': 'EPOCH_ION_DELTA'
    dtidx = where(tag_names(hopel3) eq dtidx)
    type_dt = 2*hopel3.(dtidx)[rec]*1e-3
    
    ; level 3.
    datl3 = reform((hopel3.(i0))[rec,*,*])            ; in [nen,npa].
    enl3s = reform(hopel3.(enidx)[rec,*])
    pal3s = reform(hopel3.pitch_angle)
    npal3 = n_elements(pal3s)
    idx = where(datl3 eq -1e31, cnt)
    if cnt ne 0 then datl3[idx] = !values.d_nan
    
    ; remove duplicated energy bins.
;    enl3s = enl3s[0:*:2]
    idx = uniq(enl3s,sort(enl3s))
    enl3s = enl3s[idx]
    nenl3 = n_elements(enl3s)
    datl3 = datl3[idx,*]
;    datl3 = datl3[0:*:2,*]    
    
    

; **** the data for polar contour.
    tdat = transpose([[datl3]])     ; in [npa,nen].
    tang = [pal3s]
    case unit of
        'energy': begin
            tdis = enl3s
            xtitl = 'E (eV)'
            end
        'velocity': begin
            tdis = sqrt(2*enl3s/mass0)*1e-3
            xtitl = 'V (km/s)'
            end
    endcase
    if keyword_set(log) then begin
        tdis = alog10(tdis)
        xtitl = 'Log!D10!N '+xtitl
    endif
        
;    tang = tang # ((bytarr(nenl3)+1)+smkarthm(0,0.001,nenl3,'n'))
;    tdis = tdis ## (bytarr(2*npal3)+1)
    
    tang = tang*rad
    ndis = n_elements(tdis)
    cangs = [0,9,27,45,63,81,99,117,135,153,171,180]*!dpi/180
    cdiss = 0.5*(tdis[0:ndis-2]+tdis[1:ndis-1])
    cdiss = [cdiss[0]*2-cdiss[1],cdiss,cdiss[ndis-2]*2-cdiss[ndis-3]]


    ; remove nan.
    idx = where(finite(tdat,/nan))
    tdat[idx] = 0
    

    idx = where(tdat ne 0)
    min0 = min(tdat[idx],/nan)
    max0 = max(tdat[idx],/nan)
    nztick = 10
    if n_elements(zrng) eq 0 then $
        zrng = [floor(alog10(min0)),ceil(alog10(max0))-2]
    tpos = [0.15,0.15,0.85,0.85]
    titl = 'RBSP-'+strupcase(probe)+' HOPE L3 '+type+' Eflux!C'+$
        time_string(ut0-0.5*type_dt)+' - '+time_string(ut0+0.5*type_dt,tformat='hh:mm:ss')
    
    rootdir = shomedir()+'/rbsp_hope'
    ofn = rootdir+'/rbsp'+probe+'/'+type+'/hope_l3_pitch2d_'+type+'_'+ $
        time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.pdf'
    ;ofn = 0
    sgopen, ofn, xsize = 5, ysize = 5, /inch
;    sgindexcolor, 43, file = 'ct2'
    sgtruecolor
    sgdistr2d_polygon, tdat, tang, tdis, cangs, cdiss, position = tpos, $
        zrange = zrng, title = titl, xtitle = xtitl, ncolor = 10
    sgclose
    
;    ofn = shomedir()+'/hope_l3_pitch_'+type+'_'+time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.eps'
;    sgpsopen, ofn, xsize = 5, ysize = 4, /inch
;    ;sgwopen, 0, xsize = 5, ysize = 4, /inch
;    sgindexcolor, 43
;    
;    clrs = smkarthm(10,250,nenl3,'n')
;    xrng = [0,360]
;    yrng = [1e3,1e7]
;    zrng = minmax(enl3s)
;    zttl = 'Energy (eV)'
;    tpos = [0.15,0.15,0.85,0.85]
;    
;    pos1 = tpos
;    plot, xrng, yrng, position = pos1, /nodata, $
;        color = 0, background = 255, $
;        yrange = yrng, ystyle = 1, ylog = 1, ytitle = 'Flux (s!E-1!Ncm!E-2!Nsr!E-1!NeV!E-1!N)', $
;        xrange = xrng, xstyle = 1, xlog = 0, xtitle = 'Sector Angle (Deg)'
;    for i = 0, nenl3-1 do begin
;        idx = sort(pal3s)
;        txs = [pal3s[idx],360-pal3s[idx]]
;        tys = reform(datl3[i,idx])
;        tys = [tys,tys]
;        oplot, txs, tys, color = clrs[i], psym = -4, symsize = 0.2
;    endfor
;    
;    emsz = double(!d.x_ch_size)/!d.x_size
;    exsz = double(!d.y_ch_size)/!d.y_size
;    pos2 = tpos
;    pos2[0] = pos2[2]+emsz & pos2[2] = pos2[0]+emsz
;    sgcolorbar, clrs, position = pos2, $
;        zrange = zrng, ztitle = zttl, zcharsize = 0.8
;        
;    sgpsclose, /pdf
end

utr = time_double(['2013-05-01/07:35','2013-05-01/07:50'])
uts = smkarthm(utr[0],utr[1],12,'dx')
;uts = time_double(['2013-05-01/07:36:43'])
log = 0
unit = 'velocity'
types = ['electron','proton','oxygen','helium']
probes = ['b']

;uts = time_double(['2013-01-17/02:07:53','2013-07-07/13:26:07','2014-04-14/07:48:15'])
;log = 1
;unit = 'velocity'
;types = ['electron','proton','oxygen','helium']
;
;utr = time_double(['2013-01-17/02:05','2013-01-17/02:10'])
;uts = smkarthm(utr[0],utr[1],12,'dx')
;log = 1
;unit = 'energy'
;types = ['electron','proton','oxygen','helium']
;
;utr = time_double(['2013-07-07/13:24','2013-07-07/13:29'])
;uts = smkarthm(utr[0],utr[1],12,'dx')
;log = 1
;unit = 'energy'
;types = ['electron','proton','oxygen','helium']
;
;utr = time_double(['2014-04-14/07:46','2014-04-14/07:51'])
;uts = smkarthm(utr[0],utr[1],12,'dx')
;log = 1
;unit = 'energy'
;types = ['electron','proton','oxygen','helium']


;;utr = time_double(['2015-12-14/13:00','2015-12-14/13:00:30'])
;utr = time_double(['2015-12-14/13:00','2015-12-14/14:00:00'])
;uts = smkarthm(utr[0],utr[1],12,'dx')
;log = 1
;unit = 'energy'
;types = ['electron','proton','oxygen','helium']
;probes = ['a','b']

for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l3_pitch2d_polygon, uts[i], types[j], unit = unit, $
                log = log, hopel3 = hopel3, probe = probes[k]
end
