;+
; Type: crib.
; Purpose: Plot HOPE L3 2d pitch distributtype at certain time (polar_contour).
; Parameters:
;   date, in, string, req. 'YYYY-MM-DD'.
;   ut0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   probe, in, string, req. 'a','b'.
;   type, in, string, req. 'proton','oxygen','helium','electron'.
; Keywords: none.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-12, Sheng Tian, create.
;-

pro plot_hope_l3_pitch2d, t0, the_type, probe = probe, $
    unit = unit0, log = log, hopel3 = hopel3, zrange = zrng0

    ; date and time.
    ut0 = time_double(t0)
    date = time_string(ut0,tformat='YYYY-MM-DD')
    utr = time_double(date)+[0,86400d]

    if n_elements(probe) eq 0 then probe = ['a']
    if n_elements(the_type) eq 0 then type = 'proton' else type = strlowcase(the_type)
    if n_elements(unit0) eq 0 then unit = 'velocity' else unit = strlowcase(unit0)

    case type of
        'electron':zrng = [4,8]
        'proton':zrng = [3.5,6]
        'oxygen':zrng = [3,5]
        'helium':zrng = [2,5]
    endcase
    if n_elements(zrng0) ne 0 then zrng = zrng0

    print, 'rbsp'+probe+', '+the_type+', '+time_string(ut0)

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


    ; load hope l3 data.
    load = 0
    if n_elements(hopel3) eq 0 then load = 1 else begin
        get_data, 'hopel3', t0, info
        if size(info,/type) ne 8 then load = 1 else if info.probe ne probe then load = 1
    endelse

    if load eq 1 then begin
        timespan, date, 1, /day
        hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
        store_data, 'hopel3', t0, {probe:probe}
    endif


    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel3) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel3) eq enidx)
    tmp = min(hopel3.(epidx)-stoepoch(ut0,'unix'),rec,/absolute)
    ut0 = sfmepoch(hopel3.(epidx)[rec],'unix')  ; time in the middle of a spin.
    i0 = where(tag_names(hopel3) eq type0)
    dtidx = (type eq 'electron')? 'EPOCH_ELE_DELTA': 'EPOCH_ION_DELTA'
    dtidx = where(tag_names(hopel3) eq dtidx)
    type_dt = hopel3.(dtidx)[rec]*1e-3

    ; level 3.
    datl3 = reform((hopel3.(i0))[rec,*,*])            ; in [nen,npa].
    enl3s = reform(hopel3.(enidx)[rec,*])
    pal3s = reform(hopel3.pitch_angle)
    npal3 = n_elements(pal3s)
    idx = where(datl3 eq -1e31, cnt)
    if cnt ne 0 then datl3[idx] = !values.d_nan
    idx = where(datl3 eq 0, cnt)
    if cnt ne 0 then datl3[idx] = !values.d_nan

;    ; remove duplicated energy bins.
;    enl3s = enl3s[0:*:2]
    idx = uniq(enl3s,sort(enl3s))
    enl3s = enl3s[idx]
    nenl3 = n_elements(enl3s)
    datl3 = datl3[idx,*]
;    datl3 = datl3[0:*:2,*]


    ; the data for polar contour.
    tdat = transpose([[datl3],[datl3]])     ; in [2*npa,nen].
    tang = [pal3s,360-pal3s]
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

    tang = tang # ((bytarr(nenl3)+1)+smkarthm(0,0.001,nenl3,'n'))
    tdis = tdis ## (bytarr(2*npal3)+1)

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
    titl = 'RBSP-'+strupcase(probe)+' HOPE L3 Pitch Angle!C'+type+' Eflux!C'+$
        time_string(ut0-type_dt)+' - '+time_string(ut0+type_dt,tformat='hh:mm:ss')

    rootdir = shomedir()+'/rbsp_hope'
    ofn = rootdir+'/rbsp'+probe+'/'+type+'/hope_l3_pitch2d_'+type+'_'+ $
        time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.pdf'
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    sgindexcolor, 43, file = 'ct2'
    sgdistr2d, tdat, tang, tdis, position = tpos, zrange = zrng, $
        title = titl, xtitle = xtitl, ncolor = 10
    sgclose

end

utr = time_double(['2013-05-01/07:35','2013-05-01/07:50'])
uts = smkarthm(utr[0],utr[1],12,'dx')
;uts = time_double(['2013-05-01/07:36:43'])
log = 0
unit = 'velocity'
types = ['electron','proton','oxygen','helium']

;uts = time_double(['2013-01-17/02:07:53','2013-07-07/13:26:07','2014-04-14/07:48:15'])
;log = 1
;unit = 'velocity'
;types = ['electron','proton','oxygen','helium']

utr = time_double(['2013-01-17/02:05','2013-01-17/02:10'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']

utr = time_double(['2013-07-07/13:24','2013-07-07/13:29'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']

utr = time_double(['2014-04-14/07:46','2014-04-14/07:51'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']


;utr = time_double(['2015-12-14/13:00','2015-12-14/13:00:30'])
utr = time_double(['2015-12-14/13:00','2015-12-14/14:00:00'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']
probes = ['a']

utr = time_double(['2013-01-17/01:47','2013-01-17/01:53'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']
probes = ['a','b']

utr = time_double(['2013-01-17/03:01','2013-01-17/03:07'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']
probes = ['a','b']

utr = time_double(['2013-01-17/04:11','2013-01-17/04:17'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']
probes = ['a','b']


; Mark's events.
utr = time_double(['2015-12-14/12:18','2015-12-14/12:38'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['electron','proton','oxygen']
probes = ['a']


; 2013-06-07 events.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
utr = time_double(['2013-06-07/05:02','2013-06-07/05:20'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'velocity'
types = ['proton']
probes = ['a','b']



for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l3_pitch2d_polygon, uts[i], types[j], unit = unit, log = log, hopel3 = hopel3, probe = probes[k]


end
