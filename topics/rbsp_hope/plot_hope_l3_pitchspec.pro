;+
; Type: crib.
; Purpose: Plot HOPE L3 time-pitch spectrogram.
; Parameters:
;   t0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   type0, in, string, opt. 'proton','oxygen','helium','electron'.
;       Default is 'proton'.
; Keywords:
;   probe, in, string, opt. 'a','b'. Default is 'a'
;   hopel3, in/out, struct, opt. Pass among plotting routines to fast plot.
;   energy, in, float or fltarr[2], opt. Energy in eV. Plot single energy bin
;       or certain energy range.  Default is to sum all energy bins.
; Notes: none.
; Dependence: slib.
; History:
;   2016-04-08, Sheng Tian, create.
;-

pro plot_hope_l3_pitchspec, t0, type0, probe = probe, $
    hopel3 = hopel3, energy = en0

    ; date and time.
    utr = time_double(t0)
    date = time_string(utr[0],tformat='YYYY-MM-DD')
    
    if n_elements(probe) eq 0 then probe = 'a'
    if n_elements(type0) eq 0 then type = 'proton' else $
        type = strlowcase(type0)
    if n_elements(unit0) eq 0 then unit = 'velocity' else $
        unit = strlowcase(unit0)

    case type of
        'electron':zrng = [4,8]
        'proton':zrng = [3.5,6]
        'oxygen':zrng = [3,5]
        'helium':zrng = [2,5]
    endcase
        
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
        if max(info.utr-utr) ne 0 then load = 1
    endelse

    if load eq 1 then begin
        timespan, date, 1, /day
        hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
        store_data, 'hopel3', t0, {probe:probe,utr:utr}
    endif
    
    case type of
        'electron': mass0 = 1d/1836
        'proton': mass0 = 1d
        'oxygen': mass0 = 16d
        'helium': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)    ; E in eV, mass in kg.
    
    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel3) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel3) eq enidx)
    i0 = where(tag_names(hopel3) eq type0)
    
    ; level 3.
    uts = sfmepoch(hopel3.(epidx),'unix')
    datl3 = hopel3.(i0)                 ; in [nrec,nen,npa].
    enl3s = hopel3.(enidx)              ; in [nrec,nen].
    pal3s = hopel3.pitch_angle          ; in [npa].
    npal3 = n_elements(pal3s)
    idx = where(datl3 eq -1e31, cnt)
    if cnt ne 0 then datl3[idx] = !values.d_nan
    
    ; remove duplicated energy bins.
    enl3s = reform(enl3s[0,*])
;    enl3s = enl3s[0:*:2]
    nenl3 = n_elements(enl3s)
;    datl3 = datl3[idx,*]
;    datl3 = datl3[0:*:2,*]
    
    

; **** the data for spectrogram.
    case n_elements(en0) of
        0: enidx = findgen(nenl3)
        1: tmp = min(enl3s-en0, enidx, /absolute)
        2: enidx = where(enl3s ge min(en0) and enl3s le max(en0))
    endcase
    tdat = total(datl3[*,enidx,*],2, /nan)
    tang = pal3s

    case n_elements(en0) of
        0: ztitl = 'All Energy'
        1: ztitl = sgnum2str(en0)+ ' eV'
        2: ztitl = sgnum2str(min(en0))+'-'+sgnum2str(max(en0))+' eV'
    endcase

    ; remove nan.
    idx = where(finite(tdat,/nan))
    tdat[idx] = 0
    

    idx = where(tdat ne 0)
    min0 = min(tdat[idx],/nan)
    max0 = max(tdat[idx],/nan)
    nztick = 256
    if n_elements(zrng) eq 0 then $
        zrng = [floor(alog10(min0)),ceil(alog10(max0))-2]
    tpos = [0.15,0.15,0.85,0.85]
    titl = 'RBSP-'+strupcase(probe)+' HOPE L3 '+type+' Eflux!C'+$
        time_string(uts[0])+' - '+time_string(uts[1],tformat='hh:mm:ss')

    tvar = 'rb'+probe+'_paspec'
    store_data, tvar, uts, tdat, tang
    options, tvar, 'ytitle', 'Pitch (deg)'
    options, tvar, 'yrange', [0,180]
    options, tvar, 'ystyle', 1
    options, tvar, 'labels', ztitl
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    options, tvar, 'zlog', 1
    
    rootdir = shomedir()+'/rbsp_hope'
    ofn = rootdir+'/rbsp'+probe+'/'+type+'/hope_l3_paspec_'+type+'_'+ $
        time_string(uts[0],tformat='YYYY_MMDD_hh')+'.pdf'
    ofn = 0
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    sgindexcolor, 43, file = 'ct2'
    tplot, tvar, trange = utr
    sgclose
    
end

utr = time_double(['2013-05-01/07:00','2013-05-01/08:30'])
log = 0
types = ['proton']
probes = ['b']


for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
            plot_hope_l3_pitchspec, utr, types[j], probe = probes[k], $
                hopel3 = hopel3
end

