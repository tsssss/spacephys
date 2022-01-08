;+
; Type: crib.
; Purpose: Plot HOPE L3 time-energy spectrogram.
; Parameters:
;   t0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   type0, in, string, opt. 'proton','oxygen','helium','electron'.
;       Default is 'proton'.
; Keywords:
;   probe, in, string, opt. 'a','b'. Default is 'a'
;   hopel3, in/out, struct, opt. Pass among plotting routines to fast plot.
;   pitch, in, float or fltarr[2], opt. Pich angle in deg. Plot single angle
;       or certain angle range. Default is to sum all angle bins.
; Notes: none.
; Dependence: slib.
; History:
;   2016-04-08, Sheng Tian, create.
;-

pro plot_hope_l3_energyspec, t0, type0, probe = probe, $
    hopel3 = hopel3, pitch = pa0, position = tpos

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
    enl3s = enl3s
;    enl3s = enl3s[0:*:2]
    nenl3 = n_elements(enl3s[0,*])
;    datl3 = datl3[idx,*]
;    datl3 = datl3[0:*:2,*]
    
    

; **** the data for spectrogram.
    case n_elements(pa0) of
        0: paidx = findgen(npal3)
        1: tmp = min(pal3s-pa0, paidx, /absolute)
        2: paidx = where(pal3s ge min(pa0) and pal3s le max(pa0))
    endcase
    if n_elements(paidx) eq 1 then tdat = reform(datl3[*,*,paidx]) else $
        tdat = total(datl3[*,*,paidx],3, /nan)
    tang = pal3s
    tdis = enl3s

    case n_elements(pa0) of
        0: ytitl = 'All Pictch'
        1: ytitl = 'Pitch '+sgnum2str(pa0)+ ' deg'
        2: ytitl = 'Pitch '+sgnum2str(min(pa0))+'-'+sgnum2str(max(pa0))+' deg'
    endcase
    ytitl = 'Energy (eV)!C'+ytitl
    
    case type of
        'proton': ztitl = 'H+'
        'oxygen': ztitl = 'O+'
        'heilum': ztitl = 'He+'
        'electron': ztitl = 'e+'
    endcase
    ztitl = ztitl+' flux (s!U-1!Ncm!U-2!Nster!U-1!NkeV)'

    ; remove nan.
;    idx = where(finite(tdat,/nan))
;    tdat[idx] = 0
    

    idx = where(tdat ne 0)
    min0 = min(tdat[idx],/nan)
    max0 = max(tdat[idx],/nan)
    nztick = 256
    if n_elements(zrng) eq 0 then $
        zrng = [floor(alog10(min0)),ceil(alog10(max0))-2]
    titl = 'RBSP-'+strupcase(probe)+' HOPE L3 '+type+' Eflux!C'+$
        time_string(uts[0])+' - '+time_string(uts[1],tformat='hh:mm:ss')

    tvar = 'rb'+probe+'_'+type+'_enspec'
    store_data, tvar, uts, tdat, tdis
    options, tvar, 'ytitle', ytitl
    options, tvar, 'yrange', [100,50000]
    options, tvar, 'ylog', 1
    options, tvar, 'ystyle', 1
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    options, tvar, 'zrange', [3e4,3e5]
    options, tvar, 'zlog', 1
    options, tvar, 'ztitle', ztitl
    options, tvar, 'zticks', 1
    options, tvar, 'xticklen', -0.02
    options, tvar, 'yticklen', -0.01
    options, tvar, 'xminor', 6
    
end

utr = time_double(['2013-05-01/07:00','2013-05-01/08:30'])
log = 0
types = ['proton','oxygen']
probes = ['b']

rootdir = shomedir()+'/rbsp_hope'

poss = sgcalcpos(n_elements(types))
poss[2,*] = 0.85
poss[0,*] = 0.15

for k = 0, n_elements(probes)-1 do begin
    ofn = rootdir+'/rbsp'+probes[k]+'/hope_l3_enspec_'+ $
        time_string(utr[0],tformat='YYYY_MMDD_hh')+'.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 10, ysize = 5, /inch
    sgindexcolor
    loadct2, 43
    !p.color = 0
    !p.background = 0

    for j = 0, n_elements(types)-1 do $
        plot_hope_l3_energyspec, utr, types[j], probe = probes[k], $
            hopel3 = hopel3, pitch = 162, position = poss[*,j]
    
    
    vars = 'rb'+probes[k]+'_'+types+'_enspec'
    tplot, vars, trange = utr, position = poss, /noerase, /novtitle

    sgclose

endfor

    
end


