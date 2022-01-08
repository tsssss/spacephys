;+
; Type: crib.
; Purpose: Plot HOPE L2 2d distribution at certain time.
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

pro plot_hope_l2_pitch2d, t0, type0, probe = probe, $
    unit = unit0, log = log, hopel2 = hopel2, zrange = zrng0

    ; date and time.
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
    if n_elements(zrng0) ne 0 then zrng = zrng0

    print, 'rbsp'+probe+', '+type+', '+time_string(ut0)

    ; suppress math exception.
    !except = 0

    ; constants.
    npxl = 5    ; # of polar pixels.
    nsec = 16   ; # of azimuthal sectors.
    nenb = 72   ; # of energy bins.

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
    load = 0
    if n_elements(hopel2) eq 0 then load = 1 else begin
        get_data, 'hopel2', t0, info
        if size(info,/type) ne 8 then load = 1 else if info.probe ne probe then load = 1
    endelse
    
    if load eq 1 then begin
        timespan, date, 1, /day
        hopel2 = sread_rbsp_hope_l2(utr, probes = probe)
        store_data, 'hopel2', t0, {probe:probe}
        
        ; load b_uvw.
        emfisis = sread_rbsp_emfisis(ut0, probe = probe)
        buts = sfmepoch(emfisis.epoch,'unix',/tt2000)
        bgse = scalcbg(emfisis.mag)
        store_data, 'b_gse', buts, bgse
        
        ; load spice kernal.
        rbsp_load_spice_kernels, trange = minmax(t0), probes = probe
        rbsp_uvw2gse, 'b_gse', newname = 'b_uvw', probe = probe, /inverse, /no_spice_load
    endif

    
    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel2) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel2) eq enidx)
    tmp = min(hopel2.(epidx)-stoepoch(ut0,'unix'),rec,/absolute)
    i0 = where(tag_names(hopel2) eq type0)
    dtidx = (type eq 'electron')? 'EPOCH_ELE_DELTA': 'EPOCH_ION_DELTA'
    dtidx = where(tag_names(hopel2) eq dtidx)
    t_midl = sfmepoch(hopel2.(epidx)[rec],'unix')   ; middle time of a spin.
    t_half = hopel2.(dtidx)[rec]*1e-3               ; half of spin duration.
    t_step = 10.417e-3                              ; 10.417 ms, esa sweep time.
    

    ; level 2.
    datl2 = double(reform((hopel2.(i0))[rec,*,*,*]))    ; in [72,16,5].
    idx = where(datl2 eq -1e31, cnt)
    if cnt ne 0 then datl2[idx] = !values.d_nan
    idx = where(datl2 eq 0, cnt)
    if cnt ne 0 then datl2[idx] = !values.d_nan
    
    enl2s = reform(hopel2.(enidx)[rec,*])               ; in [72], energy bins.
    nenl2 = n_elements(enl2s)
    
    pal3s = [4.5,18,36,54,72,90,108,126,144,162,175.5]  ; L3 pitch angles.
    npal3 = n_elements(pal3s)
    dpas = [4.5,9,9,9,9,9,9,9,9,9,4.5]

        
    ; get uniform PA from L2 data.
    get_data, 'b_uvw', buts, buvw
    buvw = sunitvec(reform(sinterpol(buvw,buts,t_midl)))
    
    scmode = reform(hopel2.sector_collapse_cntr[rec,*]) ; sector collapse mode.
    
    angas = 90-(findgen(npxl)*36+18)        ; polar angle, in deg.
    pal2s = dblarr(nenb,nsec,npxl)
    
    for i = 0, npxl-1 do begin
        tnsec = nsec/scmode[i]
        angbs = (findgen(nsec)-nsec/2)*360d/nsec; azimuthal angle, in deg.
        for j = 0, tnsec-1 do begin
            idx1 = j*scmode[i]              ; treat the collpaed sectors.
            idx2 = idx1+scmode[i]-1
            for k = 0, nenb-1 do begin
                tmp = 360d/nsec/nenb*k      ; adjust the angle for ebin.
                if (j mod 2) eq 1 then tmp = 360d/nsec-tmp
                tangb = (mean(angbs[idx1:idx2])+tmp)*rad
                tanga = angas[i]*rad
                puvw = [cos(tanga)*[cos(tangb),sin(tangb)],sin(tanga)]
                pal2s[k,idx1:idx2,i] = sang(puvw, buvw, /deg)
            endfor
        endfor
    endfor
    
;    for i = 0, nsec-1 do begin
;        for j = 0, nenb-1 do begin
;            tmp = dangb/nenb*j  ; the angle at this energy bin.
;            if (i mod 2) eq 1 then tmp = dangb-tmp  ; see Fig. 22.
;            tangb = (angbs[i]+tmp)*rad   ; adjust to the angle of this energy bin.
;            for k = 0, npxl-1 do begin
;                tanga = angas[k]*rad
;                puvw = [cos(tanga)*[cos(tangb),sin(tangb)],sin(tanga)]
;                pal2s[j,i,k] = sang(puvw, buvw, /deg)
;            endfor
;        endfor
;    endfor

    ; bin the L2 pitch angles and take the mean value.
    datl3 = dblarr(nenl2,npal3)
    datl2 = reform(datl2,[nenl2,nsec*npxl])
    pal2s = reform(pal2s,[nenl2,nsec*npxl])
    for i = 0, nenb-1 do begin
        for j = 0, npal3-1 do begin
            tmp = pal2s[i,*]
            idx = uniq(tmp,sort(tmp))
            tmp = datl2[i,idx]
            idx = where(abs(pal2s[i,idx]-pal3s[j]) le dpas[j], cnt)
            if cnt eq 0 then datl3[i,j] = !values.d_nan $
            else datl3[i,j] = mean(tmp[idx],/nan)
        endfor
    endfor
    enl3s = enl2s
    nenl3 = nenl2


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
    titl = 'RBSP-'+strupcase(probe)+' HOPE L2 Pitch Angle!C'+type+' Eflux!C'+$
        time_string(t_midl-t_half)+' - '+time_string(t_midl+t_half,tformat='hh:mm:ss')
        
    rootdir = shomedir()+'/rbsp_hope'
    ofn = rootdir+'/rbsp'+probe+'/'+type+'/hope_l2_pitch2d_'+type+'_'+ $
        time_string(t_midl,tformat='YYYY_MMDD_hhmm_ss')+'.pdf'
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    sgindexcolor, 43, file = 'ct2'
    sgdistr2d, tdat, tang, tdis, position = tpos, zrange = zrng, $
        title = titl, xtitle = xtitl, ncolor = 10
    sgclose

end


; Mark's events.
utr = time_double(['2015-12-14/12:18','2015-12-14/12:38'])
uts = smkarthm(utr[0],utr[1],12,'dx')
log = 1
unit = 'energy'
types = ['proton']
probes = ['a']

;; Wygant's events.
;utr = time_double(['2013-05-01/07:38','2013-05-01/07:45'])
;uts = smkarthm(utr[0],utr[1],12,'dx')
;log = 1
;unit = 'energy'
;types = ['proton']
;probes = ['b']

;for k = 0, n_elements(probes)-1 do $
;    for j = 0, n_elements(types)-1 do $
;        for i = 0, n_elements(uts)-1 do $
;            plot_hope_l3_pitch2d, uts[i], types[j], unit = unit, log = log, $
;                hopel3 = hopel3, probe = probes[k]

for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l2_pitch2d, uts[i], types[j], unit = unit, log = log, $
                hopel2 = hopel2, probe = probes[k]
end
