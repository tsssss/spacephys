
pro plot_hope_l3_keflux, utr, type, probe = probe, hopel3 = hopel3, $
    max_energy = maxen0, min_energy = minen0, pre0=pre00

    ; date and time.
    date = time_string(utr[0],tformat='YYYY-MM-DD')
    utr0 = time_double(date)+[0,86400d]

    if n_elements(probe) eq 0 then probe = ['a']
    if n_elements(type) eq 0 then type = 'proton' else type = strlowcase(type)


    ; min energy will be used in all moments. Energy bins below may be affected by Vsc.
    ; max energy will be used in KE flux.
    if n_elements(maxen0) ne 0 then maxens = maxen0
    if n_elements(minen0) ne 0 then minens = minen0


; **** constants.
    npxl = 5d   ; # of pixels.
    
    rad = !dpi/180d
    deg = 180d/!dpi
    re = 6378d & re1 = 1d/re
    pre0 = 'rbsp'+probe+'_'
    if n_elements(pre00) ne 0 then pre0=pre00

    ; var name for each species.
    case type of
        'electron': type0 = 'FEDU'
        'proton': type0 = 'FPDU'
        'oxygen': type0 = 'FODU'
        'helium': type0 = 'FHEDU'
    endcase

    ; mass but contains unit conversion coef.
    case type of
        'electron': mass0 = 1d/1836
        'proton': mass0 = 1d
        'oxygen': mass0 = 16d
        'helium': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)    ; E in eV, mass in kg.


; **** load hope l3 data.
    load = 0
    if n_elements(hopel3) eq 0 then load = 1 else begin
        if tnames(pre0+'hopel3') ne 0 then begin
            get_data, pre0+'hopel3', tmp, info
            if info.probe ne probe then load = 1
            if max(info.utr-utr) ne 0 then load = 1
        endif
    endelse

    if load eq 1 then begin
        timespan, date, 1, /day
        hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
        store_data, pre0+'hopel3', utr0[0], {probe:probe,utr:utr}
    endif
    

; **** read data. 
    epvname = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    envname = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    
    epidx = where(tag_names(hopel3) eq epvname)
    enidx = where(tag_names(hopel3) eq envname)
    uts = sfmepoch(hopel3.(epidx),'unix')   ; time in the middle of a spin.
    nrec = n_elements(uts)

    i0 = where(tag_names(hopel3) eq type0)  ; idx for spectrogram in the structure.
    datl3 = reform((hopel3.(i0)))           ; in [nrec,nen,npa].
    enl3s = reform(hopel3.(enidx))          ; in [nrec,nen], in eV.
    pal3s = reform(hopel3.pitch_angle)      ; in [npa], in deg.
    pal3s = pal3s*rad                       ; convert to rad.
    npal3 = n_elements(pal3s)
    dpas = [9d,dblarr(npal3-2)+18,9d]       ; dangle around each pitch angle.
    dpas = dpas*rad
    idx = where(datl3 eq -1e31, cnt)        ; remove fill value.
    if cnt ne 0 then datl3[idx] = !values.d_nan    
    idx = where(finite(datl3,/nan), cnt)
    if cnt ne 0 then datl3[idx] = 0d
    
    if n_elements(maxens) eq 1 then maxens = dblarr(nrec)+maxens[0]
    if n_elements(minens) eq 1 then minens = dblarr(nrec)+minens[0]
    
    

; **** density.
    ns = dblarr(nrec)                   ; in #/cm^3.
    
    for i = 0, nrec-1 do begin
        tes = reform(enl3s[i,*])        ; the energy bins, in eV.
        nen = n_elements(tes)           ; # of energy bins.
        des = [tes[0],tes[1:nen-1]-tes[0:nen-2]]; dE at each energy bin, in eV.
        nflux = reform(datl3[i,*,*])    ; in #/s-cm^2-sr-keV.
        vs = sqrt(2*tes/mass0)*1e3      ; energy converted to velocity, in km/s.
        for j = 0, nen-1 do begin       ; integrate over the energy bins.
            if n_elements(minens) ne 0 then $
                if tes[j] lt minens[i] then continue
            for k = 0, npal3-1 do $     ; integrate each pitch angle.
                ns[i]+= nflux[j,k]/vs[j]*des[j]*1d-3* $
                    2*!dpi*sin(pal3s[k])*dpas[k]
        endfor
    endfor
    store_data, pre0+'n_'+type, uts, ns, $
        limits = {ytitle:'RBSP-'+strupcase(probe)+'!CDensity!C(cm!U-3!N)', $
        labels:type+'!C  density', ylog:1, yrange:[0.01,10]}
    
    
; **** kinetic energy flux.
    kes = dblarr(nrec)                  ; in mW/m^2.
    minpa = 40d*rad
    maxpa = 140d*rad
    
    for i = 0, nrec-1 do begin
        tes = reform(enl3s[i,*])        ; the energy bins, in eV.
        nen = n_elements(tes)           ; # of energy bins.
        des = [tes[0],tes[1:nen-1]-tes[0:nen-2]]; dE at each energy bin, in eV.
        nflux = reform(datl3[i,*,*])    ; in #/s-cm^2-sr-keV.
        for j = 0, nen-1 do begin       ; integrate over the energy bins.
            if n_elements(maxens) ne 0 then $
                if tes[j] gt maxens[i] then continue
            if n_elements(minens) ne 0 then $
                if tes[j] lt minens[i] then continue
            for k = 0, npal3-1 do begin ; integrate each pich angle.
                if pal3s[k] gt minpa and pal3s[k] le maxpa then continue    ; skip perp flux.
                kes[i]+= nflux[j,k]*des[j]*tes[j]*1d-3* $
                    2*!dpi*sin(pal3s[k])*dpas[k]*cos(pal3s[k])
            endfor
        endfor
    endfor
    kes*= 1.6e-12   ; from eV/s-cm^2 to mW/m^2.
    
    store_data, pre0+'keflux_'+type, uts, kes, $
        limits = {ytitle:'RBSP-'+strupcase(probe)+'!CKEflux!C(mW/m!U2!N)', $
            labels:type+'!C  KEflux!C  in situ'}


; **** temperature. f(v,x,t) = c*exp(-mv^2/2kT). f*v*v^2*dv = J*dW,
; so f = J*dW/v^3/dv, dW = m*v*dv, f = J*m/v^2 = 2*J*m^2/W
; f(v,x,t) = d*J/W = c*exp(-W/T), so log(J/W) = C-W/T.
    tperps = fltarr(nrec)
    tparas = fltarr(nrec)
    idxperp = (npal3-1)/2+[-1,0,1]
    idxpara = [0,1,npal3-2,npal3-1]
    case type of
        'proton': tminen = 2e3
        'oxygen': tminen = 2e3
        'electron': tminen = 5e2
        else: tminen = 1e2
    endcase
    case type of
        'proton': tmaxen = 10e4
        'oxygen': tmaxen = 1e4
        'electron': tmaxen = 1e4
        else: tmaxen = 1e2
    endcase
    minflux = 1.5e2
    for i = 0, nrec-1 do begin
        ; perp slice.
        tes = reform(enl3s[i,*])        ; the energy bins, in eV.
        nflux = reform(datl3[i,*,*])    ; in #/s-cm^2-sr-keV.
        tys = total(nflux[*,idxperp],2,/nan)/n_elements(idxperp)
        idx = where(finite(tys) and tys ge minflux and tes ge tminen, cnt)
        if cnt le 2 then continue
        tes = tes[idx]
        tys = alog(tys[idx]/tes)
        res = linfit(tes,tys)
        tperps[i] = -1d/res[1]
;        print, -1d/res[1]
;        plot, tes, tys, psym = 1
;        oplot, tes, res[0]+res[1]*tes, linestyle = 3, color = 6
;        stop

        ; para slice.
        tes = reform(enl3s[i,*])        ; the energy bins, in eV.
        nflux = reform(datl3[i,*,*])    ; in #/s-cm^2-sr-keV.
        tys = total(nflux[*,idxpara],2,/nan)/n_elements(idxpara)
        idx = where(finite(tys) and tys ge minflux and tes ge tminen, cnt)
        if cnt le 2 then continue
        tes = tes[idx]
        tys = alog(tys[idx]/tes)
        res = linfit(tes,tys)
        tparas[i] = -1d/res[1]
;        print, -1d/res[1]
;        plot, tes, tys, psym = 1
;        oplot, tes, res[0]+res[1]*tes, linestyle = 3, color = 6
;        stop
    endfor

    store_data, 'tperp', uts, tperps*0.5, limits = {yrange:[1e2,1e5], ylog:1}
    store_data, 'tpara', uts, tparas*0.5, limits = {yrange:[1e2,1e5], ylog:1}
    tplot, ['tperp','tpara'], trange = utr

end

utr2 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr2 = time_double(['2013-05-01/07:20','2013-05-01/07:50'])

plot_hope_l3_keflux, utr2, 'oxygen', probe = 'a', min_energy = 30
plot_hope_l3_keflux, utr2, 'proton', probe = 'a', min_energy = 30
plot_hope_l3_keflux, utr2, 'electron', probe = 'a', min_energy = 200


stop



device, decomposed = 0
loadct2, 43

tplot_options, 'labflag', -1
tplot_options, 'constant', 0
tplot_options, 'num_lab_min', 6

types = ['electron','proton','oxygen']
probes = ['a','b']

utr = time_double(['2012-11-14/04:00','2012-11-14/05:00'])

for i = 0, n_elements(probes)-1 do begin
    
    prb = probes[i]
    pre0 = 'rbsp'+prb+'_'
    ; en spec.    
    plot_hope_l3_enspec, utr, probe = prb, 'electron', hopel3 = hopel3
    plot_hope_l3_enspec, utr, probe = prb, 'oxygen', hopel3 = hopel3
    plot_hope_l3_enspec, utr, probe = prb, 'proton', hopel3 = hopel3
    
;    vars = [pre0+types]+'_enspec'
;    ofn = shomedir()+'/rbsp'+prb+'_enspec.pdf'
;    sgopen, ofn, xsize = 8, ysize = 11.5, /inch
;    tplot, vars, trange = utr
;    sgclose
    
    ; electron.
    plot_hope_l3_keflux, utr, probe = prb, 'electron', hopel3 = hopel3
    
    ; ions.
    if tnames(pre0+'maxen') eq '' then begin
        tplot, pre0+'oxygen_enspec', trange = utr
        stplot_add_line, vxs, vys
        store_data, pre0+'maxen', vxs, vys
    endif
    
    plot_hope_l3_keflux, utr, probe = prb, 'oxygen', hopel3 = hopel3, max_energy = vys
    plot_hope_l3_keflux, utr, probe = prb, 'proton', hopel3 = hopel3, max_energy = vys
    
    get_data, pre0+'oxygen_enspec', limits = lim
    store_data, pre0+'oxygen_enspec_comb', data = pre0+['oxygen_enspec','maxen'], limits = lim
    get_data, pre0+'proton_enspec', limits = lim
    store_data, pre0+'proton_enspec_comb', data = pre0+['proton_enspec','maxen'], limits = lim

    ; map.
    ; load pos.
    efwl3 = sread_rbsp_efw_l3(utr, probes = probe)
    if size(efwl3,/type) ne 8 then return
    uts = sfmepoch(efwl3.epoch,'unix',/epoch16)
    
    store_data, pre0+'pos_gse', uts, efwl3.pos_gse/6378d
        
    ; load b.
    emfisis = sread_rbsp_emfisis(utr, probes = probe)
    if size(emfisis,/type) ne 8 then return
    uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B magnitude!C(nT)'}
    
    scalc_map_coef, pre0+'pos_gse', pre0+'b', model = 't89', coord = 'gse', $
        prefix = pre0
    
    get_data, pre0+'map_coef', t0, coef
    vars = pre0+['keflux_'+types]
    for j = 0, n_elements(vars)-1 do begin
        get_data, vars[j], t0, dat, limits = lim
        store_data, vars[j]+'_map', t0, dat*coef, limits = lim
        options, vars[j]+'_map', 'labels', types[j]+'!C  KE flux'
    endfor
    
    posl = [0.15,0.1,0.4,0.9]
    posr = [0.65,0.1,0.9,0.9]
    labs = pre0+'fpt_mlat'
    
    ofn = shomedir()+'/'+pre0+'keflux.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 11.5, ysize = 8, /inch
    
    vars = [pre0+['keflux_'+types]]+'_map'
    tpos = sgcalcpos(n_elements(vars), position = posl)
    tplot, vars, trange = utr, var_label = labs, position = tpos, /noerase, $
        title = 'RBSP-'+strupcase(prb)+' KE Flux @100 km'
    
    vars = pre0+[types+'_enspec']
    tpos = sgcalcpos(n_elements(vars), position = posr)
    tplot, vars, trange = utr, var_label = labs, position = tpos, /noerase, $
        title = 'HOPE Energy Spec e-, H+, O+'
    
    vars = pre0+[['oxygen','proton']+'_enspec']
    for j = 0, n_elements(vars)-1 do begin
        get_data, vars[j], limits = lim
        plot, utr, lim.yrange, /nodata, /noerase, position = tpos[*,j+1], $
            xstyle = 5, ystyle = 5, /ylog
        tvar = pre0+'maxen'
        get_data, tvar, vxs, vys
        oplot, vxs, vys;, psym = -1
    endfor
    
    sgclose
    
endfor

end
