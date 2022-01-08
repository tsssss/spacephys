;+
; Type: procedure.
; Purpose: FAST, sdt interface, save ESA data to *.tplot file.
;   Default interface {out:'ion_en_spec','ion_pa_spec','ion_nflux',
;   'ion_j','ion_eflux','ion_n','ion_para_spec','ion_perp_spec','ion_anti_spec',
;   electron counterparts,'dis','ilat','mlt','alt','fa_pos','fa_vel'}
; Parameters: none.
; Keywords:
;   fn, in, string, optional. Output filename.
;   t1, in, double, required. Start time for data.
;   t2, in, double, required. End time for data.
;   ionen, in, dblarr[2], optional. Ion energy limit for integration, 
;       used to eliminate contaminated low energy ions.
;   noshift, in, boolean, optional. Set PA in [0,360], default [-90,270].
;   omni, in, boolean, optional. Set to get omni-dir energy spectrogram only.
; Notes: The tplot in /local/usr/sdt/ will NOT work. There can be 2 or more 
;   versions of tdas. We need the FAST part in the tdas at /local/usr/sdt/ 
;   and the tplot part in the regular tdas6 or tdas7. Thus, I copied tdas7 
;   in slib at /home/shengt/works/code/slib/, overwriting tdas6 at 
;   /local/usr/TDAS/, then run fast_sdt_save_esa_data, several runs are
;   needed to discover conflicting FAST func/pro in tdas7, remove them all.
;       Set t1 and t2 as near as possible, otherwise long time range is very
;   slow to run.
;-

pro fast_sdt_save_esa_data, filename = fn, t1 = t1, t2 = t2, $
    ionen = ionen, noshift = noshift, omni = omni   ; TODO omni.
    
    if n_elements(t1) eq 0 or n_elements(t2) eq 0 then $
        message, 'no start or end time ...'
  
    shift90 = n_elements(noshift) ? 0: 1
    shift90 = ~keyword_set(noshift)

    mu = '!9'+string(109b)+'!X'
    gamma = '!9'+string(71b)+'!X'
    para = '||'
    perp = '!9'+string(94b)+'!X'
  
    ; init names.
    ionenspec = 'ion_en_spec'
    ionpaspec = 'ion_pa_spec'
    eleenspec = 'ele_en_spec'
    elepaspec = 'ele_pa_spec'
    ionnflux  = 'ion_nflux'
    ioneflux  = 'ion_eflux'
    ionj      = 'ion_j'
    ionn      = 'ion_n'
    ionp      = 'ion_p'
    elenflux  = 'ele_nflux'
    eleeflux  = 'ele_eflux'
    elej      = 'ele_j'
    elen      = 'ele_n'
    elep      = 'ele_p'

    ionparaspec = 'ion_para_spec'   ; [0,30].
    ionperpspec = 'ion_perp_spec'   ; [40-140].
    ionantispec = 'ion_anti_spec'   ; [150-180].
    eleparaspec = 'ele_para_spec'   ; [0,30].
    eleperpspec = 'ele_perp_spec'   ; [60-120].
    eleantispec = 'ele_anti_spec'   ; [150-180].
    ionparaang = [0,30]           ; from fast berkeley website.
    ionperpang = [40,140]
    ionantiang = [150,180]
    eleparaang = [0,30]
    eleperpang = [60,120]
    eleantiang = [150,180]

    ; init scales.
    ionspecscale = [1e4, 1e8]
    elespecscale = [1e6, 1e10]
  
    ; ==============
    ; particle data.
    ; ==============
    ; ion energy spectrum.
    get_en_spec, 'fa_ies', unit = 'eflux', name = ionenspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2
    options, ionenspec, 'ytitle', 'Energy!C(eV)!Cion'
    get_data, ionenspec, data = tmp
    times = tmp.x
    get_en_spec, 'fa_ies', unit = 'eflux', name = ionparaspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = ionparaang
    options, ionparaspec, 'ytitle', 'Energy!C(eV)!Cion para'
    get_en_spec, 'fa_ies', unit = 'eflux', name = ionperpspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = ionperpang
    options, ionperpspec, 'ytitle', 'Energy!C(eV)!Cion perp'
    get_en_spec, 'fa_ies', unit = 'eflux', name = ionantispec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = ionantiang
    options, ionantispec, 'ytitle', 'Energy!C(eV)!Cion anti'
  
    ; ion pitch angle spectrum.
    get_pa_spec, 'fa_ies', unit = 'eflux', name = ionpaspec, $
        energy = ionen, shift90 = shift90, $
        retrace = 1, /calib, t1 = t1, t2 = t2
    options, ionpaspec, 'ytitle', 'Pitch Angle!C(deg)!Cion'
  
    ; ele energy spectrum.
    get_en_spec, 'fa_ees', unit = 'eflux', name = eleenspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2
    options, eleenspec, 'ytitle', 'Energy!C(eV)!Cele'
    get_en_spec, 'fa_ees', unit = 'eflux', name = eleparaspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = eleparaang
    options, eleparaspec, 'ytitle', 'Energy!C(eV)!Cele para'
    get_en_spec, 'fa_ees', unit = 'eflux', name = eleperpspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = eleperpang
    options, eleperpspec, 'ytitle', 'Energy!C(eV)!Cele perp'
    get_en_spec, 'fa_ees', unit = 'eflux', name = eleantispec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, angle = eleantiang
    options, eleantispec, 'ytitle', 'Energy!C(eV)!Cele anti'
  
    ; ele pitch angle spectrum.
    get_pa_spec, 'fa_ees', unit = 'eflux', name = elepaspec, $
        retrace = 1, /calib, t1 = t1, t2 = t2, shift90 = shift90
    options, elepaspec, 'ytitle', 'Pitch Angle!C(deg)!Cele'
  
    ; for all energy spectrum.
    vars = [ionenspec, ionparaspec, ionperpspec, ionantispec, $
         eleenspec, eleparaspec, eleperpspec, eleantispec]
    ylim, vars, 4, 40000, 1
    ; for all pitch angle spectrum.
    vars = [ionpaspec, elepaspec]
    if keyword_set(noshift) then ylim, vars, 0, 360 $
    else ylim, vars, -90, 270
  
    ; for all ion spectrum.
    vars = [ionenspec, ionpaspec, ionparaspec, ionperpspec, ionantispec]
    zlim, vars, ionspecscale[0], ionspecscale[1], 1
  
    ; for all ele spectrum.
    vars = [eleenspec, elepaspec, eleparaspec, eleperpspec, eleantispec]
    zlim, vars, elespecscale[0], elespecscale[1], 1
  
    ; for all spectrum.
    vars = [ionenspec, ionparaspec, ionperpspec, ionantispec, ionpaspec, $
         eleenspec, eleparaspec, eleperpspec, eleantispec, elepaspec]
    options, vars, 'ztitle', 'Log Eflux!C(eV/cm!U2!N-s-sr-eV)'
    options, vars, 'spec', 1
    options, vars, 'no_interp', 1
  
    coef = 1.6e-9   ; convert #/cm^2-s to uA/m^2.
    ; ion number flux.
    get_2dt, 'j_2d', 'fa_ies', name = ionnflux, energy = ionen, $
        t1 = t1, t2 = t2
    get_data, ionnflux, data = tmp, limits = {ytitle:'(#/cm!U2!N-s)', $
        labels:'ion nflux'}
    tmp.y = tmp.y * coef
    store_data, ionj, data = tmp, limits = {ytitle:'(!9'+mu+'A/m!U2!N)', $
        labels:'nqUi!D||!N'}
  
    ; ion energy flux.
    get_2dt, 'je_2d', 'fa_ies', name = ioneflux, energy = ionen, $
        t1 = t1, t2 = t2
    options, ioneflux, 'ytitle', '(ergs/cm!U2!N-s)'
    options, ioneflux, 'labels', 'KEi'
    
    ; ion density.
    get_2dt, 'n_2d', 'fa_ies', name = ionn, energy = ionen, t1 = t1, t2 = t2
    options, ionn, 'ytitle', '(1/cm!U3!N)'
    options, ionn, 'labels', 'Ni'

    ; ion pressure.
    get_2dt, 'p_2d', 'fa_ies', name = ionp, energy = ionen, t1 = t1, t2 = t2
    get_data, ionp, t0, tmp
    store_data, ionp, t0, tmp[*,[1,2]], limits = {ytitle:'(nPa)', $
        labels:'Pi!D'+[perp,para]+'!N'}
  
    ; ele number flux.
    get_2dt, 'j_2d', 'fa_ees', name = elenflux, $
        t1 = t1, t2 = t2
    options, elenflux, limits = {ytitle:'(#/cm!U2!N-s)', $
        labels:'ele nflux'}
    get_data, elenflux, data = tmp
    tmp.y = tmp.y * coef
    store_data, elej, data = tmp, limits = {ytitle:'(!9'+mu+'A/m!U2!N)', $
        labels:'nqUe!D||!N'}
  
    ; ele energy flux.
    get_2dt, 'je_2d', 'fa_ees', name = eleeflux, $ 
        t1 = t1, t2 = t2
    options, eleeflux, 'ytitle', '(ergs/cm!U2!N-s)'
    options, eleeflux, 'labels', 'KEe'
    
    ; ele density.
    get_2dt, 'n_2d', 'fa_ees', name = elen, t1 = t1, t2 = t2
    options, elen, 'ytitle', '(1/cm!U3!N)'
    options, elen, 'labels', 'Ne'

    ; ele pressure.
    get_2dt, 'p_2d', 'fa_ees', name = elep, t1 = t1, t2 = t2
    get_data, elep, t0, tmp
    store_data, elep, t0, tmp[*,[1,2]], limits = {ytitle:'(nPa)', $
        labels:'Pe!D'+[perp,para]+'!N'}
  
    ; get orbit data.
    get_fa_orbit, t1, t2
    get_data, 'ILAT', data = tmp
    store_data, 'fa_ilat', data = tmp, limits = {ytitle:'ILat'}
    get_data, 'MLT', data = tmp
    store_data, 'fa_mlt', data = tmp, limits = {ytitle:'MLT'}
    get_data, 'ALT', data = tmp
    store_data, 'fa_alt', data = tmp, limits = {ytitle:'ALT'}
    store_data, 'fa_dis', data = {x:tmp.x, y:tmp.y/6356.7523+1}, $
        limits = {ytitle:'Dist (Re)'}
    get_data, 'ORBIT', data = tmp & orbstr = string(tmp.y[0],format='(I05)')
    vars = ['ORBIT','ALT','ILAT','ILNG','MLT']
    for i = 0, n_elements(vars)-1 do store_data, vars[i], /delete

    ; prepare filename.
    if n_elements(fn) eq 0 then begin
        t0str = time_string(tmp.x[0], format = 2, precision = -3)
        t0str = strmid(t0str,0,4)+'_'+strmid(t0str,4,4)
        fn = '~/fa_sdt_esa_'+t0str+'_'+orbstr
    endif

    vars = [ionenspec,ionparaspec,ionperpspec,ionantispec,ionpaspec, $
        eleenspec,eleparaspec,eleperpspec,eleantispec,elepaspec, $
        ionn,ionj,ionp,ionnflux,ioneflux, $
        elen,elej,elep,elenflux,eleeflux, $
        'fa_pos','fa_vel','fa_ilat','fa_mlt','fa_alt','fa_dis']

    tplot_save, vars, filename = fn
end
