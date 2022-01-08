;+
; Type: procedure.
; Purpose: FAST, sdt interface, save E, B field to *.tplot file.
;   Default interface {out:'dB_fac_v','B_model','E_ALONG_V','E_NEAR_B',
;   'ilat','mlt','alt','dis','fa_pos','fa_vel'}
; Parameters: none.
; Keywords:
;   fn, in, string, optional. Output filename.
;   tc, in, double, optional. Suggestive time center.
;   t1, in, double, optional. Mandatory start time for E.
;   t2, in, double, optional. Mandatory end time for E.
;   nodelete = nodelete, in, boolean, optionsl. Set to store unsued vars.
;   vars = vars, in, strarr[n], optional. Varnames to save.
; Notes: The tplot in /local/usr/sdt/ will NOT work. cd to tdas6 at
;   /local/usr/TDAS/ or tdas7 in slib at /home/shengt/works/code/slib/.
;       Need SDT data loaded. For E field: 'FastFieldsMode_1032', 'V4_S', 
;   'V8_S', 'V1-V2_S', 'V5-V8_S', 'V1-V4_S'. 'FastFieldsMode_1032' contains
;   flag for E field, 255 means no data, otherwise have data. Use 'V1-V2_S'
;   by default, when it is bad, use 'V5-V8' instead. For B field: 'MagX', 
;   'MagY', 'MagZ', '1032-spinPhase'. See ucla_mag_despin, fa_fields_despin.
;       B spans all orbit, E spans to when flag = 255. Set tc to suggest 
;   center time for E to span. Because E field has very high time resolution,
;   data file can be large, set t1/t2 for E time range is useful.
;       There are FAST orbits have several sectors of E field, if so, set tc
;   to indicate which one you want.
;       Default to load E12, if no E12, load E158, see fa_fields_despin.
; Dependence: tdas, sdt.
; Author: Sheng Tian.
; History: 2013-11-15, Sheng Tian, create.
;-
;
pro fast_sdt_save_field_data, filename = fn, tc = tc, t1 = t1, t2 = t2, $
    vars = var0, nodelete = nodelete
    
    ; load flag from 'DataHdr_1032'.
    ; interface: {out: etr}
    tmp = get_fa_fields('DataHdr_1032')
    fmode = reform(tmp.comp1[13,*]) & t0 = tmp.time
    idx = where(fmode ne 255, cnt)
    if cnt eq 0 then message, 'no valid E ...'
    del = idx[1:*]-idx[0:cnt-2]
    secs = where(del gt 1, nsec)    ; sectors.
    if nsec eq 0 then etr = [t0[idx[0]],t0[idx[cnt-1]]] $    ; 1 sector.
    else begin      ; multiple sectors.
        sec = secs
        if secs[0] gt idx[0] then sec = [idx[0],sec]
        if secs[nsec-1] lt idx[cnt-1] then sec = [sec,idx[cnt-1]]
        nsec = n_elements(sec)-1
        secs = lonarr(nsec,2)       ; [idx0,idx1].
        secs[*,0] = sec[0:nsec-1]+1 & secs[0] = sec[0]
        secs[*,1] = sec[1:nsec]
        secs = idx[secs]
        stop
        if n_elements(tc) ne 0 then begin   ; tc is set.
            if tc lt nsec then etr = [t0[secs[tc,0]],t0[secs[tc,1]]] $  ; tc is #.
            else begin                      ; tc is section time.
                for i = 0, nsec-1 do $
                    if t0[secs[i,0]] le tc and t0[secs[i,1]] ge tc then break
                etr = [t0[secs[i,0]],t0[secs[i,1]]]
            endelse
        endif else begin  ; user define sec.
            tc = 1000
            while tc gt nsec do begin
                print, 'which sector out of '+string(nsec,format='(I0)')+'? 0,1,...'
                read, tc
            endwhile
            etr = [t0[secs[tc,0]],t0[secs[tc,1]]]
        endelse
    endelse
    ; implement t1, t2.
    if n_elements(t1) ne 0 then etr[0] = t1 else t1 = etr[0]
    if n_elements(t2) ne 0 then etr[1] = t2 else t2 = etr[1]
    print, 'E time range: ', time_string(etr[0]), ' to ', time_string(etr[1])
    
    ; get B field.
    ucla_mag_despin
    
    ; truncate b to [t1,t2].
    get_data, 'dB_fac_v', data = tmp
    idx = where(tmp.x ge t1 and tmp.x le t2)
    store_data, 'dB_fac_v', data = {x:tmp.x[idx],y:tmp.y[idx,*]}
    
    ; get model b in gei.
    get_data, 'B_model', data = tmp
    idx = where(tmp.x ge t1 and tmp.x le t2)
    store_data, 'B_model', data = {x:tmp.x[idx],y:tmp.y[idx,*]}, $
        limits = {ytitle:'B model GEI!C(nT)', $
        colors:[6,4,2], labels:['x','y','z']}

    ; prepare filename.
    if n_elements(fn) eq 0 then begin
        get_data, 'ORBIT', data = tmp
        orbstr = string(tmp.y[0], format = '(I05)')
        tstr = time_string(tmp.x[0], format = 2, precision = -3)
        tstr = strmid(tstr,0,4)+'_'+strmid(tstr,4,4)
        fn = '~/fa_sdt_fld_'+tstr+'_'+orbstr
    endif

    ; get orbit data.
    get_fa_orbit, t1, t2
    get_data, 'ILAT', data = tmp
    store_data, 'ilat', data = tmp, limits = {ytitle:'ILat'}
    get_data, 'MLT', data = tmp
    store_data, 'mlt', data = tmp, limits = {ytitle:'MLT'}
    get_data, 'ALT', data = tmp
    store_data, 'alt', data = tmp, limits = {ytitle:'ALT'}
    store_data, 'dis', data = {x:tmp.x, y:tmp.y/6356.7523+1}, $
        limits = {ytitle:'Dist (Re)'}
    get_data, 'ORBIT', data = tmp & orbstr = string(tmp.y[0],format='(I05)')
    vars = ['ORBIT','ALT','ILAT','ILNG','MLT']
    for i = 0, n_elements(vars)-1 do store_data, vars[i], /delete

    ; get E field.
    ; interface {out: 'E_ALONG_V', 'E_NEAR_B'}
    fa_fields_despin
    get_data, 'E_ALONG_V', data = tmp
    idx = where(tmp.x ge etr[0] and tmp.x le etr[1])
    dum = where(finite(tmp.y[idx], /nan) eq 0, cnt)	; 1 means nan.
    if cnt eq 0 then begin    ; load v158.
        printf, -1, 'no default E field ...'
        fa_fields_despin, /use_v158
        get_data, 'E_ALONG_V', data = tmp
        idx = where(tmp.x ge etr[0] and tmp.x le etr[1])
        dum = where(finite(tmp.y[idx], /nan) eq 0, cnt)	; 1 means nan.
        if cnt eq 0 then message, 'no v158 E field ...'
    endif
    ; trim data to etr.
    vname = 'E_ALONG_V' & get_data, vname, data = tmp
    store_data, vname, data = {x:tmp.x[idx], y:tmp.y[idx]}
    vname = 'E_NEAR_B' & get_data, vname, data = tmp
    store_data, vname, data = {x:tmp.x[idx], y:tmp.y[idx]}
    tmp = ['EFIT_NEAR_B','EFIT_ALONG_V']
    if not keyword_set(nodelete) then store_data, tmp, /delete

    if n_elements(var0) eq 0 then $
        var0 = ['E_ALONG_V','E_NEAR_B','dB_fac_v','B_model', $
        'fa_pos','fa_vel','ilat','mlt','alt','dis']
    tplot_save, var0, filename = fn

    ; delete unsued vars.
    tmp = ['Bx_sp','By_sp','Bz_sp','Bx_sc','By_sc','Bz_sc', $
        'Bx_sp_sm','By_sp_sm','Bz_sp_sm','B_gei','B_sm', $
        'dB_sc','dB_gei','dB_sm','dB_fac','MAG_FLAGS', $
        'spin_freq','spin_phase','BX_DEL','BY_DEL','BZ_DEL', $
        'TW_ZX','TW_ZY','TW_YY','TW_YX','O_X','O_Y', $
        'BFOOT','LAT','LNG','FLAT','FLNG','ORBIT','ILNG', $
        'despun_to_gei','gei_to_sm','gei_to_fac','gei_to_fac_v', $
        'EFIT_NEAR_B','EFIT_ALONG_V']
    if not keyword_set(nodelete) then store_data, tmp, /delete
end
