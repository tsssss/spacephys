;+
; test different mapping models.
;-

pro test_mapping, model

    probe = 'b'
    pre0 = 'rbsp'+probe+'_'
    utr = time_double(['2013-05-01/07:30','2013-05-01/07:45'])
    rgb = [6,4,2]
    re = 6378d & re1 = 1d/re
    deg = 180d/!dpi
    rad = !dpi/180d
    

; **** load some data.

    efwl3 = sread_rbsp_efw_l3(utr, probes = probe)
    uts = sfmepoch(efwl3.epoch,'unix',/epoch16)
    store_data, pre0+'pos_gse', uts, efwl3.pos_gse, $
        limits = {colors:rgb, labels:['x','y','z']}
    
; **** map footprint.
    ; prepare mapping.
    r0 = 1+110*re1  ; 110 km altitude.
    dir = -1        ; always north hem, b/c conjugate to thm_asi.
    sgeopack_par, utr, model, /delete  ; get tplot var <model>_par.
    t89 = 0 & t96 = 0 & t01 = 0 & ts04 = 0
    case model of
        't89': t89 = 1
        't96': t96 = 1
        't01': t01 = 1
        't04s': ts04 = 1
    endcase
    ; map pos to fpt, convert to mag.
    get_data, pre0+'pos_gse', data = tmp
    uts = tmp.x & ets = 1000D*uts+62167219200000D
    pos0 = tmp.y*re1 & pos1 = pos0      ; in re.
    ; interpolate par.
    if model ne 't89' then begin
        get_data, model+'_par', data = tmp
        pars = sinterpol(tmp.y, tmp.x, uts)
    endif
    ; loop for each time.
    for j = 0, n_elements(uts)-1 do begin
        ; set geopack.
        geopack_epoch, ets[j], yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001D, /date, tilt = tilt
        ; pos in gse, which is the mapping coord.
        x0 = pos0[j,0] & y0 = pos0[j,1] & z0 = pos0[j,2]
        ; convert from gse to gsm.
        geopack_conv_coord, x0, y0, z0, /from_gse, $
            x1, y1, z1, /to_gsm
        dir = (z1 gt 0)? -1: 1
        if model ne 't89' then par = reform(pars[j,*]) else par = 2
        geopack_trace, x1, y1, z1, dir, par, xf, yf, zf, $
            epoch = ets[j], tilt = tilt, /refine, /ionosphere, $
            t89 = t89, t96 = t96, t01 = t01, ts04 = ts04
        ; convert from gse to mag.
        geopack_conv_coord, xf, yf, zf, /from_gsm, $
            x1, y1, z1, /to_mag
        pos1[j,*] = [x1,y1,z1]
    endfor
    mlat = asin(pos1[*,2]*(1/r0))*deg
    mlon = atan(pos1[*,1],pos1[*,0])*deg
    store_data, pre0+'fpt_mlon_'+model, uts, mlon, $
        limits = {ytitle:'MLon/fpt (deg)'}
    store_data, pre0+'fpt_mlat_'+model, uts, mlat, $
        limits = {ytitle:'MLat/fpt (deg)'}
    store_data, pre0+'fpt_mlt_'+model, uts, slon2lt(mlon, stoepoch(uts,'unix'), /mag, /deg)/15, $
        limits = {ytitle:'MLT/fpt (hr)'}
    
    vars = pre0+'fpt_mlat_'+model
    options, vars, 'psym', -4
    tplot, vars, trange = utr
end

models = ['t89','t96','t01','t04s']
nmodel = n_elements(models)
;foreach model, models do test_mapping, model

vars = 'rbspb_fpt_mlat_'+models
tvar = 'rbspb_fpt_mlat'
tmp = []
for i = 0, nmodel-1 do begin
    get_data, vars[i], t0, dat
    tmp = [tmp,dat]
endfor
store_data, tvar, t0, reform(tmp, [n_elements(tmp)/nmodel,nmodel]), $
    limits = {colors:findgen(nmodel), labels:models, psym:-4, $
    ytitle: 'RBSP-B!CFootpoint!CMLat (deg)'}
    
vars = 'rbspb_fpt_mlon_'+models
tvar = 'rbspb_fpt_mlon'
tmp = []
for i = 0, nmodel-1 do begin
    get_data, vars[i], t0, dat
    tmp = [tmp,dat]
endfor
store_data, tvar, t0, reform(tmp, [n_elements(tmp)/nmodel,nmodel]), $
    limits = {colors:findgen(nmodel), labels:models, psym:-4, $
    ytitle: 'RBSP-B!CFootpoint!CMLon (deg)'}


ofn = shomedir()+'/test_models_trace_footpoint.pdf'
sgopen, ofn, xsize = 8, ysize = 6, /inch
device, decomposed = 0
loadct2, 43
vars = ['rbspb_fpt_mlat','rbspb_fpt_mlon']
tplot, vars

sgclose

end
