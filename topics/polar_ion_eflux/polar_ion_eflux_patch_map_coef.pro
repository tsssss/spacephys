
pro polar_ion_eflux_patch_map_coef, utr0

;---Constants.
    cns = -1
    re = 6378d & re1 = 1/re
    r0 = 100d*re1+1
    
    
;---Settings.
    par = 2
    posvars = ['Epoch','GSM_POS']
    tmsvars = ['utsec']


    fn = sdiskdir('Research')+'/sdata/polar/timas/'+$
        time_string(utr0[0],tformat='YYYY')+'/po_tim_moments_'+$
        time_string(utr0[0],tformat='YYYY_MMDD')+'.cdf'
    if file_test(fn) eq 0 then return
    tms = scdfread(fn, tmsvars, skt=skt)
    vnames = tms.name
    idx = where(vnames eq 'mapcoef', cnt)
    if cnt ne 0 then return

    ; add map_coef to cdf.
    pos = sread_polar_orbit(utr0, flag='def', vars=posvars, /local_only)
    if size(pos,/type) ne 8 then $
        pos = sread_polar_orbit(utr0, flag='pre', vars=posvars, /local_only)
    if size(pos,/type) ne 8 then begin
        printf, cns, 'no position data ...'
        return
    endif
    if size(pos,/type) ne 8 then return
    
    uts = *tms[0].value
    tuts = sfmepoch(pos.epoch, 'unix')
    r0gsm = sinterpol(pos.gsm_pos*re1, tuts, uts)
    ets = stoepoch(uts, 'unix')
    nrec = n_elements(uts)
    
    mapc = dblarr(nrec)
    b0gsm = dblarr(nrec,3)
    b1gsm = dblarr(nrec,3)
    r1gsm = dblarr(nrec,3)
    
    
    for i=0, nrec-1 do begin
        geopack_epoch, ets[i], year, mo, dy, hr, mi, sc, /breakdown_epoch
        geopack_recalc, year, mo, dy, hr, mi, sc, /date, tilt = tilt
        
        ; find in-situ B (IGRF as bg + t89 modification).
        geopack_igrf_gsm, r0gsm[i,0], r0gsm[i,1], r0gsm[i,2], bx, by, bz
        b0gsm[i,*] = [bx,by,bz]
        geopack_t89, par, r0gsm[i,0], r0gsm[i,1], r0gsm[i,2], bx, by, bz
        b0gsm[i,*] += [bx,by,bz]
        
        ; trace to find footpoint pos in gsm.
        dir = (r0gsm[i,2] ge 0)? -1: 1
        geopack_trace, r0gsm[i,0], r0gsm[i,1], r0gsm[i,2], $
            dir, par, xf, yf, zf, /t89, r0=r0, /refine, /ionosphere
        r1gsm[i,*] = [xf,yf,zf]
        
        ; footpoint B.
        geopack_igrf_gsm, xf, yf, zf, bx, by, bz
        b1gsm[i,*] = [bx,by,bz]
    endfor
    
    mapcoef = snorm(b1gsm)/snorm(b0gsm)
    ainfo = {$
        FIELDNAM:'T89 map coefficient',$
        DEPEND_0:'utsec',$
        VAR_TYPE:'support_data'}
    scdfwrite, fn, 'mapcoef', value=mapcoef, attribute=ainfo

    ainfo = {$
        FIELDNAM:'T89 footpoint in GSM in Re',$
        DEPEND_0:'utsec',$
        UNIT:'(Re)',$
        VAR_TYPE:'support_data'}
    scdfwrite, fn, 'foot_gsm', value=r1gsm, attribute=ainfo

    
;    store_data, 'mapcoef', uts, mapcoef
;    store_data, 'r0gsm', uts, r0gsm, limits={labels:['x','y','z'], colors:[6,4,2]}
;    store_data, 'r1gsm', uts, r1gsm, limits={labels:['x','y','z'], colors:[6,4,2]}
;    store_data, 'b0gsm', uts, b0gsm, limits={labels:['x','y','z'], colors:[6,4,2]}
;    store_data, 'b1gsm', uts, b1gsm, limits={labels:['x','y','z'], colors:[6,4,2]}
;
;    device, decomposed=0
;    loadct2, 43
;    
;    vars = ['r0gsm','r1gsm','b0gsm','b1gsm','mapcoef']
;    tplot, vars, trange=utr0
    

end



secofday = 86400d

utr0 = time_double(['1996-03-17','1997-12-31']) 
utr0 = time_double(['1998-01-01','1998-12-31'])
;utr0 = time_double(['1996-03-19','1996-03-20'])
nday = (utr0[1]-utr0[0])/secofday
ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
ut2s = ut1s+secofday

for i=0, nday-1 do begin
    tutr = [ut1s[i],ut2s[i]]
    polar_ion_eflux_patch_map_coef, tutr
endfor

end
