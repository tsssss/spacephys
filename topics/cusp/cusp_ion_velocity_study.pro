;+
; study the ion parallel velocity vs perpendicular velocity.
; 
; also check the possibility to perform a quick statistical study
; on the cusp velocity distribution, which resembles the Cluster one.
; 
; also check if pressure can be obtained.
;-
if tnames('cusp_scidat') eq '' then begin
    ids = cusp_id('south_imf')
    fns = shomedir()+'/Google Drive/works/data/cusp/'+ids+'_all_data.tplot'
    nfn = n_elements(fns)
    tinfo = cusp_info_struct()
    infos = replicate(tinfo, nfn)
    store_data, '*', /delete
    
    for i = 0, nfn-1 do begin
        printf, -1, fns[i]
        tplot_restore, filename = fns[i]
        get_data, 'scidat', tmp, info
        infos[i] = info
    endfor
    
    store_data, 'cusp_scidat', 0, infos
endif
get_data, 'cusp_scidat', 0, infos
ninfo = n_elements(infos)
vparas = dblarr(ninfo)
vperps = dblarr(ninfo)
dens = dblarr(ninfo)
diss = dblarr(ninfo)

for i = 0, ninfo-1 do begin
    tinfo = infos[i].polar.cusp
    cusptr = [tinfo.entry.ut, tinfo.exit.ut]
    hyd = sread_polar_hydra_moment(cusptr)
    
    etr = stoepoch(cusptr, 'unix')
    idx = where(hyd.epoch ge etr[0] and hyd.epoch le etr[1], cnt)
    if cnt eq 0 then continue

    ; times, density, velocity in GSM.
    uts = sfmepoch(hyd.epoch[idx], 'unix')
    eden = hyd.density_ele[idx]
    vgsm = hyd.bulk_velocity_ion[idx,*]

    ; load b gsm, pos gsm.
    mfe = sread_polar_mfe(cusptr)
    tmp = sfmepoch(mfe.epoch,'unix')
    bgsm = sinterpol(mfe.b_gsm,tmp, uts)
    rgsm = sinterpol(mfe.pos_gsm,tmp, uts)

    vpara = sdot(vgsm,sunitvec(bgsm))
    vperp = vgsm[*,1]*cos(atan(bgsm[*,1],bgsm[*,0]))
    vperp = sqrt(total(vgsm^2,2)-vpara^2-vperp^2)

    hem = infos[i].polar.hem
    idx = (hem eq -1)? where(vpara gt 0, cnt): where(vpara lt 0, cnt)
    vpara*= hem
    if cnt eq 0 then begin
        vparas[i] = !values.d_nan
        vperps[i] = !values.d_nan
    endif else begin
        vparas[i] = mean(vpara[idx],/nan)
        vperps[i] = mean(vperp[idx],/nan)        
    endelse

    
;    vparas[i] = vpara[-1]
;    vperps[i] = vperp[-1]
    diss[i] = 0.5*(infos[i].polar.cusp.entry.dis+infos[i].polar.cusp.exit.dis)
    dens[i] = mean(eden,/nan)
endfor


stop

end
