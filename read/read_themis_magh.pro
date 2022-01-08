;+
; Read themis ground B field (H component) into tplot.
; set varlist to add the produced tplot vars.
; set times to interpolate all vars to that time.
;-

pro read_themis_magh, utr0, sites=site0, newname=newname, $
    mlat_range = mlatrng, mlt_range = mltrng, $
    addto=varlist, times=times

    compile_opt idl2

;---Prepare sites. 
    if n_elements(site0) ne 0 then begin
        sites = site0
        infos = sread_thg_mag_location(sites)
    endif else begin
        ; filter sites according to mlat and mlt.
        infos = sread_thg_mag_location()
        if keyword_set(mlatrng) then begin
            idx = where(infos.mlat ge mlatrng[0] and infos.mlat le mlatrng[1], cnt)
            if cnt eq 0 then message, 'no site found in MLat range: '+string(mlatrng)+' ...'
            infos = infos[idx]
            sites = infos.id
        endif
        if keyword_set(mltrng) then begin
            ets = stoepoch(mean(utr0),'unix')
            mlts = slon2lt(infos.mlon, ets[0], /mag, /degree)
            idx = where(mlts ge mltrng[0] and mlts le mltrng[1], cnt)
            if cnt eq 0 then message, 'no site found in MLT range: '+string(mltrng)+' ...'
            infos = infos[idx]
            sites = infos.id
        endif
    endelse
    
;---Load data for given sites.
    store_data, 'thg_mag_*', /delete
    thm_load_gmag, site=sites, trange=utr0
    ; the sites that have data.
    vars = tnames('thg_mag_*')
    sites = strmid(vars, 8)
    nsite = n_elements(sites)
    sinfos = []
    oldsites = strlowcase(infos.id)
    for i=0, nsite-1 do sinfos = [sinfos, infos[where(oldsites eq sites[i])]]
    
    ; sort sites according to mlon (i.e., mlt).
    idx = sort(sinfos.mlon)
    sinfos = sinfos[idx]
    sites = strlowcase(sinfos.id)
    nsite = n_elements(sites)
    vars = 'thg_mag_'+sites
    mlts = slon2lt(sinfos.mlon, /mag, /deg, stoepoch(mean(utr0),'unix'))
    
    ; pick out h component, interpolate to uniform time.
    dt = 0.5
    ut0s = smkarthm(utr0[0], utr0[1], dt, 'dx')
    for i=0, nsite-1 do begin
        get_data, 'thg_mag_'+sites[i], uts, dat
        dt1 = 1d/sdatarate(uts)
        store_data, 'thg_magh_'+sites[i], ut0s, interpol(dat[*,0]-dat[0,0],uts,ut0s), $
            limits={ytitle:'(nT)', $
            labels:strupcase(sites[i])+' B!DH!N!C  '+$
            string(mlts[i],format='(F5.2)')+' MLT!C  '+$
            string(sinfos[i].mlat,format='(F4.1)')+' MLat'}
    endfor


    vars = 'thg_magh_'+sites
    stplot_merge, vars, newname = 'thg_magh', colors=indgen(n_elements(vars))+1
    
    get_data, 'thg_magh', ut0s, dat, limits=lim
    for i=0, nsite-1 do dat[*,i] = smooth(deriv(dat[*,i]),30/dt)
    store_data, 'thg_mag_dh', ut0s, dat, limits=lim
    
    stplot_split, 'thg_mag_dh', newname='thg_mag_dh_'+sites
    options, 'thg_mag_dh_'+sites, 'constant', 0
    
    
    vars = 'thg_'+['magh','mag_dh','mag_ddh']
    options, vars, 'constant', 0
    tplot, vars

    
    gmag_uts = time_double(['2014-08-28/10:07:52','2014-08-28/10:10:10', $
        '2014-08-28/10:11:24','2014-08-28/10:16:32','2014-08-28/10:17:08'])
    gmag_mlts = dblarr(nsite)
    for i=0, nsite-1 do begin
        gmag_mlts[i] = slon2lt(sinfos[i].mlon, /mag, /degree, stoepoch(gmag_uts[i],'unix'))
    endfor
    
    plot, gmag_uts-gmag_uts[0], gmag_mlts, psym=1, xstyle=1, ystyle=1, xrange=[-100,600],yrange=[-1,4]
    
    res = linfit(gmag_uts-gmag_uts[0], gmag_mlts)
    tx = [-100,600]
    plots, tx, tx*res[1]+res[0], linestyle=1
    
    stop


;---prepare var names.
    pre0 = 'th'+prb+'_'
    var0s = pre0+['efs_dot0_time','efs_dot0_gsm']   ; 3 sec resolution.
    var1s = idl_validname(var0s)

    
    ; make a uniform time.
    dr0 = 3d    ; sec.
    utr1 = utr0-(utr0 mod dr0)+[1,0]*dr0
    ut0s = smkarthm(utr1[0], utr1[1], dr0, 'dx')
    if keyword_set(times) then ut0s = times
    ; interpol to that uniform time.
    egsm = sinterpol(dat.(1), dat.(0), ut0s)

    myvars = pre0+['edot0_gsm']
    if not keyword_set(newname) then newname = 'edot0_gsm'
    store_data, pre0+newname, ut0s, egsm, $
        limits = {ytitle:'(mV/m)', colors:[6,4,2], labels:'GSM Edot0'+['x','y','z'], labflag:-1}
    if not keyword_set(varlist) then varlist=[]
    varlist = [varlist, myvars]
end

store_data, '*', /delete
utr = time_double(['2014-08-28/09:00','2014-08-28/12:00'])
sites = ['whit','fsim','fsmi','gill','fcc']
read_themis_magh, utr, site=sites, mlat_range = [63,68], mlt_range = [-6,8]

;utr = time_double(['2013-06-07/05:00','2013-06-07/06:00'])
;read_themis_magh, utr, mlat_range = [60,70], mlt_range = [-4,-2]
end

