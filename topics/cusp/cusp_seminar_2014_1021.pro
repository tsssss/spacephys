
pro cusp_seminar_2014_1021, eventid, lun = lun

    if n_elements(lun) eq 0 then lun = -1   ; console.
    
    rootdir = sdiskdir('Works')
    datroot = sdiskdir('Research')
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun.log'
    info = cusp_read_conjun_list(logfile, event = eventid)
    
    ; prepare.
    potr = info.polar.plot_time
    fatr = info.fast.plot_time
    potrcusp = info.polar.cusp_time
    fatrcusp = info.fast.cusp_time
    id = info.id
    
    ; **** polar keflux.
    tinfo = info.polar
    fn = rootdir+'/works/cusp/cusp list conjun/'+ $
        id+'/'+strmid(id,0,9)+'_ke_special.svg'
    tmp = tinfo.keflux    ; [ele,ion].
    polar_read_ke_flux, fn, potr, tmp[0]*[-1,1], tmp[1]*[-1,1]
    sdt = ssdtread(rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt') & pre = 'po_'
    ; ilat, mlt, dis.
    t1  = sdt.var.polarinvariantlatitude.depend_0
    mlat = sdt.var.polarmaglat.value
    store_data, pre+'mlat', data = {x:t1, y:mlat}, limits = {ytitle:'MLat'}
    ilat = sdt.var.polarinvariantlatitude.value
    tmp = where(mlat lt 0)
    if tmp[0] ne -1 then ilat[tmp] *= -1
    store_data, pre+'ilat', data = {x:t1, y:ilat}, limits = {ytitle:'ILat'}
    mlt  = sdt.var.polarmlt.value
    store_data, pre+'mlt', data = {x:t1, y:mlt}, limits = {ytitle:'MLT'}
    dis  = sdt.var.polarspcraftdist.value
    store_data, pre+'dis', data = {x:t1, y:dis}, limits = {ytitle:'Dist (Re)'}
    ; map and integrate.
    vars = ['po_'+['ele_keflux','ion_keflux']]
    for j = 0, n_elements(vars)-1 do begin
        smap2iono, vars[j], 'po_dis', newname = vars[j]+'_map'
        cusp_int_eflux, vars[j]+'_map', 'po_ilat', tr = tinfo.cusp_time
    endfor
    get_data, 'po_ilat', tmp, dat & hem = (dat[0] ge 0)*2-1 ; 1:N,-1:S.
    get_data, 'po_ion_keflux_map', tmp, dat
    poionyr = minmax(dat)*hem   ; [d,u].
    get_data, 'po_ele_keflux_map', tmp, dat
    poeleyr = minmax(dat)*hem   ; [d,u].
    
    ; **** fast keflux.
    tinfo = info.fast
    fn = rootdir+'/data/cusp/fa_sdt_esa_'+strmid(id,0,9)+'_'+ $
        string(tinfo.orbit,format='(I05)')+'.tplot'
    tplot_restore, filename = fn
    stplot_renew, 'ele_eflux', newname = 'fa_ele_keflux', /delete
    stplot_renew, 'ion_eflux', newname = 'fa_ion_keflux', /delete
    vars = ['ion_*','ele_*']
    store_data, vars, /delete
    ; map and integrate.
    vars = ['fa_'+['ele_keflux','ion_keflux']]
    for j = 0, n_elements(vars)-1 do begin
        smap2iono, vars[j], 'fa_dis', newname = vars[j]+'_map'
        cusp_int_eflux, vars[j]+'_map', 'fa_ilat', tr = tinfo.cusp_time
    endfor
    get_data, 'fa_ilat', tmp, dat & hem = (dat[0] ge 0)*2-1 ; 1:N,-1:S.
    get_data, 'fa_ion_keflux_map', tmp, dat
    faionyr = minmax(dat)*hem   ; [d,u].
    get_data, 'fa_ele_keflux_map', tmp, dat
    faeleyr = minmax(dat)*hem   ; [d,u].
    
    ; generate info line.
    trmin = min([potrcusp,fatrcusp], max = trmax)
    tmp = sread_omni(t = stoepoch([trmin,trmax],'unix'), rootdir = datroot)
    t0 = sfmepoch(tmp.epoch,'unix')
    idx = where(t0 ge trmin and t0 le trmax)
    store_data, 'ae', t0[idx], tmp.ae[idx]
    store_data, 'symh', t0[idx], tmp.symh[idx]
    store_data, 'imf', t0[idx], tmp.bgse[idx,*]
    aemax = max(tmp.ae[idx])
    symhmin = min(tmp.symh[idx])
    get_data, 'po_dis', tmp, dat
    podisyr = minmax(dat[where(tmp ge potrcusp[0] and tmp le potrcusp[1])])
    
    get_data, 'po_ion_keflux_map', tmp, tmp, poionint
    get_data, 'po_ele_keflux_map', tmp, tmp, poeleint
    get_data, 'fa_ion_keflux_map', tmp, tmp, faionint
    get_data, 'fa_ele_keflux_map', tmp, tmp, faeleint
    
    vars = ['po_'+['ele','ion']+'_keflux_map']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, dat
        store_data, vars[i]+'_abs', t0, abs(dat)
        cusp_int_eflux, vars[i]+'_abs', 'po_ilat', tr = potrcusp
    endfor
    vars = ['fa_'+['ele','ion']+'_keflux_map']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, dat
        store_data, vars[i]+'_abs', t0, abs(dat)
        cusp_int_eflux, vars[i]+'_abs', 'fa_ilat', tr = fatrcusp
    endfor
    get_data, 'po_ion_keflux_map_abs', tmp, tmp, poionabs
    get_data, 'po_ele_keflux_map_abs', tmp, tmp, poeleabs
;    get_data, 'fa_ion_keflux_map_abs', tmp, tmp, faionabs
;    get_data, 'fa_ele_keflux_map_abs', tmp, tmp, faeleabs
    
    ratio = string(poionint/poionabs,format='(F5.2)')
    tmp = poionint*hem & tmp1 = (tmp le 0)? 'U':'D'
    poionint = string(abs(tmp),format='(I4)')+tmp1
    tmp = poeleint*hem & tmp1 = (tmp le 0)? 'U':'D'
    poeleint = string(abs(tmp),format='(I4)')+tmp1
    tmp = faionint*hem & tmp1 = (tmp le 0)? 'U':'D'
    faionint = string(abs(tmp),format='(I4)')+tmp1
    tmp = faeleint*hem & tmp1 = (tmp le 0)? 'U':'D'
    faeleint = string(abs(tmp),format='(I4)')+tmp1
    poionabs = string(poionabs,format='(I5)')
    
    dt = double(mean(potrcusp)-mean(fatrcusp))/3600d  ; in hour.
    sep = '  '
    infostr = id+sep+ $
        string(mean(podisyr),format='(F3.1)')+sep+ $
        string(dt,format='(F4.1)')+sep+ $
        string(symhmin,format='(I4)')+sep+ $
        string(aemax,format='(I4)')+sep+ $
        string(abs(poionyr[0]),format='(F4.1)')+sep+ $  ; polar ion max up.
        string(abs(poionyr[1]),format='(F4.1)')+sep+ $  ; polar ion max down.
        poionint+sep+ $                                 ; polar ion integrated.
        poeleint+sep+ $                                 ; polar ele integrated.
        faionint+sep+ $                                 ; fast ion integrated.
        faeleint+sep+ $                                 ; fast ion integrated.
        poionabs+sep+ $                                 ; polar ion abs integrated.
        ratio                                           ; polar ele abs integrated.
    
    if lun eq -1 then begin
        printf, lun, '  event id    dis   dt   symh   ae   po ion map   po int map    fa int map    po int abs'
        printf, lun, '------------  ---  ----  ----  ----  ----------  ------------  ------------  ------------'
        printf, lun, 'yyyy_mmdd_hh  (R)  (hr)  (nT)  (nT)  umax  dmax   ion    ele    ion    ele    ion   ratio'
        printf, lun, '------------  ---  ----  ----  ----  ----  ----  -----  -----  -----  -----  -----  -----'
    endif
    printf, lun, infostr
end

events = ['1998_0914_19','1998_0925_05','1998_1001_02','1998_1001_03', $
    '1998_1002_16', '1998_1019_18','1998_1021_05','1998_1022_18', $
    '1999_1001_15','1999_1009_20','1999_1010_15','1999_1022_14','1999_1028_15']
;events = ['1998_1001_02']
openw, lun, shomedir()+'/cusp_seminar_summary.txt', /get_lun
for i = 0, n_elements(events)-1 do cusp_seminar_2014_1021, events[i], lun = lun
free_lun, lun
end