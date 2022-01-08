; Compare Hydra to Timas for given event id.

; look for all available events before 1999.
ids = cusp_id('all')
idx = where(stregex(ids, '199[678]') ne -1)
ids = ids[idx]

rootdir = '/Volumes/GoogleDrive/My Drive/works/data/cusp'
outdir = shomedir()+'/cusp_hydra'
if file_test(outdir,/directory) eq 0 then file_mkdir, outdir
pre0 = 'po_'
unit = '(mW/m!U2!N)'

tplot_options, 'ystyle', 1

foreach id, ids do begin
    ; load first order data.
    datafiles = rootdir+'/'+id+['_1st_order_data','_scinfo']+'.tplot'
    foreach file, datafiles do tplot_restore, filename = file
    
    get_data, 'scinfo', tmp, info
    utr0 = info.polar.plot_time

    timas = sread_polar_timas_sheng_moment(utr0, vars=['ut_sec',['h','o','he_1','he_2']+'_energy_flux'])
    idx = where(timas.ut_sec ge utr0[0] and timas.ut_sec le utr0[1], nrec)
    store_data, pre0+'h_eflux', timas.ut_sec[idx], timas.h_energy_flux[idx], limits = {labels:'Timas H+'}
    store_data, pre0+'o_eflux', timas.ut_sec[idx], timas.o_energy_flux[idx], limits = {labels:'Timas O+'}
    store_data, pre0+'he_1_eflux', timas.ut_sec[idx], timas.he_1_energy_flux[idx], limits = {labels:'Timas He+'}
    store_data, pre0+'he_2_eflux', timas.ut_sec[idx], timas.he_2_energy_flux[idx], limits = {labels:'Timas He++'}
    timas_eflux = fltarr(nrec)
    
    ; add up selected efluxes.
    var0 = pre0+['h','o','he_1','he_2']+'_eflux'
    nvar = n_elements(var0)
    for i=0, nvar-1 do begin
        get_data, vars[i], uts, data
        timas_eflux += data
    endfor
    store_data, pre0+'timas_eflux', uts, data
    
    get_data, pre0+'ion_keflux', tuts, hydra_eflux
    store_data, pre0+'eflux', uts, [[interpol(hydra_eflux,tuts,uts)],[timas_eflux]], $
        limits = {colors:[0,6], labels:['Hydra','Timas']}
    
    vars = [var0,pre0+'eflux']
    options, vars, 'ytitle', unit
    options, vars, 'yrange', minmax([timas_eflux,hydra_eflux])
    options, vars, 'yrange', sg_autolim([timas_eflux,hydra_eflux])
    
    ofn = 0
    ofn = outdir+'/cusp_hydra_compare_to_timas_'+id+'.pdf'
    sgopen, ofn, xsize=5, ysize=6, /inch
    device, decomposed=0
    loadct2, 43
    
    label0 = ['a','b','c','d','e','f','g']
    nvar = n_elements(vars)
    figlabs = label0[0:nvar-1]+'.'
    
    tplot, vars, trange=utr0, title='Polar Ion energy flux, Hydra vs Timas', figlabel=figlabs
    sgclose
    
endforeach

end
