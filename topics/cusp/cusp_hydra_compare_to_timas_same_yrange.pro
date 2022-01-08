; Compare Hydra to Timas for given event id.
; force same yrange.

; look for all available events before 1999.
ids = cusp_id('all')
idx = where(stregex(ids, '199[678]') ne -1)
ids = ids[idx]

rootdir = '/Volumes/GoogleDrive/My Drive/works/data/cusp'
outdir = shomedir()+'/cusp_hydra'
if file_test(outdir,/directory) eq 0 then file_mkdir, outdir
pre0 = 'po_'
unit = '(mW/m!U2!N)'

xticklen = -0.02
yticklen = -0.01
tplot_options, 'ystyle', 1
tplot_options, 'labflag', -1
tplot_options, 'xticklen', xticklen
tplot_options, 'yticklen', yticklen

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
        get_data, var0[i], uts, data
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
    
    yrange1 = [0.1,0.3]    
    yrange2 = [-0.1,0.1]
    yrange3 = [-0.3,-0.1]
    poss = sgcalcpos(nvar, ypad=1)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    half = 0.3
    linestyle = 1
    for i=0, nvar-1 do begin
        nouttick = 1
        tpos = sgcalcpos(4, position=poss[*,i], ypad=0)
        pos1 = tpos[*,0]
        pos2 = tpos[*,1] & pos2[1] = tpos[1,2]
        pos3 = tpos[*,3]
        get_data, vars[i], limits=options
        options, vars[i], 'ytickformat', '(A1)'
        options, vars[i], 'ystyle', 1
        options, vars[i], 'labels', ''
        options, vars[i], 'ytitle', ''
        options, vars[i], 'xstyle', 5
        options, vars[i], 'yticks', 1
        options, vars[i], 'yminor', 2
        options, vars[i], 'yrange', yrange1
        tplot, vars[i], trange=utr0, position=pos1, /noerase, nouttick=nouttick
        plots, pos1[[0,2]], pos1[3]+[0,0], /normal
        xyouts, pos1[0]-xchsz*1, pos1[3]-ychsz*half, /normal, alignment=1, sgnum2str(max(yrange1))
        
        options, vars[i], 'yticks', 2
        options, vars[i], 'labels', options.labels
        options, vars[i], 'yrange', yrange2
        tplot, vars[i], trange=utr0, position=pos2, /noerase, nouttick=nouttick
        xyouts, pos2[0]-xchsz*1, pos2[3]-ychsz*half, /normal, alignment=1, sgnum2str(max(yrange2))
        xyouts, pos2[0]-xchsz*1, pos2[1]-ychsz*half, /normal, alignment=1, sgnum2str(min(yrange2))
        xyouts, pos2[0]-xchsz*1, (pos2[1]+pos2[3])*0.5-ychsz*half, /normal, alignment=1, '0'
        plots, pos2[[0,2]], pos2[1]+[0,0], /normal, linestyle=linestyle
        plots, pos2[[0,2]], pos2[3]+[0,0], /normal, linestyle=linestyle
        plots, pos2[[0,2]], (pos2[1]+pos2[3])*0.5, /normal, linestyle=linestyle
        xyouts, pos2[0]-xchsz*5, (pos2[1]+pos2[3])*0.5, /normal, alignment=0.5, orientation=90, options.ytitle
        
        if i eq nvar-1 then nouttick = 0
        options, vars[i], 'yticks', 1
        options, vars[i], 'labels', ''
        options, vars[i], 'xstyle', 8
        options, vars[i], 'xticklen', xticklen*4
        options, vars[i], 'yrange', yrange3
        tplot, vars[i], trange=utr0, position=pos3, /noerase, nouttick=nouttick
        xyouts, pos3[0]-xchsz*1, pos3[1]-ychsz*half, /normal, alignment=1, sgnum2str(min(yrange3))
        
        xyouts, poss[0,i]-xchsz*10, poss[3,i]-ychsz*0.5, /normal, figlabs[i], alignment=0
    endfor
    sgclose
    
endforeach

end
