
;---Settings.
    id = '2014_0828_10'
    rootdir = sparentdir(srootdir())
    datadir = rootdir+'/data' & if file_test(datadir,/directory) eq 0 then file_mkdir, datadir
    plotdir = rootdir+'/plot' & if file_test(plotdir,/directory) eq 0 then file_mkdir, plotdir
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    deg = 180d/!dpi
    rad = !dpi/180d
    cns = -1

    
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero'
    tplot_options, 'ystyle', 1
    tplot_options, 'yticklen', -0.01
    tplot_options, 'xticklen', -0.03

    
    utr0 = time_double(['2014-08-28/09:50','2014-08-28/11:00'])     ; time to load B field.
    utr1 = time_double(['2014-08-28/10:00','2014-08-28/10:30'])     ; time to load asf.
    
    
    scs = ['thd','g15','the','tha','rbb','g13'] ; ordered according to mlt.
    nsc = n_elements(scs)
    pres = scs+'_'
    scsyms = [6,6,6,6,6,6,6]
    scolors = sgcolor(['red','tomato','orange','gold','lime','cyan','deep_sky_blue'])


    _2014_0828_load_data
    
    
    ; ewogram.
    ewo = 0
if keyword_set(ewo) then begin
    get_data, 'asf_mos', uts, mos, pxidx
    get_data, 'asf_info', 0, asinfo
    mlts = asinfo.mlts
    mlats = asinfo.mlats
    imgsz = asinfo.imgsz
    
    utr2 = time_double(['2014-08-28/10:05','2014-08-28/10:20'])
    idx = where(uts ge utr2[0] and uts le utr2[1])
    uts = uts[idx]
    mos = mos[idx,*]
    
    nrec = n_elements(uts)
    latvals = [60d,70]
    latvals = [63d,68]
    ;latvals = [68d,72]
    mltrng = [-1,2]
    mincnt = 20
    dmltval = 0.05d  ; deg.
    mltvals = smkarthm(mltrng[0],mltrng[1],dmltval,'dx')
    nmltval = n_elements(mltvals)-1
    mlts = mlts/15      ; convert to hr.
    ewos = fltarr(nrec,nmltval)
    for i=0, nrec-1 do begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        
        for j=0, nmltval-1 do begin
            idx = where(mlts ge mltvals[j] and mlts le mltvals[j+1] $
                and mlats ge latvals[0] and mlats le latvals[1] and timg ge mincnt, cnt)
            tvals = (timg[idx])[*]
            ewos[i,j] = total(tvals)/cnt
            ;                ewos[i,j] = median(tvals)
        endfor
    endfor
    
    store_data, 'ewo', uts, ewos, mltvals[0:nmltval-2], limits=$
        {ytitle:'MLT (hr)', ztitle:'Photon Count', spec:1, no_interp:1, zrange:[mincnt,200], yrange:mltrng, yticks:2, yminor:7}
endif

    mltrng = [-1,6]
    options, 'ewo', 'yrange', mltrng
    
    figfn = plotdir+'/fig_df_prop.pdf'
    ;figfn = 0
    
    vars = [pres+'bmag_tmp','ewo']
    nvar = n_elements(vars)
        
    sgopen, figfn, xsize=5, ysize=6, /inch
        
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    poss = sgcalcpos(nvar, lmargin=9, rmargin=8, tmargin=2, bmargin=8)
    poss[1,nvar-1] = 4*ychsz
        
    device, decomposed=0
    loadct2, 43
    tplot, vars, trange=utr0, position=poss, /novtitle, figlab=figlab

    device, decomposed=1
    thick = (size(ofn,/type) eq 7)? 4:2
    
    pos0 = poss[*,nvar-1]
    get_data, 'df_info', dfuts, dfmlts & dfmlts *= 15
    get_data, 'ewo', limits=lim0
    yrng = lim0.yrange
    
    for i=0, nsc-1 do begin
        get_data, vars[i], limits=lim
        plot, utr0, lim.yrange, position=poss[*,i], xstyle=5, ystyle=5, /nodata, /noerase
        plots, dfuts[i]+[0,0], lim.yrange, color=scolors[i]
    endfor
    
    plot, utr0, yrng, /noerase, /nodata, position=pos0, xstyle=5, ystyle=5
    npair = (size(ut0s,/dimensions))[0]
    res = linfit(dfuts, dfmlts)
    for i=0, nsc-1 do begin
        oplot, utr0, utr0*res[1]+res[0], linestyle=1, color=white, thick=thick*0.5
        plots, dfuts[i], dfmlts[i], psym=1, symsize=0.5, color=scolors[i], thick=thick
        tx = (yrng[0]-res[0])/res[1]
        ty = yrng[0]
        tmp = convert_coord([tx,ty,0], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.6
    endfor

    xyouts, pos0[2]-xchsz, pos0[1]+ychsz*0.2, /normal, alignment=1, $
        'DF eastward @ '+string(res[1]*3600,format='(F4.2)')+' MLT/hr'
    sgclose


end
