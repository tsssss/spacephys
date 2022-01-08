
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

    
    utr0 = time_double(['2014-08-28/09:30','2014-08-28/12:00'])     ; time to load B field.
    utr1 = time_double(['2014-08-28/10:00','2014-08-28/10:30'])     ; time to load asf.
    mltrng = [-1,6]
    
    
;---Load particle injection data.

    scs = 'rbb'
    foreach tsc, scs do begin
    endforeach


    scs = ['13','15']
    labels = ['40','75','150','275','475']+'keV'
    nlabel = n_elements(labels)
    foreach tsc, scs do begin
        pre0 = 'g'+tsc+'_'
        goes_load_data, trange=utr0, datatype='maged', probes=tsc, /noephem
        vars = pre0+'maged_'+labels+'_dtc_uncor_flux'
        get_data, vars[0], uts
        nrec = n_elements(uts)
        dat = dblarr(nrec,nlabel)
        for i=0, nlabel-1 do begin
            get_data, vars[i], tuts, tdat
            dat[*,i] = interpol(tdat[*,0], tuts, uts, /nan)
        endfor
        store_data, pre0+'ele_flux', uts, dat, limits=$
            {ytitle:'(#/cm!U2!N-s-sr-keV)', labels:labels, colors:findgen(nlabel), ylog:1}
    endforeach
    
    sgopen, plotdir+'/fig_drift_echo.pdf', xsize=10, ysize=5, /inch
    device, decomposed=0
    loadct2, 43
    
    
    tplot, ['g13_ele_flux','g15_ele_flux'], figlab = ['GOES-13','GOES-15'], /novtitle
    sgclose
    

end
