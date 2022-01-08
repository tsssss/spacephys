;+
; Study the different between E and Edot0.
;-

test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project['pflux_calc_setting']
    probes = ['a','b']
    plot_dir = join_path([project.plot_dir,'pflux_survey','figures'])
    orbit_time_step = 60.

    de_vars = ['de','dedot0']
    fac_labels = ['earth','west','out']
    data_file = join_path([project.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    if file_test(data_file) eq 1 then tplot_restore, filename=data_file


;---Pick out times when Edot0 are not NaNs.
    short_name = 'E'
    de_unit = 'mV/m'
    de_log = 0
    de_labels = ['earth','west','out']+'ward'
    de_range = [-1,1]*200
    xtickv = [-1,0,1]*150
    xticks = n_elements(xtickv)-1
    xminor = 3
    ytickv = xtickv
    yticks = xticks
    yminor = xminor
    ndim = 3
    fig_letters = letters(ndim)
    de_ticklen = -0.02

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        get_data, prefix+'dedot0_fac_interp', times, dedot0
        get_data, prefix+'de_fac_interp', times, de

        index = where(finite(dedot0[*,0]), ntime)
        times = times[index]
        de = de[index,*]
        dedot0 = dedot0[index,*]

;        ; Too many data points.
;        step = 5
;        times = times[0:*:step]
;        de = de[0:*:step,*]
;        dedot0 = dedot0[0:*:step,*]

        mlt = azim_df_calc_pseudo_mlt(get_var_data(prefix+'r_gsm', at=times))
        ct = 6  ; circular color table.
        ncolor = 24
        color_min = 5
        color_max = 245
        index_colors = smkarthm(color_min, color_max, ncolor, 'n')
        colors = index_colors
        foreach color, index_colors, ii do colors[ii] = sgcolor(color, ct=ct)
        mlt_min = -12
        mlt_max = 12
        dmlt = 24/ncolor
        mlt_bins = smkarthm(mlt_min+dmlt*0.5, mlt_max-dmlt*0.5, dmlt, 'dx')
        mlt_colors = fltarr(ntime)+colors[0]
        for ii=0, n_elements(mlt_bins)-2 do begin
            index = lazy_where(mlt, '[)', mlt_bins[ii:ii+1], count=count)
            if count eq 0 then continue
            mlt_colors[index] = colors[ii+1]
        endfor

        fig_xsize = 8
        fig_ysize = 3
        plot_file = join_path([project.plot_dir,'pflux_survey','figures','fig_de_vs_dot0_rbsp'+probe+'_vs_mlt.pdf'])
        if keyword_set(test) then plot_file = test
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize

        xpad = 1
        poss = sgcalcpos(1,3, xchsz=xchsz, ychsz=ychsz, xpad=xpad, rmargin=8, tmargin=2)

        ;---Colorbar.
        cbpos = poss[*,ndim-1]
        cbpos[0] = poss[2,ndim-1]+xchsz*xpad
        cbpos[2] = cbpos[0]+xchsz*1
        ztickv = [-6,0,6,12]
        zticks = n_elements(ztickv)
        ztitle = 'MLT (hr)'
        sgcolorbar, index_colors, ct=ct, position=cbpos, $
            zrange=[mlt_min, mlt_max]-dmlt*0.5, zticks=zticks, ztickv=ztickv, ztitle=ztitle


        for ii=0,ndim-1 do begin
            tpos = poss[*,ii]
            txs = de[*,ii]
            tys = dedot0[*,ii]
            xtitle = short_name+' ('+de_unit+')'
            ytitle = short_name+'!Ddot0!N ('+de_unit+')'
            if ii ne 0 then ytitle = ' '
            if ii ne 0 then ytickformat='(A1)' else ytickformat=''

            plot, de_range,de_range, /iso, $
                xstyle=1, xrange=de_range, xlog=de_log, xtitle=xtitle, xticklen=de_ticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
                ystyle=1, yrange=de_range, ylog=de_log, ytitle=ytitle, yticklen=de_ticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
                position=tpos, /noerase, ytickformat=ytickformat, /nodata

            plots, de_range, [0,0], linestyle=1, color=sgcolor('silver')
            plots, [0,0], de_range, linestyle=1, color=sgcolor('silver')
            plots, de_range, de_range, linestyle=1, color=sgcolor('silver')

            for jj=0,ncolor-1 do begin
                index = where(mlt_colors eq colors[jj], count)
                if count eq 0 then continue
                plots, txs[index], tys[index], color=colors[jj], psym=6, symsize=0.1
            endfor
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            fig_label = fig_letters[ii]+'. FAC '+short_name+' '+labels[ii]
            xyouts, tx,ty,/normal, fig_label
        endfor

        if keyword_set(test) then stop
        sgclose

    endforeach

end
