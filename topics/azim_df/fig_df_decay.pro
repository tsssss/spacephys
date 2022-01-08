;+
; Plot DF decay in MLT for all DFs:
;-

test = 0
    if n_elements(project) eq 0 then project = azim_df_load_project()


    time_range = time_double(['2007-01-01','2017-10-01'])
    all_probes = project.all_probes
    ramps = list()
    foreach probe, all_probes do begin
        ramps.add, azim_df_read_ramp(time_range, probe=probe, project=project), /extract
    endforeach
    
    target_nramp = 40000.
    nramp = ramps.length
    index = round(findgen(target_nramp)*nramp/(target_nramp-1))
    ramps = ramps[index]

    nramp = ramps.length
    df_widths = fltarr(nramp)
    df_heights = fltarr(nramp)
    df_rxys = fltarr(nramp)
    df_mlts = fltarr(nramp)
    foreach df, ramps, ii do begin
        df_widths[ii] = df.width
        df_heights[ii] = df.height
        df_rxys[ii] = df.obs_rxy
        df_mlts[ii] = df.obs_mlt
    endforeach

;    regions = list()
;    regions.add, dictionary($
;        'name', 'pre_midn', $
;        'mlt_range', [-9,0])
;    regions.add, dictionary($
;        'name', 'post_midn', $
;        'mlt_range', [0,9])
;;    regions.add, dictionary($
;;        'name', 'around_midn', $
;;        'mlt_range', [-3,3])
;
;    file = join_path([project.data_dir,'azim_df_roi_ramps.txt'])
;    if file_test(file) eq 0 then begin
;        ftouch, file
;        
;        candidates = list()
;        foreach region, regions do begin
;            search_setting = dictionary()
;            search_step = 'search_roi'
;            search_info = dictionary($
;                'file_suffix', 'azim_df_search_'+region.name+'_'+search_step+'.txt')
;            if search_step eq 'search_roi' then begin
;                search_info['name'] = region.name
;                search_info['mlt_range'] = region.mlt_range
;            endif
;            search_setting[search_step] = search_info
;            candidates.add, azim_df_search_roi(search_setting, project=project), /extract
;        endforeach
;        
;        
;        ncandidate = candidates.length
;        times = dblarr(ncandidate)
;        foreach candidate, candidates, ii do times[ii] = candidate.time_range[0]
;        index = sort(times)
;        candidates = candidates[index]
;        foreach candidate, candidates, ii do candidates[ii].id = ii+1
;        
;        
;        foreach candidate, candidates do begin
;            lprmsg, 'Processing candidate '+string(candidate.id,format='(I0)')+' ...'
;            ; Collect DFs.
;            probe_list = candidate.probe_list
;            time_range_list = candidate.time_range_list
;            foreach sector_time_range, time_range_list, sector_id do begin
;                probes = probe_list[sector_id]
;                foreach probe, probes do begin
;                    the_ramps = azim_df_read_ramp(sector_time_range, probe=probe, project=project)
;                    foreach ramp, the_ramps, id do azim_df_vertex_write, ramp, filename=file
;                endforeach
;            endforeach
;        endforeach
;    endif
;
;    lines = read_all_lines(file)
;    nramp = n_elements(lines)
;    target_nramp = 40000.
;    index = round(findgen(target_nramp)*nramp/(target_nramp-1))
;    lines = lines[index]
;    nramp = n_elements(lines)
;    df_widths = fltarr(nramp)
;    df_heights = fltarr(nramp)
;    df_rxys = fltarr(nramp)
;    df_mlts = fltarr(nramp)
;    
;    foreach line, lines, ii do begin
;        df = azim_df_vertex_read(line)
;        df_widths[ii] = df.width
;        df_heights[ii] = df.height
;        df_rxys[ii] = df.obs_rxy
;        df_mlts[ii] = df.obs_mlt
;    endforeach
    
    


;if keyword_set(test) then stop


    file = join_path([project.plot_dir, 'fig_df_decay.pdf'])
    if keyword_set(test) then file = 0
    fig_xsize = 3
    fig_ysize = 3
    magnify = keyword_set(test)? 2:1
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch, magnify=magnify

    margins = [8,4,2,1]
    tpos = sgcalcpos(1,1, margins=margins, xchsz=xchsz, ychsz=ychsz)

    xticklen = -0.02
    yticklen = -0.02

    mlt_range = [-1,1]*12
    mlt_step = 6
    mlt_minor = 3
    mlt_title = 'MLT (hr)'
    height_range = [1e-2,8e2]
    height_title = tex2str('Delta')+tex2str('theta')+' (deg)'
    height_color = sgcolor('wheat')
    label_size = constant('label_size')



;---Panel a. Height vs MLT.
    xrange = mlt_range
    xstep = mlt_step
    xminor = mlt_minor
    xtickv = make_bins(xrange,xstep)
    xticks = n_elements(xtickv)-1
    xlog = 0
    xtitle = mlt_title
    xtickformat = ''

    yrange = height_range
    ylog = 1
    ytitle = height_title
    ytickformat = ''

    bin_psym = 6
    bin_symsize = label_size*0.5
    symbol_color = height_color


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=xlog, $
        ystyle=5, yrange=yrange, ylog=ylog, $
        /nodata, /noerase, position=tpos

    xxs = df_mlts
    yys = df_heights
;    index = where(df_rxys le 8)
;    yys = df_heights[index]
;    xxs = df_mlts[index]

    oplot, xxs, yys, psym=1, symsize=bin_symsize, color=symbol_color
    ; Do binning.
    bins = make_bins(xrange,xstep/xminor)
    nbin = n_elements(bins)-1
    levels = smkarthm(0.1,0.9, 0.4, 'dx')
    ;levels = smkarthm(0.25,0.75, 0.25, 'dx')
    nlevel = n_elements(levels)
    bin_xxs = (bins[1:nbin]+bins[0:nbin-1])*0.5
    bin_yys = fltarr(nbin,nlevel)+!values.f_nan
    min_count = 10
    for ii=0, nbin-1 do begin
        index = lazy_where(xxs, '[]', bins[ii:ii+1], count=count)
        if count le min_count then continue
        the_data = yys[index]
        sorted_data = the_data[sort(the_data)]
        ndata = n_elements(sorted_data)
        bin_yys[ii,*] = sorted_data[ndata*levels]
    endfor

    colors = smkarthm(100,250,nlevel, 'n')
    level_ct = 49
    for jj=0, nlevel-1 do colors[jj] = sgcolor(colors[jj], ct=level_ct)
    fit_linestyle = 0
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*label_size
    msg = 'Percentile contours fitted by exp(-z!U2!N/2)'
    xyouts, tx,ty,/normal, msg, charsize=label_size
    for jj=0, nlevel-1 do begin
        color = colors[jj]
        plots, bin_xxs, bin_yys[*,jj], psym=-bin_psym, symsize=bin_symsize, color=color
        index = where(finite(bin_yys[*,jj]))
        fit_xxs = bin_xxs[index]
        fit_yys = gaussfit(fit_xxs, bin_yys[index,jj], coefs, nterm=3)
        oplot, fit_xxs, fit_yys, linestyle=fit_linestyle, color=color

    ;---Label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]+ychsz*(jj-nlevel-1)*label_size
        msg = string(levels[jj]*100,format='(I0)')+'%  z = (MLT'
        if coefs[1] ge 0 then msg += '-'+sgnum2str(coefs[1],ndec=1) else msg += '+'+sgnum2str(-coefs[1],ndec=1)
        msg += ')/'+sgnum2str(abs(coefs[2]),ndec=1)

        xyouts, tx,ty,/normal, msg, color=color, charsize=label_size
    endfor

    tmp_color = sgcolor('red')
    tmp_size = 0.4
    settings = azim_df_search_event_settings()
    setting = (settings[0]).search_large_df
    theta_midn = 8
    theta_scale = project.scale_width
    fit_yys = exp(-(fit_xxs/theta_scale)^2*0.5)*theta_midn
    plots, fit_xxs, fit_yys, linestyle=1, color=tmp_color
    tmp = convert_coord(0,theta_midn, /data, /to_normal)
    tx = tmp[0]-xchsz*0
    ty = tmp[1]+ychsz*0.5
    msg = 'Use DF > '+sgnum2str(theta_midn)+$
        ' exp(-z!U2!N/2),!C where z=MLT/'+sgnum2str(theta_scale)
    xyouts, tx,ty,/normal, msg, color=tmp_color, charsize=tmp_size, alignment=0.5


    ; Plot axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xlog=xlog, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ylog=ylog, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos


    if keyword_set(test) then stop
    sgclose


end
