;+
; Plot keV H+.
;-

pro fig_case2_kev_h_flux

;---Load electron flux.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1

;---Settings.
    einfo = get_var_data(info_var)
    time_range_plot = time_double(['2014-08-28/10:00','2014-08-28/10:40'])

    probes = ['LANL-02A','LANL-04A','LANL-97A','1994-084','rbspb','g13','1991-080','tha','the','g15','thd','LANL-01A']
    nprobe = n_elements(probes)
    color_start = 80
    color_end = 250
    color_table = 40
    colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
    for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)
    colors[where(probes eq 'the')] = sgcolor('gold')    ; the yellow is too bright.


    the_probes = ['thd','LANL-01A','LANL-02A','LANL-04A','LANL-97A']
    plot_info = dictionary()
    plot_info['probes'] = the_probes
    plot_info['colors'] = []
    plot_info['energy_range'] = [70,500]
    foreach probe, the_probes do begin
        plot_info.colors = [plot_info.colors,colors[where(probes eq probe)]]
    endforeach

    foreach probe, plot_info.probes do begin
        get_data, probe+'_kev_h_flux', times, flux, energy_bins
        index = where_pro(energy_bins, plot_info.energy_range, count=nenergy_bin)
        flux = flux[*,index]
        energy_bins = energy_bins[index]

        if nenergy_bin gt 6 then begin
            flux = flux[*,0:*:2]
            energy_bins = energy_bins[0:*:2]
        endif

        yrange = alog10(minmax(flux))
        yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]

        var = probe+'_kev_h'
        store_data, var, times, flux, energy_bins
        add_setting, var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            yrange: yrange, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'H!U+!N flux'}
    endforeach


;---Plot.
    vars = plot_info.probes+'_kev_h'
    nvar = n_elements(vars)
    fig_xsize = 4
    aspect_ratio = 0.3
    fig_ysize = nvar*fig_xsize*aspect_ratio

    if keyword_set(test) then begin
        file = test
        magnify = 1.2
    endif else begin
        file = join_path([sparentdir(srootdir()),'plot','fig_case2_kev_h_flux.pdf'])
        magnify = 1
    endelse

    ; x-labels.
    options, vars, 'xticklen', -0.02
    tplot_options, 'version', 2

    ; y-labels.
    foreach var, vars do begin
        get_data, var, limits=lim
        yrange = lim.yrange
        ytickv = make_bins(alog10(yrange),1)
        ytickn = '10!U'+string(ytickv,format='(I0)')
        index = where(ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '10'
        index = where(ytickv eq 0, count)
        if count ne 0 then ytickn[index] = '1'
        index = where(ytickv eq -1, count)
        if count ne 0 then ytickn[index] = '0.1'

        yticks = n_elements(ytickv)-1
        if yticks ge 4 then ytickn[1:*:2] = ' '

        options, var, 'ytickv', 10^ytickv
        options, var, 'ytickname', ytickn
        options, var, 'yticks', yticks
        options, var, 'yminor', 10
    endforeach

    ; color.
    zrange = plot_info.energy_range
    foreach var, vars do begin
        get_data, var, xxs, yys, zzs, limits=lim
        ct = lim.color_table
        index_colors = bytscl(zzs, min=zrange[0], max=zrange[1], top=color_end-color_start)+color_start
        ncolor = n_elements(index_colors)
        colors = fltarr(ncolor)
        for ii=0,ncolor-1 do colors[ii] = sgcolor(index_colors[ii],ct=ct)
        options, var, 'colors', reverse(colors)
    endforeach



    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    poss = sgcalcpos(nvar, margins=[9,4,8,1], xchsz=xchsz, ychsz=ychsz)
    tplot, vars, position=poss, trange=time_range_plot

    fig_labels = letters(nvar)
    full_ysz = 0.8
    foreach probe, plot_info.probes, ii do begin
        tpos = poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        label = fig_labels[ii]+'.  '+strupcase(probe)
        xyouts, tx,ty, /normal, label
    endforeach


    if keyword_set(test) then stop
    sgclose

end
