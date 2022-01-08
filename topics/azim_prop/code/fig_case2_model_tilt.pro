;+
; A plot showing the measured field and model fields.
;-


;---Constants.
    rgb = [6,4,2]
    xyz = ['x','y','z']
    deg = 180d/!dpi
    rad = !dpi/180d
    tplot_options, 'ystyle', 1
    tplot_options, 'labflag', -1
    tplot_options, 'num_lab_min', 4
    strpm = '!M'+string(177b)+'!X'
    strsqrt = '!M'+string(214b)+'!X'
    spinrate = 12d
    label_size = 0.8


;---Settings.
    test = 0

    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
    einfo = get_var_data(info_var)
    time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    models = ['t89','t96','t01','t04s']
    model_colors = sgcolor(['red','green','blue','purple'])
    nmodel = n_elements(models)
    probes = ['thd','the']


    trace_dir = -1  ; to northern-hemisphere.
    foreach probe, probes do begin
        rvar = probe+'_r_gsm'
        foreach model, models do begin
            model_bvar = probe+'_bmod_gsm_'+model
            if tnames(model_bvar) eq '' then begin
                read_geopack_info, rvar, model=model, direction=trace_dir
            endif
            model_tiltvar = probe+'_bmod_sm_tilt_'+model
            if tnames(model_tiltvar) eq '' then begin
                get_data, model_bvar, times, bgsm
                bsm = cotran(bgsm, times, 'gsm2sm')
                tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
                store_data, model_tiltvar, times, tilt, limits={$
                    ytitle:'(deg)', labels:strupcase(model)}
            endif
        endforeach

        ; Combine variables.
        tilt_var = probe+'_sm_tilt'
        colors = [sgcolor('black'),model_colors]
        labels = ['Measured',strupcase(models)]
        if tnames(tilt_var) eq '' then begin
            get_data, probe+'_bmod_sm_tilt_'+models[0], times
            ntime = n_elements(times)
            data = fltarr(ntime,nmodel+1)
            bgsm = get_var_data(probe+'_b_gsm', at=times)
            bsm = cotran(bgsm, times, 'gsm2sm')
            data[*,0] = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
            for ii=0, nmodel-1 do data[*,ii+1] = get_var_data(probe+'_bmod_sm_tilt_'+models[ii], at=times)
            store_data, tilt_var, times, data, limits={ytitle:'(deg)', colors:colors, labels:labels}
        endif
    endforeach



;---Plot options.
    ticklen = -0.01
    tplot_options, 'xticklen', ticklen
    tplot_options, 'yticklen', ticklen*2
    tplot_options, 'labflag', -2


;---Plot.
    file = sparentdir(srootdir())+'/plot/fig_case2_model_vs_measured_b.pdf'
    if test eq 1 then begin
        file = 0
        magnify = 2
    endif else magnify = 1
    sgopen, file, xsize=4, ysize=3, magnify=magnify

    pres = ['thd','the']+'_'
    vars = pres+'sm_tilt'
    nvar = n_elements(vars)
    options, vars, 'yticks', 2
    options, vars, 'yminor', 5
    
    options, 'thd_sm_tilt', 'yrange', [20,80]
    options, 'the_sm_tilt', 'yrange', [0,60]

    poss = sgcalcpos(nvar, lmargin=11, tmargin=1, bmargin=4, rmargin=8)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    tplot, vars, trange=time_range, /novtitle, position=poss

    labels = ['a.','b.']+'TH-'+strupcase(strmid(pres,2,1))+'!C    B!DSM!N tilt'
    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        xyouts, xchsz*1.5, tpos[3]-ychsz*0.8, /normal, labels[ii], color=0
    endfor

    if keyword_set(test) then stop
    sgclose

end
