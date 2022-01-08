;+
; Generate EWOgram for case 2.
;-


;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1

;---Settings.
    einfo = get_var_data(info_var)
    time_range_plot = time_string(['2014-08-28/10:02','2014-08-28/10:17'])
    count_range = [0,300]

    aurora_color_table = 49


;---Calculate EWOgram.
    mlonimg_var = 'thg_mlonimg'
    ewo_var = mlonimg_var+'_ewo'

    ; Set settings.
    newo_info = 2
    ewo_infos = replicate({mlat_range:[0d,0], var_name:'', label:''},newo_info)
    ; the lower part, capture the initial eastward motion.
    ewo_infos[0].mlat_range = [63.5,65.5]
    ewo_infos[0].var_name = ewo_var+'_1'
    ; the upper part, capture the expansion.
    ewo_infos[1].mlat_range = [67,69]
    ewo_infos[1].var_name = ewo_var+'_2'

    foreach tinfo, ewo_infos, ii do begin
        if tnames(tinfo.var_name) ne '' then continue
        themis_read_mlonimg_ewo, mlonimg_var, to=tinfo.var_name, mlat_range=tinfo.mlat_range
        ewo_infos[ii].label = strjoin(string(tinfo.mlat_range,format='(F4.1)'),'-')
    endforeach


;---Prepare for plot.
    ofn = join_path([sparentdir(srootdir()),'plot','fig_case2_ewo.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=5, ysize=4

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    vars = ewo_infos.var_name
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, tmargin=4, bmargin=4, rmargin=3, lmargin=8)
    colorbar_pos = poss[[0,3,2,3],0]+[0,0.5,0,1]*ychsz

    options, vars, 'zrange', count_range
    options, vars, 'zposition', colorbar_pos
    options, vars, 'ztitle', 'Photon count (#)'
    options, vars, 'zhorizontal', 1
    options, vars, 'zticks', 3ß

    ytickv = [-100,-80,-60]
    yticks = n_elements(ytickv)-1
    options, vars, 'yrange', [-100,-60]
    options, vars, 'ytickv', ytickv
    options, vars, 'yticks', yticks
    options, vars, 'yminor', 4

    white = sgcolor('black')

    device, decomposed=0
    loadct2, aurora_color_table
    tplot, vars, trange=time_range_plot, position=poss, /novtitle
    device, decomposed=1

    ; Add labels.
    fig_labels = ['a','b']+'.'
    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1, /normal, fig_labels[ii]+' EWOgram in '+ewo_infos[ii].label+' deg MLat', color=white
    endfor

    sgclose
end
