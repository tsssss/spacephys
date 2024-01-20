;+
; Plot ewo, pflux, dp.
;-

function fig_2015_0416_0800_rbsp_v02_rbsp_panel, rbsp_poss, event_info=event_info

;---Load data and settings.
    version = 'v02'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    ;time_range = event_info['time_range']
    dp_times = event_info['dp_times']
    
    power = 1d/4
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1


;---B spec.
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        probe = info['probe']

        ; B spec.
        b_fac_var = prefix+'b1_fac'
        b_vars = stplot_split(b_fac_var)
        the_b_var = b_vars[1]
        b_mor_var = stplot_mor_new(the_b_var, scale_info=scale_info)
        var = b_mor_var
        zrange = [1,1e6]
        log_ztickv = make_bins(alog10(zrange),1,inner=1)
        ztickv = 10d^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(log_ztickv eq 0, count)
        if count ne 0 then ztickn[index] = '1'
        index = where(log_ztickv eq 1, count)
        if count ne 0 then ztickn[index] = '10'
        ztickn[0:*:2] = ' '
        
        time_step = info['field_time_step']
        b0_window = info['b0_window']
        yrange = minmax(1d3/[b0_window*0.5,time_step*4])
        log_ytickv = [1,2,3]
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 9
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(log_ytickv eq 0, count)
        if count ne 0 then ytickn[index] = '1'
        index = where(log_ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '10'
        options, var, 'display_type', 'spec'
        options, var, 'ztitle', 'B!D'+fac_labels[1]+'!N (nT)!U2!N'
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'ztickname', ztickn
        options, var, 'zcharsize', label_size
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'ytickname', ytickn
        options, var, 'yminor', yminor
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', ytickv
        options, var, 'zcharsize', label_size
        
        
        ; Pflux.
        pf_spinfit_var = prefix+'pfdot0_fac_spinfit_phasef_map'
        var = pf_spinfit_var
        options, var, 'yrange', [-60,30]
        options, var, 'ytickv', [-60,-20,20]
        options, var, 'yticks', 2
        options, var, 'yminor', 2
        options, var, 'constant', [0]
        pf_survey_var = prefix+'pfdot0_fac_survey_map'
        if probe eq 'a' then pf_survey_var = prefix+'pfdot0_fac_v24_map'
        
        var = prefix+'pfdot0_fac_spinfit_map'
        copy_data, pf_spinfit_var, var
        var = prefix+'pfdot0_fac_survey_map'
        copy_data, pf_survey_var, var
        options, var, 'yrange', [-360,-60]
        options, var, 'ytickv', [-360,-180]
        options, var, 'yticks', 1
        options, var, 'yminor', 3
        
        vars = prefix+'pfdot0_fac_'+['spinfit','survey']+'_map'
        options, vars, 'ytitle', ' '
        options, vars, 'labels', [' ',' ',' ']
        
        var = prefix+'pfdot0_fac_survey_map'
        pfdot0 = get_var_data(var, times=times, limits=lim)
        pfpara = pfdot0[*,0]
        if probe eq 'a' then pfpara[*] = !values.f_nan
        var = prefix+'pfdot0_fac_survey_map_para'
        store_data, var, times, pfpara, limits=lim
        options, var, 'labels', ' '
        options, var, 'colors', sgcolor('light_salmon')
    endforeach


;---Plot.
    rbsp_letters = ['b','c','d','e']
    foreach key, sort_uniq(event_info.rbsp.keys()), rbsp_id do begin
        sc_info = (event_info.rbsp)[key]
        rbsp_letter = rbsp_letters[rbsp_id*2]
        rbsp_letter2 = rbsp_letters[rbsp_id*2+1]
        prefix = sc_info['prefix']
        probe = sc_info['probe']
        sc_color = sc_info['sc_color']
        
        plot_info = orderedhash()
        plot_info[prefix+'b1_fac_comp2_mor'] = dictionary($
            'fig_label', 'Wave', $
            'ypan', 0.8 )
        plot_info[prefix+'pfdot0_fac_spinfit_map'] = dictionary($
            'fig_label', 'S FAC', $
            'ypan', 0.6, $
            'ypad', 0. )
        plot_info[prefix+'pfdot0_fac_survey_map_para'] = dictionary($
            'fig_label', ' ', $
            'ypan', 0.4 )
        
        vars = 'rbspb_'+['b1_fac_comp2_mor','pfdot0_fac_spinfit_map','pfdot0_fac_survey_map_para']
        options, vars, 'ytickformat', '(A1)'
        options, vars, 'ytitle', ''
        
        vars = 'rbspa_'+['b1_fac_comp2_mor']
        options, vars, 'no_color_scale', 1
        
        
        plot_vars = (plot_info.keys()).toarray()
        nvar = n_elements(plot_vars)
        fig_labels = strarr(nvar)
        ypans = fltarr(nvar)
        ypads = fltarr(nvar)
        foreach key, plot_info.keys(), var_id do begin
            info = plot_info[key]
            fig_labels[var_id] = info['fig_label']
            if probe eq 'b' then fig_labels[var_id] = ' '
            if ~info.haskey('ypan') then info['ypan'] = 1.
            if ~info.haskey('ypad') then info['ypad'] = 0.4
            ypans[var_id] = info['ypan']
            ypads[var_id] = info['ypad']
        endforeach
        ypads = ypads[0:nvar-2]
        fig_labels = [[rbsp_letter,rbsp_letter2]+')',''];+fig_labels,'']

        my_pos = rbsp_poss[*,rbsp_id]
        poss = sgcalcpos(nvar, margins=margins, ypans=[ypans], ypad=ypads, position=my_pos)

        fig_size = double([!d.x_size,!d.y_size])
        tmp = get_charsize()
        xchsz = tmp[0]
        ychsz = tmp[1]
        uniform_ticklen = -ychsz*0.15*fig_size[1]
    
        ; ticklen.
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            options, plot_vars[pid], 'xticklen', xticklen
            options, plot_vars[pid], 'yticklen', yticklen
            var = plot_vars[pid]
            display_type = get_setting(var, 'display_type')
            if display_type eq 'spec' then begin
                zticklen = uniform_ticklen/xchsz*1/fig_size[0]
                options, var, 'zticklen', zticklen
            endif
        endfor
        symsize = 0.4
        
        ; Add labels for pflux.
        pf_var1 = prefix+'pfdot0_fac_spinfit_map'
        pf_var2 = prefix+'pfdot0_fac_survey_map'
        pid = where(plot_vars eq pf_var1, count)
        if count ne 0 then begin
            pos1 = poss[*,pid]
            pos2 = poss[*,pid+1]
            tpos = pos1
            tpos[1] = pos2[1]


            colors = constant('rgb')
            labels = 'S!D'+fac_labels+'!N'
            nlabel = n_elements(labels)
            dy = (total(tpos[[1,3]]*[-1,1])-ychsz*1)/(nlabel)
            label_poss = tpos[3]-findgen(nlabel+1)*dy


            if probe eq 'b' then begin
                for ii=0,nlabel-1 do begin
                    tx = tpos[2]+xchsz*0.5
                    ty = label_poss[ii]-ychsz*0.8
                    msg = labels[ii]
                    xyouts, tx,ty,msg, normal=1, color=colors[ii]
                endfor
                
                get_data, 'rbspb_pfdot0_fac_survey_map', times, pf_fac
                pf_para = pf_fac[*,0]
                the_color = colors[0]
                the_color = sgcolor('light_salmon')
                the_linestyle = 0

                get_data, pf_var1, limits=lim
                xrange = time_range
                yrange = lim.yrange
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, $
                    ystyle=5, yrange=yrange, $
                    position=pos1, nodata=1, noerase=1
                oplot, times, pf_para, color=the_color, linestyle=the_linestyle

                get_data, pf_var2, limits=lim
                xrange = time_range
                yrange = lim.yrange
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, $
                    ystyle=5, yrange=yrange, $
                    position=pos2, nodata=1, noerase=1
                oplot, times, pf_para, color=the_color, linestyle=the_linestyle
                msg = 'S!D||!N hires'
                tx = pos2[2]+xchsz*0.5
                ty = label_poss[-1]-ychsz*0.8
                xyouts, tx,ty,msg, normal=1, color=the_color
            endif
        endif

        tplot_options, 'version', 2
        tplot, plot_vars, position=poss, trange=time_range, noerase=1, nouttick=1
        
        ; manual unit.
        msg = '(mW/m!U2!N)'
        tpos = poss[*,1]
        tpos[1] = poss[1,2]
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4.5
        xyouts, tx,ty,msg, normal=1, alignment=0.5, orientation=90

        ; manual xticknames.
        tpos = poss[*,-1]
        xstep = 5*60d
        xrange = time_range
        xtickv = make_bins(xrange, xstep, inner=1)
        xtickn = time_string(xtickv, tformat='hhmm')
        xtickn[0] = time_string(xtickv[0], tformat='YYYY MTH DD          ')
        if probe eq 'b' then xtickn[0] = ' '
        
        plot, xrange, [0,1], $
            xstyle=5, ystyle=5, $
            nodata=1, noerase=1, position=tpos
        foreach msg, xtickn, id do begin
            tx = xtickv[id]
            tmp = convert_coord(tx,0, data=1, to_normal=1)
            tx = tmp[0]
            ty = tpos[1]-ychsz*1.1
            xyouts, tx,ty,normal=1, msg, alignment=0.5
        endforeach

        dp_tr = dp_times[rbsp_id]+[0,1.5]*60
        timebar, dp_tr, linestyle=1, color=sc_color
        label_yshift = -ychsz*0.6
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*5
            if probe eq 'b' then tx = tpos[0]-xchsz*2
            ty = tpos[3]+label_yshift
            xyouts, tx,ty, normal=1, fig_labels[pid]
        endfor
        if n_elements(bar_times) ne 0 then timebar, bar_times, color=sgcolor('red'), linestyle=1
        
        tpos = poss[*,0]
        msg = strupcase('RBSP-'+probe)
        tx = total(tpos[[0,2]])*0.5
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,normal=1, msg, color=sc_color, alignment=0.5

    endforeach

end




function fig_2015_0416_0800_rbsp_v02_mid_panel, my_pos, event_info=event_info

;---Load data and settings.
    version = 'v02'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    dp_times = event_info['dp_times']
    
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1


;---Plot ewo and dp.
    ; ewogram.
    ewo_var = 'thg_asf_ewo'
    mlat_range = [62,64]
    mlt_range = [-4,-2]
    ewo_zlog = 0
    if ewo_zlog eq 0 then begin
        zrange = [-1,1]*5e3
        ct = 70
        ztickv = [-1,0,1]*abs(max(zrange))
        zticks = n_elements(ztickv)-1
        zminor = 10
        ztickn = string(ztickv,format='(I0)')
    endif else begin
        zrange = [1e1,1e4]
        ct = 49
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(ztickn eq '10!U1', count)
        if count ne 0 then ztickn[index] = '10'
        index = where(ztickn eq '10!U0', count)
        if count ne 0 then ztickn[index] = '1'
    endelse
    
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    mlt_bins = lim.mlt_bins
    mlt_index = where_pro(mlt_bins, '[]', mlt_range, count=nmlt_bin)
    mlt_bins = mlt_bins[mlt_index]
    index = where(mlt_bins le 0, count)
    if count ne 0 then mlt_bins[index] += 24
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    mlat_bins = mlat_bins[mlat_index]
    ewo = total(data[*,mlt_index,mlat_index],3)/nmlat_bin
    if ewo_zlog eq 1 then ewo = ewo>0.01
    store_data, ewo_var, times, ewo, mlt_bins
    yrange = mlt_range+24
    ystep = 1
    ytickv = make_bins(yrange,ystep,inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 4
    add_setting, ewo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLT!C(h)', $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'zrange', zrange, $
        'zlog', ewo_zlog, $
        'ztickv', ztickv, $
        'zticks', zticks, $
        'zminor', zminor, $
        'ztickname', ztickn, $
        'color_table', ct )

    

;---Plot settings.
    plot_info = orderedhash()
    plot_info[ewo_var] = dictionary($
        'fig_label', 'Ewo', $
        'ypan', 1.5 )

    plot_vars = (plot_info.keys()).toarray()
    nvar = n_elements(plot_vars)
    fig_labels = 'a)'

    poss = my_pos
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    
    ; ticklen.
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
        var = plot_vars[pid]
        display_type = get_setting(var, 'display_type')
        if display_type eq 'spec' then begin
            zticklen = uniform_ticklen/xchsz*1/fig_size[0]
            options, var, 'zticklen', zticklen
        endif
    endfor
    symsize = 0.4
    
    tplot_options, 'version', 2
    tplot, plot_vars, position=poss, trange=time_range, noerase=1, nouttick=1

    ; manual xticknames.
    tpos = poss[*,-1]
    xstep = 5*60d
    xrange = time_range
    xtickv = make_bins(xrange, xstep, inner=1)
    xtickn = time_string(xtickv, tformat='hhmm')
    xtickn[0] = time_string(xtickv[0], tformat='YYYY MTH DD        ')
    
    plot, xrange, [0,1], $
        xstyle=5, ystyle=5, $
        nodata=1, noerase=1, position=tpos
    foreach msg, xtickn, id do begin
        tx = xtickv[id]
        tmp = convert_coord(tx,0, data=1, to_normal=1)
        tx = tmp[0]
        ty = tpos[1]-ychsz*1.1
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach


    label_yshift = -ychsz*0.6
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*5
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor


    ; Add FMLT.
    ewo_var = 'thg_asf_ewo'
    pid = where(plot_vars eq ewo_var, count)
    
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    model_index = where(all_models eq external_model)
    
    if count ne 0 then begin
        tpos = poss[*,pid]
        get_data, ewo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        
        dp_mlts = []
        foreach key, sort_uniq(event_info.rbsp.keys()), probe_id do begin
            info = (event_info.rbsp)[key]
            prefix = info['prefix']
            probe = info['probe']
            sc_color = info['sc_color']
            fmlt_var = prefix+'fmlt_'+internal_model+'_north'
            get_data, fmlt_var, times, fmlts
            fmlts = fmlts[*,model_index]+24
            oplot, times, fmlts, color=sc_color, linestyle=2
            ; Add sc.
            tx = time_range[0]
            ty = interpol(fmlts,times,tx)
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]+ychsz*0.3
            sc_name = info['sc_name']
            msg = strupcase(sc_name+'-'+probe)
            xyouts, tx,ty,normal=1, msg, charsize=label_size, alignment=0, color=sc_color
            ; Add dp time.
            tx = dp_times[probe_id]
            ty = interpol(fmlts,times,tx)
            plots, tx,ty, psym=8, symsize=0.5, color=sc_color
            dp_mlts = [dp_mlts,ty]
            
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            txs = tmp[0]+[0,0]
            tys = [tmp[1],poss[1,-1]]
            plots, txs,tys, normal=1, linestyle=1, color=sc_color
        endforeach
        
        ; overall slope.
        plots, dp_times,dp_mlts, linestyle=1
        dp_omega = total(dp_mlts*[-1,1])/total(dp_times*[-1,1])*60*15
        event_info['dp_omega'] = dp_omega
        tx = mean(dp_times)
        ty = mean(dp_mlts)
        tmp = convert_coord(tx, ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*2
        ty = tmp[1]+ychsz*0.2
        msg = tex2str('omega')+' = '+string(abs(dp_omega),format='(F3.1)')+' deg/min'
        xyouts, tx,ty,normal=1, alignment=0, msg, charsize=label_size
    endif


    return, poss
end



function fig_2015_0416_0800_rbsp_v02, event_info=event_info, test=test

;---Load data and settings.
    version = 'v02'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    event_info['dp_times'] = time_double('2015-04-16/'+['08:07:42','08:10:00'])

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_rbsp_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0


;---Figure out panel size.
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]
    uniform_ticklen = -abs_ychsz*0.15

    asi_times = event_info['dp_times']
    nasi_panel = n_elements(asi_times)
    top_xsize = 5.5
    top_margins = [10,4,12,1]+[1,0,1,0]*3
    top_ysize = 1.6

    margins = [0,0,0,0]
    ypad = 0
    bottom_ysize = 1.8

    all_poss = panel_pos(plot_file, ypans=[top_ysize,bottom_ysize], $
        margins=margins, panid=[0,0], pansize=[top_xsize,top_ysize], ypad=ypad, fig_size=fig_size)
    
    top_pos = all_poss[*,0]
    mid_pos = sgcalcpos(1, $
        margins=top_margins, region=top_pos, $
        xchsz=xchsz, ychsz=ychsz)

    bottom_pos = all_poss[*,1]
    bottom_margins = [8,2,7,0]
    bottom_poss = sgcalcpos(1,2, xpad=3, margins=bottom_margins, region=bottom_pos)
    left_pos = bottom_poss[*,0]
    right_pos = bottom_poss[*,1]
    rbsp_poss = [[left_pos],[right_pos]]
 
    
;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, mid_pos
        panel_list.add, left_pos
        panel_list.add, right_pos
        for ii=0,nasi_panel-1 do panel_list.add, reform(asi_poss[*,ii])

        sgopen, 0, size=fig_size

        foreach tpos, panel_list do begin
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, $
                xtickformat='(A1)', ytickformat='(A1)', $
                xticklen=xticklen, yticklen=yticklen, $
                nodata=1, noerase=1, position=tpos
        endforeach
        stop
    endif
    
    

;--Panels.
    sgopen, plot_file, size=fig_size
    tpos = fig_2015_0416_0800_rbsp_v02_mid_panel(mid_pos, event_info=event_info)
    tpos = fig_2015_0416_0800_rbsp_v02_rbsp_panel(rbsp_poss, event_info=event_info)

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 0
print, fig_2015_0416_0800_rbsp_v02(event_info=event_info, test=test)
end