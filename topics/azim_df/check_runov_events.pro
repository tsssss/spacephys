;+
; Use the timing in Table 1 to do 3-sc timing whereever available.
; 
; Give up, no point to prove Runov is correct or wrong.
;-

infos = list()
infos.add, [$
    'thb','2009-02-27/07:51:25.60','-20.1,-0.6,-1.5', $
    'thc','2009-02-27/07:52:34.40','-16.7,-1.6,-2.2', $
    'thd','2009-02-27/07:54:06.60','-11.1,-2.7,-2.1', $
    'the','2009-02-27/07:54:10.30','-11.1,-1.8,-2.4']
infos.add, [$
    'thc','2009-03-05/03:13:04.10','-17.9, 1.4,-1.6', $
    'thd','2009-03-05/03:14:26.77','-10.3, 1.5,-1.7', $
    'tha','2009-03-05/03:14:33.75',' -9.1, 2.4,-2.3', $
    'the','2009-03-05/03:14:38.50',' -9.2, 2.4,-1.5']
infos.add, [$
    'thc','2009-03-09/09:06:49.00','-14.3,-0.8,-1.2', $
    'the','2009-03-09/09:07:45.10','-11.4,-1.2,-1.6', $
    'thd','2009-03-09/09:07:58.70','-11.1,-2.1,-1.3']
infos.add, [$
    'thc','2009-03-15/08:49:05.08','-13.7, 0.1,-0.9', $
    'thb','2009-03-15/08:49:11.35','-12.6, 0.1,-0.2', $
    'tha','2009-03-15/08:49:19.30','-11.5,-0.2,-2.3', $
    'the','2009-03-15/08:49:28.00','-11.5,-0.2,-1.3', $
    'thd','2009-03-15/08:49:38.05','-11.3,-1.1,-1.0']
infos.add, [$
    'thc','2009-03-19/08:25:11.50','-13.4, 0.7,-0.6', $
    'thb','2009-03-19/08:25:19.20','-12.3, 0.7,-0.0', $
    'the','2009-03-19/08:25:30.10','-11.5, 0.6,-1.1', $
    'thd','2009-03-19/08:25:33.10','-11.4,-0.3,-0.9']
infos.add, [$
    'thd','2009-03-31/08:25:47.50','-11.2, 1.2,-0.1', $
    'thc','2009-03-31/08:25:50.00','-11.1, 1.5, 0.0', $
    'the','2009-03-31/08:25:58.25','-11.3, 2.2,-0.4', $
    'thb','2009-03-31/08:26:09.17',' -9.6, 1.2, 0.4']
ncol = 3
all_probes = letters('e')
nall_probe = n_elements(all_probes)

df_group_list = list()
foreach info, infos, info_id do begin
    nprobe = n_elements(info)/ncol
    info = reform(info, [ncol,nprobe])

;---Basic info.
    probes = strarr(nprobe)
    obs_times = dblarr(nprobe)
    r_gsm = fltarr(nprobe,3)
    
    for ii=0,nprobe-1 do begin
        probes[ii] = strmid(info[0,ii],2,1)
        obs_times[ii] = time_double(info[1,ii],tformat='YYYY-MM-DD/hh:mm:ss.ff')
        r_gsm = float(strsplit(info[2,ii],',',/extract))
    endfor
    r_sm = cotran(r_gsm, obs_times, 'gsm2sm')
    
    
;---Check all probes.
    sorted_probes = list([probes,all_probes],/extract)
    foreach probe, all_probes do begin
        index = sorted_probes.where(probe, count=count)
        if count eq 1 then continue
        for ii=1,count-1 do sorted_probes.remove, index[ii]
    endforeach
    sorted_probes = sorted_probes.toarray()
    
    time_range = minmax(obs_times)+[-0.5,1]*300
    foreach probe, sorted_probes do begin
        themis_read_orbit, time_range, probe=probe
        
    endforeach

        
    
    
    sgopen, 0, xsize=8, ysize=8, xchsz=xchsz, ychsz=ychsz
    plot_time_range = time_range
    xrange = plot_time_range
    for ii=1,10 do begin
        xtickv = make_bins(xrange,60*ii, /inner)
        if n_elements(xtickv) le 4 then break
    endfor
    xticks = n_elements(xtickv)-1
    xtickn = time_string(xtickv, tformat='hh:mm')
    xtickn[0] += '!C'+time_string(xtickv[0],tformat='YYYY-MTH-DD')
    xminor = total(xtickv[0:1]*[-1,1])/60
    if xminor le 1 then xminor = 6
    
    pos_xrange = [1,-40]
    pos_yrange = [-1,1]*15
    
    
    foreach probe, letters('e') do begin
        themis_read_bfield, time_range, probe=probe
        prefix = 'th'+probe+'_'
        b_gsm = get_var_data(prefix+'b_gsm', times=times)
        b_sm = cotran(b_gsm, times, 'gsm2sm')
        tilt = azim_df_calc_tilt(b_sm)
        store_data, prefix+'tilt_sm', times, tilt
        add_setting, prefix+'tilt_sm', /smart, dictionary($
            'display_type', 'scalar', $
            'unit', 'deg', $
            'short_name', tex2str('theta') )
        store_data, prefix+'bz_gsm', times, b_gsm[*,2]
        add_setting, prefix+'bz_gsm', /smart, dictionary($
            'display_type', 'scalar', $
            'unit', 'nT', $
            'short_name', 'Bz' )
    endforeach
    
    xticklen_chsz = -0.30
    yticklen_chsz = -0.40
    
    colors = sgcolor(['purple','red','green','cyan','blue'])
    probe_colors = dictionary()
    foreach probe, all_probes, ii do probe_colors[probe] = colors[ii]
    
    
    poss = sgcalcpos(nall_probe, position=[0.1,0.35,0.45,0.95])
    foreach probe, sorted_probes, panel_id do begin
        tpos = poss[*,panel_id]
        prefix = 'th'+probe+'_'
        
        xtickformat = (panel_id eq nall_probe-1)? '': '(A1)'
        yys = get_var_data(prefix+'bz_gsm', times=xxs, in=xrange)
        yrange = sg_autolim(minmax(yys))
        ytitle = '(nT)'
        
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=1, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
            position=tpos
        oplot, xxs, yys
        tx = tpos[2]+xchsz*1
        ty = (tpos[1]+tpos[3])*0.5-ychsz*0.3
        label = 'TH'+strupcase(probe)
        xyouts, tx,ty,/normal, label, color=probe_colors[probe]
        
        index = where(probes eq probe, count)
        if count ne 0 then begin
            tx = obs_times[index[0]]
            oplot, tx+[0,0],yrange, linestyle=1
        endif
    endforeach
    
    
    tpos = sgcalcpos(position=[0.1,0.1,0.45,0.25])
    plot, pos_xrange, pos_yrange, /nodata, /noerase, $
        xstyle=1, xrange=pos_xrange, xtitle='GSM X (Re)', $
        ystyle=1, yrange=pos_yrange, ytitle='GSM Y (Re)', $
        position=tpos, /isotrop
    foreach probe, sorted_probes, panel_id do begin
        prefix = 'th'+probe+'_'

        index = where(probes eq probe, count)
        if count ne 0 then begin
            the_time = obs_times[index[0]]
        endif else begin
            the_time = mean(time_range)
        endelse
        r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
        tx = r_gsm[0]
        ty = r_gsm[1]
        plots, tx,ty, psym=1, color=probe_colors[probe]
    endforeach
    
    
    poss = sgcalcpos(nall_probe, position=[0.6,0.35,0.95,0.95])
    foreach probe, sorted_probes, panel_id do begin
        tpos = poss[*,panel_id]
        prefix = 'th'+probe+'_'

        xtickformat = (panel_id eq nall_probe-1)? '': '(A1)'
        yys = get_var_data(prefix+'tilt_sm', times=xxs, in=xrange)
        yrange = sg_autolim(minmax(yys))
        ytitle = '(deg)'

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=1, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
            position=tpos
        oplot, xxs, yys
        tx = tpos[2]+xchsz*1
        ty = (tpos[1]+tpos[3])*0.5-ychsz*0.3
        label = 'TH'+strupcase(probe)
        xyouts, tx,ty,/normal, label, color=probe_colors[probe]

        index = where(probes eq probe, count)
        if count ne 0 then begin
            tx = obs_times[index[0]]
            oplot, tx+[0,0],yrange, linestyle=1
        endif
    endforeach
    
    
    tpos = sgcalcpos(position=[0.6,0.1,0.95,0.25])
    plot, pos_xrange, pos_yrange, /nodata, /noerase, $
        xstyle=1, xrange=pos_xrange, xtitle='SM X (Re)', $
        ystyle=1, yrange=pos_yrange, ytitle='SM Y (Re)', $
        position=tpos, /isotrop
    foreach probe, sorted_probes, panel_id do begin
        prefix = 'th'+probe+'_'
        
        index = where(probes eq probe, count)
        if count ne 0 then begin
            the_time = obs_times[index[0]]
        endif else begin
            the_time = mean(time_range)
        endelse
        r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
        r_sm = cotran(r_gsm, the_time, 'gsm2sm')
        tx = r_sm[0]
        ty = r_sm[1]
        plots, tx,ty, psym=1, color=probe_colors[probe]
    endforeach
;        
;
;        xtickformat = (panel_id eq nall_probe-1)? '': '(A1)'
;        yys = get_var_data(prefix+'tilt_sm', times=xxs, in=xrange)
;        yrange = sg_autolim(minmax(yys))
;        ytitle = '(deg)'
;
;        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
;        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
;
;        plot, xrange, yrange, /nodata, /noerase, $
;            xstyle=1, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
;            ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
;            position=tpos
;        oplot, xxs, yys
;        tx = tpos[2]+xchsz*1
;        ty = (tpos[1]+tpos[3])*0.5-ychsz*0.3
;        label = 'TH'+strupcase(probe)
;        xyouts, tx,ty,/normal, label, color=probe_colors[probe]
;
;        index = where(probes eq probe, count)
;        if count ne 0 then begin
;            tx = obs_times[index[0]]
;            oplot, tx+[0,0],yrange, linestyle=1
;        endif
;    endforeach
    
    
    
    

    
    df_group = dictionary('id', info_id+1)
    df_group['time_range'] = time_range
    df_group['region'] = 'around_midn'
    df_group['search_name'] = 'runov'
    df_group['probes'] = probes
    df_group['df_list'] = list()
    df_group['df_list'] = list()
    df_group['triad_list'] = list()
    df_group['edge_list'] = list()

    df_list = df_group.df_list
    for ii=0,nprobe-1 do begin
        
    endfor
    stop
endforeach

end
