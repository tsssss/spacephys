;+
;-

probes = 'th'+letters('e')
time_range = time_double(['2007-03-01','2009-10-01'])
plot_base = 20

test = 1

;probes = 'rbsp'+letters('b')
;time_range = time_double(['2012-10-01','2017-10-01'])
;plot_base = 10

;probes = 'g'+['13','14','15']
;time_range = time_double(['2012-10-01','2017-10-01'])
;plot_base = 1

ndim = 3

foreach probe, probes, probe_id do begin
    dfs = azim_df_read_df(time_range, probe=probe, /no_scaled_height)
    ndf = n_elements(dfs)
    r_sms = fltarr(ndf,ndim)
    df_widths = fltarr(ndf)
    df_heights = fltarr(ndf)
    foreach df, dfs, df_id do begin
        r_sms[df_id,*] = df.obs_r_sm
        df_widths[df_id] = df.width
        df_heights[df_id] = df.height
    endforeach
    df_mlts = pseudo_mlt(r_sms)
    df_rxys = snorm(r_sms[*,0:1])

    plot_file = join_path([project.plot_dir,'df_distr','df_distr_'+probe+'.pdf'])
    if keyword_set(test) then plot_file = probe_id+plot_base


;---Filter by height.
    df_min_height = 1.
    index = where(df_heights ge df_min_height)
    df_mlts = df_mlts[index]
    df_rxys = df_rxys[index]
    df_heights = df_heights[index]
    df_widths = df_widths[index]

;---Bin data.
    mlt_range = [-1,1]*9
    mlt_bin_size = 0.5
    mlt_bin_boundarys = make_bins(mlt_range, mlt_bin_size)
    nmlt_bin = n_elements(mlt_bin_boundarys)-1
    mlt_bin_centers = mlt_bin_boundarys[0:nmlt_bin-1]+mlt_bin_size*0.5

    rxy_range = [4.,30]
    rxy_bin_size = 1
    rxy_bin_boundarys = make_bins(rxy_range, rxy_bin_size,/inner)
    nrxy_bin = n_elements(rxy_bin_boundarys)-1
    rxy_bin_centers = rxy_bin_boundarys[0:nrxy_bin-1]+rxy_bin_size*0.5

    height_range = [0,max(df_heights)]
    width_range = [0,max(df_widths)]
    str_delta_theta = tex2str('Delta')+tex2str('theta')
    str_delta_t = tex2str('Delta')+'T'

    bin_map_index = ptrarr(nmlt_bin,nrxy_bin)
    the_mlts = df_mlts
    for ii=0, nmlt_bin-1 do begin
        the_index = lazy_where(the_mlts, '[)', mlt_bin_boundarys[ii:ii+1], count=count)
        the_rxys = df_rxys[the_index]
        if count eq 0 then continue
        for jj=0, nrxy_bin-1 do begin
            index = lazy_where(the_rxys, '[)', rxy_bin_boundarys[jj:jj+1], count=count)
            if count eq 0 then continue
            bin_map_index[ii,jj] = ptr_new(the_index[index])
        endfor
    endfor

    ; Counts.
    count_ct = 49
    binned_count = fltarr(nmlt_bin,nrxy_bin)
    for ii=0, nmlt_bin-1 do begin
        for jj=0, nrxy_bin-1 do begin
            if ~ptr_valid(bin_map_index[ii,jj]) then continue
            binned_count[ii,jj] = n_elements(*bin_map_index[ii,jj])
        endfor
    endfor


;---Figure settings.
    xrange = [10.,-35]
    xstep = 10.
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtitle = 'SM X (Re)'

    yrange = [-1,1]*22.5
    ystep = 10.
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5

    panel_aspect_ratio = abs(total(xrange*[-1,1])/total(yrange*[-1,1]))
    panel_ysize = 2.
    panel_xsize = panel_ysize*panel_aspect_ratio
    nxpanel = 3
    nypanel = 1
    xpad = 2
    ypad = 0
    margins = [8,5,2,5]
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    fig_xsize = panel_xsize*nxpanel+xpad*(nxpanel-1)*abs_xchsz+total(margins[[0,2]])*abs_xchsz
    fig_ysize = panel_ysize*nypanel+ypad*(nypanel-1)*abs_ychsz+total(margins[[1,3]])*abs_ychsz
    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    poss = sgcalcpos(nypanel, nxpanel, position=pos, $
        ypad=ypad, xpad=xpad, xchsz=xchsz, ychsz=ychsz)


    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.


    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify, /inch
    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    lineskip = constant('lineskip')

    height_color = sgcolor('silver')
    point_psym = 3  ; dot.
    point_psym = 1  ; cross.
    point_symsize = 0.2
    width_color = sgcolor('wheat')


    fit_color = sgcolor('black')
    fit_psym = 1
    fit_symsize = 0.3
    fit_linestyle = 2

    stddev_color = sgcolor('red')

    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    label_size = constant('label_size')


    ; Stat settings.
    stat_ratio = 0.2
    stat_min_count = 15
    data_mag = 5

    top_color = 255
    bottom_color = 15
    data_colors = smkarthm(bottom_color, top_color, 1, 'dx')

    ; Circle for earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)

    ; Add magnetopause.
    magn_test = fltarr(nangle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    pdyn = 10.
    tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=pdyn)
    magn_pos = magn_pos[*,[0,2,1]]
    magn_pos = [magn_pos,magn_pos]
    magn_pos[nangle:nangle*2-1,1] *= -1

    fig_labels = letters(nxpanel)+'.'
    for panel_id=0,nxpanel-1 do begin
        case panel_id of
            0: begin
                data = []
                data_range = [stat_min_count,200]
                ztitle = '# of dipolarizations'
                end
            1: begin
                data = df_heights
                data_range = [0,20]
                ztitle = str_delta_theta+' (deg)'
                end
            2: begin
                data = df_widths/60
                data_range = [10,50]
                ztitle = str_delta_t+' (min)'
                end
        endcase

        tpos = poss[*,panel_id,0]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])


        if n_elements(data) eq 0 then binned_data = binned_count else begin
            binned_data = fltarr(nmlt_bin,nrxy_bin)
            for ii=0, nmlt_bin-1 do begin
                for jj=0, nrxy_bin-1 do begin
                    the_count = binned_count[ii,jj]
                    if the_count lt stat_min_count then continue
                    the_data = data[*bin_map_index[ii,jj]]
                    the_data = the_data[sort(the_data)]
                    the_data = the_data[the_count*stat_ratio:the_count*(1-stat_ratio)]
;                    binned_data[ii,jj] = stddev(the_data)*data_mag
                    binned_data[ii,jj] = mean(the_data)*data_mag
                endfor
            endfor
        endelse
        zzs = bytscl(binned_data, min=data_range[0], max=data_range[1], top=top_color-bottom_color)+bottom_color



        ytitle = 'SM Y (Re)'
        ytickformat = ''
        if panel_id ne 0 then begin
            ytitle = ''
            ytickformat = '(A1)'
        endif else begin
            ytitle = 'SM Y (Re)'
            ytickformat = ''
        endelse

        ; Create coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            /nodata, /noerase, /isotropic, position=tpos
        ;polyfill, xrange[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=sgcolor('wheat')
        ; Draw each bin.
        rrs = rxy_bin_boundarys
        tts = (180+mlt_bin_boundarys)*15*constant('rad')
        for ii=0, nmlt_bin-1 do begin
            for jj=0, nrxy_bin-1 do begin
                if binned_count[ii,jj] lt stat_min_count then continue
                the_tt = tts[ii:ii+1]
                the_rr = rrs[jj:jj+1]
                txxs = the_rr[[0,1,1,0,0]]*cos(the_tt[[0,0,1,1,0]])
                tyys = the_rr[[0,1,1,0,0]]*sin(the_tt[[0,0,1,1,0]])
                the_color = sgcolor(zzs[ii,jj], ct=count_ct)
                polyfill, txxs, tyys, color=the_color
            endfor
        endfor

        ; Add earth.
        polyfill, circle_x<0, circle_y, color=sgcolor('silver')
        plots, circle_x, circle_y

        ; Add magnetopause.
        oplot, magn_pos[*,0], magn_pos[*,1]

        ; Add more circles.
        foreach tmp, [5,10,20] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
        foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp

        ; Add mlt_range.
        foreach tmp, mlt_range do begin
            tt = tmp*15*constant('rad')+!dpi
            oplot, rxy_range*cos(tt), rxy_range*sin(tt)
        endforeach

        ; Draw axis.
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
            /nodata, /noerase, /isotropic, position=tpos

        ; Add colorbar.
        cbpos = tpos[[0,3,2,3]]+ychsz*[0,lineskip,0,lineskip+half_ychsz]
        sgcolorbar, data_colors, /horizontal, position=cbpos, $
            zrange=data_range, ztitle=ztitle, ct=count_ct

        ; Add figure label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*(full_ychsz+lineskip)
        xyouts, tx,ty,/normal, fig_labels[panel_id]

    endfor

    if keyword_set(test) then stop
    sgclose

endforeach
end
