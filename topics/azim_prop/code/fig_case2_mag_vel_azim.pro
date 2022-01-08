;+
; Plot the mag data and calculate the azimuthal velocity.
;-

pro fig_case2_mag_vel_azim

;---Load basic data.
    _2014_0828_10_load_data, event_info=event_info

    test = 1

;---Settings.
    mag_site_mlat_range = [65d,70]
    mag_site_mlon_range = [-170,10]
    aurora_mlon_range = [-100d,-55]    ; exclude the region where aurora extends, and use it as a separator.
    label_time_range = time_double(['2014-08-28/10:00','2014-08-28/10:05'])
    time_range_plot = time_double(['2014-08-28/09:50','2014-08-28/11:00'])


;---Select mag sites.
    mag_var = 'thg_dbh'
    get_data, mag_var, times, gmag_data
    ntime = n_elements(times)
    gmag_sites = event_info.gmag.sites      ; all gmag sites that are loaded.
    gmag_info = event_info.gmag_info        ; info for all gmag sites.
    used_sites = list(gmag_sites,/extract)  ; a subset of gmag_sites.



;---Ad-hoc times.
;    fac_time_infos = hash($
;        'tik' , '2014-08-28/10:29:30', $
;        'pbk' , '2014-08-28/10:20:50', $
;        'kian', '2014-08-28/10:12:00', $
;        'fsmi', '2014-08-28/10:12:50', $
;        'gill', '2014-08-28/10:18:10', $
;        'fcc' , '2014-08-28/10:20:00')
;    foreach tinfo, mag_infos, ii do begin
;        if ~fac_time_infos.haskey(tinfo.name) then continue
;        mag_infos[ii].time = time_double(fac_time_infos[tinfo.name])
;    endforeach


    ; Filter with the overall mlat/mlon range.
    if n_elements(mag_site_mlat_range) eq 2 then begin
        foreach site, used_sites, ii do begin
            if gmag_info[site].mlat lt min(mag_site_mlat_range) then used_sites.remove, ii
            if gmag_info[site].mlat gt max(mag_site_mlat_range) then used_sites.remove, ii
        endforeach
    endif

    if n_elements(mag_site_mlon_range) eq 2 then begin
        foreach site, used_sites, ii do begin
            if gmag_info[site].mlon lt min(mag_site_mlon_range) then used_sites.remove, ii
            if gmag_info[site].mlon gt max(mag_site_mlon_range) then used_sites.remove, ii
        endforeach
    endif


    ; Filter with aurora expansion.
    if n_elements(aurora_mlon_range) eq 0 then message, 'No auroral MLon range'

    used_gmag = dictionary()
    used_gmag['east'] = dictionary('site',list())
    foreach site, used_sites do if gmag_info[site].mlon ge max(aurora_mlon_range) then used_gmag.east.site.add, site
    used_gmag['west'] = dictionary('site',list())
    foreach site, used_sites do if gmag_info[site].mlon le min(aurora_mlon_range) then used_gmag.west.site.add, site


    ; Save the used gmag data.
    parts = used_gmag.keys()
    gmag_unit = 'nT'
    foreach part, parts do begin
        sites = used_gmag[part].site
        foreach site, sites do begin
            tvar = site+'_mag'
            index = where(gmag_sites eq site)
            store_data, tvar, times, gmag_data[*,index]
            add_setting, tvar, /smart, {$
                display_type: 'scalar', $
                unit: gmag_unit, $
                short_name: strupcase(site)}
        endforeach
    endforeach


;---Combine data.
    base_value_scale = 550d     ; nT.
    base_value_shift = -400d    ; nT. The basic shift when combining mag data.
    color_start = 100
    color_end = 250
    used_gmag['east'].color_table = 60
    used_gmag['west'].color_table = 55
    foreach part, parts do begin
        ; Get the sites, sort by mlon.
        sites = used_gmag[part].site
        nsite = n_elements(sites)
        mlons = fltarr(nsite)
        foreach site, sites, ii do mlons[ii] = gmag_info[site].mlon
        sites = sites[sort(mlons)]

        ; Collect data and shift to avoid overlapping.
        gmag_data = fltarr(ntime, nsite)
        foreach site, sites, ii do begin
            gmag_data[*,ii] = get_var_data(site+'_mag')
            gmag_data[*,ii] -= gmag_data[0,ii]
        endforeach

        base_values = fltarr(nsite)
        base_values[0] = 0
        for ii=1, nsite-1 do begin
            base_values[ii] = base_value_shift*total(minmax(gmag_data[*,ii-1])*[-1,1])/base_value_scale + base_values[ii-1]
            gmag_data[*,ii] += base_values[ii]
        endfor
        used_gmag[part].base_value = base_values

        ; Save the combined var.
        used_gmag[part].var = part+'_gmag_combo'
        colors = round(smkarthm(color_start,color_end, nsite, 'n'))
        foreach color, colors, ii do colors[ii] = sgcolor(colors[ii], ct=used_gmag[part].color_table)
        store_data, used_gmag[part].var, times, gmag_data
        add_setting, used_gmag[part].var, /smart, {$
            display_type:'stack', $
            labels:strupcase(sites.toarray()), $
            colors:colors}
    endforeach
    stop

;---Plot the data.
    color['timebar'] = sgcolor('red')

    figure_size = [4,5]
    if keyword_set(test) then begin
        file = 0
        magnify = 2
    endif else begin
        file = join_path([sparentdir(srootdir()),'plot','fig_case2_mag_vel_azim.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=figure_size[0], ysize=figure_size[1], magnify=magnify

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    half_ychsz = 0.4

    poss = sgcalcpos(3, lmargin=8, tmargin=1, rmargin=3, bmargin=3)
    pos1 = poss[*,0]
    pos2 = poss[*,1]
    pos3 = poss[*,2]
    panel_yshift = ychsz*1
    pos1[1] += panel_yshift
    pos2[3] += panel_yshift
    pos2[1] -= panel_yshift
    pos3[3] -= panel_yshift

    ;---Prepare colors.
        east_color_table = 60
        west_color_table = 55
        color_start = 100
        color_end = 250

        nmag_var = n_elements(east_infos)
        east_colors = round(smkarthm(color_start, color_end, nmag_var, 'n'))
        for ii=0, nmag_var-1 do east_colors[ii] = sgcolor(east_colors[ii], ct=east_color_table)

        nmag_var = n_elements(west_infos)
        west_colors = round(smkarthm(color_end, color_start, nmag_var, 'n'))
        for ii=0, nmag_var-1 do west_colors[ii] = sgcolor(west_colors[ii], ct=west_color_table)


    ;---Combine the settings to fit a for loop.
        west_plot_infos = hash($
            'position', pos1, $
            'var_name', 'west_mag', $
            'colors', west_colors, $
            'site_info', west_infos)
        east_plot_infos = hash($
            'position', pos3, $
            'var_name', 'east_mag', $
            'colors', east_colors, $
            'site_info', east_infos)
        plot_infos = list(east_plot_infos,west_plot_infos)


    ;---Plot.
        ticklen = -0.01
        xrange = time_range_plot    ; 1 hour.
        xminor = 4          ; 5 min.
        xtickv = smkarthm(xrange[0], xrange[1], xminor*5*60, 'dx')
        xticks = n_elements(xtickv)-1
        xticklen = ticklen*2
        xtickn = time_string(xtickv,tformat='hh:mm')
        xtickn[0] = 'UT (hh:mm)     '
        yrange = reverse(mag_site_mlon_range)
        yminor = 4
        ytickv = make_bins(yrange, yminor*10, /inner)
        yticks = n_elements(ytickv)
        yticklen = ticklen
        ytitle = 'MLon (deg)'
        psym = 6
        symsize = 0.7
        label_size = 0.7
        arrow_size = (size(file,/type) eq 7)? 80:8
        fig_label_x = 2
        earth_rotation = 360d/86400     ; deg/sec.


        foreach plot_info, plot_infos, pp do begin
            tinfos = plot_info['site_info']
            tmag_var = plot_info['var_name']
            tpos = plot_info['position']
            tcolors = plot_info['colors']
            tlabel = (strsplit(tmag_var,'_',/extract))[0]

            mag_vars = tinfos.name+'_mag'
            nmag_var = n_elements(mag_vars)

            ; Sort by MLon.
            index = sort(tinfos.mlon)
            mag_vars = mag_vars[index]
            tinfos = tinfos[index]

            ; Combine data from each site into one tplot variable.
            gmag_data = fltarr(ntime, nmag_var)
            for ii=0, nmag_var-1 do begin
                get_data, mag_vars[ii], times, data
                gmag_data[*,ii] = data
            endfor

            ; Collect info for plotting.
            for ii=0, nmag_var-1 do begin
                if ii eq 0 then begin
                    base_value = data[0]
                endif else begin
                    base_value = base_value_shift*total(minmax(gmag_data[*,ii-1])*[-1,1])/base_value_scale + mag_data_infos[ii-1].base_value
                endelse
                stop
                mag_data_infos[ii].base_value = base_value
                gmag_data[*,ii] += base_value
            endfor

            store_data, tmag_var, times, gmag_data, limits={$
                colors: tcolors, labels:mag_data_infos.name, ytitle:'', ytickformat:'(A1)'}
            mag_yrange = sg_autolim(gmag_data)
            plot, xrange, mag_yrange, $
                xstyle=5, $
                ystyle=5, $
                /noerase, /nodata, position=tpos
            for ii=0, nmag_var-1 do begin
                oplot, times, gmag_data[*,ii], color=tcolors[ii]
                tx = mag_data_infos[ii].time
                ty = interpol(gmag_data[*,ii], times, tx)
                tmp = convert_coord(tx,ty, /data, /to_normal)
                tx = tmp[0]
                ty = tmp[1]-ychsz*0.2
                arrow, tx,ty-ychsz*0.5, tx,ty, /normal, /solid, hsize=arrow_size
            endfor

            xtickformat = (tlabel eq 'east')? '': '(A1)'
            ytickformat = '(A1)'
            plot, xrange, mag_yrange, $
                xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, xtickn=xtickn, $
                ystyle=1, yrange=mag_yrange, yticks=2, yminor=5, yticklen=yticklen, ytickformat=ytickformat, $
                /noerase, /nodata, position=tpos

            tx = tpos[0]-xchsz*1
            for ii=0, nmag_var-1 do begin
                ty = mag_data_infos[ii].base_value
                ty = mean(gmag_data[lazy_where(times, label_time_range),ii])
                tmp = convert_coord(xrange[0],ty, /data, /to_normal)
                ty = tmp[1]-ychsz*half_ychsz
                xyouts, tx,ty, /normal, alignment=1, mag_data_infos[ii].name, color=tcolors[ii]
            endfor

            tx = tpos[0]+xchsz*1
            tmp = convert_coord(xrange[0], yrange[0], /data, /to_normal)
            ty = tmp[1]
            tmp = convert_coord(xrange[0], yrange[0]+base_value_shift, /data, /to_normal)
            dy = abs(tmp[1]-ty)
            ty = tpos[1]+ychsz*1+dy*0.5
            plots, tx+[0,0], ty+[-1,1]*dy*0.5, /normal
            plots, tx+[-1,1]*xchsz*0.5, ty-dy*0.5, /normal
            plots, tx+[-1,1]*xchsz*0.5, ty+dy*0.5, /normal
            xyouts, tx+xchsz*0.5, ty-ychsz*half_ychsz*label_size, /normal, sgnum2str(abs(base_value_shift))+' nT', charsize=label_size
            tx = fig_label_x*xchsz
            ty = tpos[3]-ychsz*half_ychsz
            fig_label = (tlabel eq 'east')? 'c.': 'a.'
            xyouts, tx,ty, /normal, fig_label

            ; Add to the middle panel.
            xtickformat = '(A1)'
            tpos = pos2
            plot, xrange, yrange, position=tpos, $
                xstyle=5, xrange=xrange, $
                ystyle=5, yrange=yrange, $
                /nodata, /noerase

            txs = tinfos.time
            tys = tinfos.mlon
            for ii=0, nmag_var-1 do begin
                plots, txs[ii], tys[ii], psym=psym, color=tcolors[ii], symsize=symsize
            endfor

            ; Linear fit.
            fit_coef = linfit(txs, tys)
            linfit_yrange = yrange
            linfit_xrange = (linfit_yrange-fit_coef[0])/fit_coef[1]
            linfit_xrange = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
            linfit_yrange = linfit_xrange*fit_coef[1]+fit_coef[0]
            oplot, linfit_xrange, linfit_yrange, linestyle=1
            tx = tpos[2]-xchsz*15
            ty = (tlabel eq 'east')? tpos[1]+ychsz*1: tpos[3]-ychsz*2
            xyouts, tx,ty, /normal, 'V!D'+tlabel+' !N: '+sgnum2str(abs(fit_coef[1]-earth_rotation), nsgn=2)+' deg/sec', charsize=label_size

            plot, xrange, yrange, position=tpos, $
                xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
                ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
                /nodata, /noerase
            tx = fig_label_x*xchsz
            ty = tpos[3]-ychsz*half_ychsz
            xyouts, tx,ty, /normal, 'b.'

            tx = tpos[0]+xchsz*1
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty, /normal, 'Magnetometers in!C['+sgnum2str(mag_site_mlat_range[0])+','+sgnum2str(mag_site_mlat_range[1])+'] deg MLat', charsize=label_size

        endforeach


    if keyword_set(test) then stop
    sgclose

end
