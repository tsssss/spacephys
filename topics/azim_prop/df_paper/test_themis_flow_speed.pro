;+
; Plot flow to equatorial plane.
;-

pro plot_themis_flow_along_orbit, time_range, probes=probes

    time_step = 600d
    common_times = make_bins(time_range, time_step)
    xrange = [10,-30]
    yrange = [1,-1]*20
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    coef = 1d/50    ; 1 Re/50 km/s.
    ntime = n_elements(common_times)
    ndim = 3

    xspan = abs(total(xrange*[-1,1]))
    yspan = abs(total(yrange*[-1,1]))

    plot_file = 0
    tpos = panel_pos(plot_file, pansize=[xspan,yspan]*5d/xspan, fig_size=fig_size)
    sgopen, plot_file, xsize=fig_size[0], ysize=fig_size[1]
    hsize = 2
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, $
        ystyle=1, yrange=yrange, ytitle=ytitle, $
        nodata=1, noerase=1, iso=1, position=tpos
    tmp = smkarthm(0,2*!dpi,50,'n')
    txs = cos(tmp)
    tys = sin(tmp)
    polyfill, txs>0, tys, color=sgcolor('white')
    polyfill, txs<0, tys, color=sgcolor('gray')
    oplot, txs, tys
    foreach rr, [5,10,15] do begin
        oplot, txs*rr,tys*rr, linestyle=1
    endforeach

    colors = [constant('rgb'),sgcolor('purple'),sgcolor('orange')]
    rs = list()
    us = list()
    root_dir = join_path([googledir(),'works','global_efield','data'])
    foreach probe, probes, probe_id do begin
        prefix = 'th'+probe+'_'
        color = colors[probe_id]

        ; Load V and R.
;        if check_if_update(prefix+'u_gsm', time_range) then themis_read_ion_vel, time_range, probe=probe
        if check_if_update(prefix+'r_gsm', time_range) then themis_read_orbit, time_range, probe=probe
        base = prefix+'ion_vel_'+time_string(time_range[0],tformat='YYYY')+'.tplot'
        file = join_path([root_dir,'ion_vel',base])
        tplot_restore, filename=file

        ; Convert data from gsm to sm.
        vars = prefix+['u','r']
        foreach var, vars do begin
            in_var = var+'_gsm'
            out_var = var+'_sm'
            vec = get_var_data(in_var, times=times, in=time_range+[0,time_step])
            vec = cotran(vec, times, 'gsm2sm')
            store_data, out_var, times, vec
            if var eq prefix+'u' then begin
                lim = {$
                    display_type: 'vector', $
                    unit: 'km/s', $
                    short_name: 'U!S!Uion!N!R', $
                    coord: 'SM', $
                    coord_labels: constant('xyz'), $
                    colors: constant('rgb')}
            endif else if var eq prefix+'r' then begin
                lim = {$
                    display_type: 'vector', $
                    unit: 'Re', $
                    short_name: 'R', $
                    coord: 'SM', $
                    coord_labels: constant('xyz'), $
                    colors: constant('rgb')}
            endif
            add_setting, out_var, /smart, lim
            
            
            vec_avg = fltarr(ntime,ndim)
            foreach time, common_times, time_id do begin
                time_index = where_pro(times, '[]', time+[0,time_step])
                for dim_id=0,ndim-1 do vec_avg[time_id,dim_id] = mean(vec[time_index,dim_id],/nan)
            endforeach
            store_data, out_var, common_times, vec_avg
        endforeach

        r_sm = get_var_data(prefix+'r_sm')
        u_sm = get_var_data(prefix+'u_sm')

        oplot, r_sm[*,0], r_sm[*,1], color=color
        foreach time, common_times, time_id do begin
            x0 = r_sm[time_id,0]
            y0 = r_sm[time_id,1]
            x1 = x0+u_sm[time_id,0]*coef
            y1 = y0+u_sm[time_id,1]*coef
            ;plots, [x0,x1], [y0,y1], data=1
            arrow, x0,y0, x1,y1, data=1, solid=0, hsize=hsize, color=color
        endforeach

        rs.add, snorm(r_sm[*,0:1]), /extract
        us.add, snorm(u_sm[*,0:1]), /extract
    endforeach

    rs = rs.toarray()
    us = us.toarray()

    stop

end

time_range = time_double(['2014-08-28/09:00','2014-08-28/11:00'])
;time_range = time_double(['2014-08-28/09:00','2014-08-28/24:00'])
time_range = time_double(['2013-06-07/04:30','2013-06-07/06:30'])
time_range = time_double(['2013-06-07/03:00','2013-06-07/07:30'])
;time_range = time_double(['2016-10-13/12:00','2016-10-13/13:00'])
probes = ['a','d','e']

time_range = time_double(['2008-01-09/11:00','2008-01-09/12:00'])
probes = ['a','c','d','e']

time_range = time_double(['2008-03-20/10:00','2008-03-20/14:00'])
probes = ['b','c','d','e']

plot_themis_flow_along_orbit, time_range, probes=probes
end
