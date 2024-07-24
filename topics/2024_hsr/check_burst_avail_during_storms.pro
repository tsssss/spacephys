
storm_list = join_path([googledir(),'codes','idl','spacephys','topics','pflux_grant','interesting_events','alfven_arc','data','storm_list_v1.txt'])
log_file = join_path([srootdir(),'check_burst_avail_during_storms.txt'])
if file_test(log_file) eq 0 then ftouch, log_file

lines = read_all_lines(storm_list)
nstorm = n_elements(lines)
storm_trs = dblarr(nstorm,2)
b1_info = orderedhash()


probes = ['a','b']
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_'
    var = prefix+'efw_vb1_time_rate'
    rbsp_efw_phasef_read_b1_time_rate, probe=probe
endforeach




for ii=0,nstorm-1 do begin
    info = strsplit(lines[ii],' ',extract=1)
    storm_tr = time_double(info)
    
    storm_trs[ii,*] = storm_tr
    print, lines[ii]
    the_info = dictionary()
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        var = prefix+'efw_vb1_time_rate'
        b1_trs = get_var_data(var, times=times)
        index = where(b1_trs[*,0] ge storm_tr[0] and b1_trs[*,1] le storm_tr[1], nb1_tr)
        if nb1_tr eq 0 then begin
            msg = 'No B1 data ...'
        endif else begin
            duration = total(b1_trs[index,1]-b1_trs[index,0])/3600
            if duration lt 0.5 then begin
                msg = 'B1 data too short ...'
            endif else begin
                msg = 'RBSP-'+probe+'has b1 for '+string(duration,format='(F10.2)')+' h'
                the_info[probe] = duration
            endelse
        endelse
        lprmsg, msg, log_file
    endforeach

    key = lines[ii]
    if n_elements(the_info) ne 0 then b1_info[key] = the_info
endfor


foreach key, b1_info.keys(), id do begin
    lprmsg, ''
    msg = string(id,format='(I3)')+'    '+key
    lprmsg, msg
    the_info = b1_info[key]
    foreach probe, the_info.keys() do begin
        duration = the_info[probe]
        msg = probe+'    '+string(duration,format='(F10.2)')
        lprmsg, msg
    endforeach
endforeach


;---Collect orbit info.
    count_var = 'stat_counts'
    if check_if_update(count_var) then begin
        mlat_bins = make_bins([50,70],1)
        mlt_bins = make_bins([-6,6],0.5)
        nmlat_bin = n_elements(mlat_bins)-1
        nmlt_bin = n_elements(mlt_bins)-1
        mlat_centers = (mlat_bins[1:nmlat_bin-1]+mlat_bins[0:nmlat_bin-2])*0.5
        mlt_centers = (mlt_bins[1:nmlt_bin-1]+mlt_bins[0:nmlt_bin-2])*0.5
        stat_counts = fltarr(nmlt_bin,nmlat_bin)
        foreach key, b1_info.keys(), id do begin
            tr = time_double(strsplit(key,' ',extract=1))
            foreach probe, probes do begin
                ;r_gsm_var = rbsp_read_orbit(tr, probe=probe, get_name=1)
                ;del_data, r_gsm_var
                r_gsm_var = rbsp_read_orbit(tr, probe=probe)
                mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
                mlt_var = mlat_vars.mlt
                mlts = get_var_data(mlt_var, in=tr, times=times)

                f_gsm_var = lets_trace_to_ionosphere(orbit_var=r_gsm_var, $
                    internal_model='dipole', external_model='t89', hemisphere='north')
                f_vars = lets_read_mlat_vars(orbit_var=f_gsm_var, suffix='_map')
                fmlat_var = f_vars.mlat
                fmlats = get_var_data(fmlat_var, in=tr, times=times)

                for ii=0,nmlt_bin-1 do begin
                    for jj=0,nmlat_bin-1 do begin
                        mlt_range = mlt_bins[ii:ii+1]
                        mlat_range = mlat_bins[jj:jj+1]
                        index = where(fmlats ge mlat_range[0] and fmlats lt mlat_range[1] and $
                            mlts ge mlt_range[0] and mlts lt mlt_range[1], count)
                        stat_counts[ii,jj] += count
                    endfor
                endfor
            endforeach
        endforeach
        store_data, count_var, 0, stat_counts
    endif
    
test = 0
    margins = [6,4,8,1]
    plot_file = join_path([srootdir(),'2024_hsr_mlt_mlat_coverage_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    tpos = panel_pos(plot_file, fig_size=fig_size, pansize=[2,1]*1.5, margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    uniform_ticklen = -ychsz*fig_size[0]*0.15

    zrange = [0,50]
    zminor = 5
    ztitle = 'Duration (h)'
    top_color = 254
    ct = 49
    zzs = bytscl(stat_counts/60,min=zrange[0],max=zrange[1],top=top_color)
    
    
    cbpos = tpos
    cbpos[0] = cbpos[2]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.8
    zticklen = uniform_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
    sgcolorbar, findgen(top_color), $
        ztitle=ztitle, zrange=zrange, ct=ct, position=cbpos, zticklen=zticklen, zminor=zminor
    
    xrange = minmax(mlt_bins)
    xtitle = 'MLT (h)'
    
    yrange = minmax(mlat_bins)
    ytitle = 'MLat (deg)'
    
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    
    sgtv, zzs, position=tpos, ct=ct, resize=1
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, $
        ystyle=1, yrange=yrange, ytitle=ytitle, $
        position=tpos, nodata=1, noerase=1, $
        xticklen=xticklen, yticklen=yticklen
        
    
    if keyword_set(test) then stop
    sgclose

end