;+
; Use cross correlation to get the time lag. Doesn't work well because correlation falls off too slow.
;-


;---Get data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1
    mlon_range = einfo.mag_mlon_range
    mlat_range = einfo.mag_mlat_range
    
    data_rate_asi = 3d      ; sec.
    data_rate_mag = 0.5d    ; sec.
    wdith = data_rate_asi/data_rate_mag*10
    
    
;---Use 'cigo','pokr','cmo' as reference point.
    mag_var = 'thg_dbh'
    get_data, mag_var, times, magh, sites
    ntime = n_elements(times)
    sites = strlowcase(sites)
    index = where(sites eq 'cigo')
    f1 = magh[*,index[0]]
    df1 = deriv(f1)
    sdespike, times, df1
    df1 = smooth(abs(df1), width)
    store_data, 'f1', times, df1
    ;for ii=1,ntime-1 do f1[ii] += f1[ii-1]
    
    index = where(sites eq 'kian')
    f2 = magh[*,index[0]]
    df2 = deriv(f2)
    sdespike, times, df2
    df2 = smooth(abs(df2), width)
    store_data, 'f2', times, df2
stop
    ;for ii=1,ntime-1 do f2[ii] += f2[ii-1]
    test_time_range = [-1,1]*600   ; sec.
    test_rec_range = test_time_range/data_rate_mag
    ntest_rec = total(test_rec_range*[-1,1])
    correlations = fltarr(ntest_rec)
    for ii=0, ntest_rec-1 do $
        correlations[ii] = c_correlate(f1,f2,ii+test_rec_range[0])
    max_correlation = max(correlations, index)
    print, 'Max Corr:', max_correlation
    print, 'dT:', (index+test_rec_range[0])*data_rate_mag
    stop
    
    index = where(sites eq 'cmo')
    f3 = magh[*,index[0]]
    
    print, c_correlate(f1,f3,0)
    
    
    stop

    mag_var = 'thg_dbh'
    get_data, mag_var, times, magh, sites, limit=lims
    site_infos = themis_read_mag_metadata(sites=sites)
    mag_mlons = site_infos.mlon
    mag_mlats = site_infos.mlat
    
    ct = 40
    colors = fltarr(nsite)
    for ii=0, nsite-1 do colors[ii] = sgcolor(lims.colors[ii],ct=ct)
    
    ntime = n_elements(times)
    nsite = n_elements(sites)
    width = data_rate_asi/data_rate_mag
    dmagh = fltarr(ntime,nsite)
    for ii=0, nsite-1 do begin
        dmagh[*,ii] = deriv(magh[*,ii])/data_rate_mag
        dmagh[*,ii] = smooth(dmagh[*,ii], width)
    endfor
    
    
    dmag_var = 'thg_dmagh'
    store_data, dmag_var, times, dmagh, limits=lims
    
    times = make_bins(time_range, data_rate_asi)
    ntime = n_elements(times)
    dmag_stddev = fltarr(ntime, nsite)
    get_data, dmag_var, uts, dmagh
    for ii=0, ntime-1 do begin
        index = lazy_where(uts, 'in', times[ii]+[0,data_rate_asi], count=count)
        if count eq 0 then continue
        for jj=0, nsite-1 do dmag_stddev[ii,jj] = stddev(dmagh[index,jj])
    endfor
    dmag_stddev_var = 'thg_dmagh_stddev'
    store_data, dmag_stddev_var, times, dmag_stddev, limits=lims
    
    mlat_filter = [65,67.5]
    index = where(mag_mlats le mlat_filter[1] and mag_mlats ge mlat_filter[0], complement=index2)
    colors[index2] = 0
    options, dmag_stddev_var, 'colors', colors
    
    
    var = mag_var
    db_shift = 200
    get_data, var, times, data
    data = data<db_shift
    for ii=0, nsite-1 do data[*,ii] -= db_shift*ii
    dmag_var2 = 'thg_dmagh_tmp'
    store_data, dmag_var2, times, data, limits=lims
    options, dmag_var2, 'colors', colors
    
end
