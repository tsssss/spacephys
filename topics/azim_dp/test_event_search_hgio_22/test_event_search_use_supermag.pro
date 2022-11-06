;+
; Use supermag data to find APD during storms.
;-

test = 0

plot_dir = join_path([googledir(),'works','azim_dp','plots','supermag_adp'])
time_range = ['2008','2019']
if keyword_set(test) then time_range = ['2008','2009']
if n_elements(storm_times) eq 0 then storm_times = ts_load_storm_list(time_range)
nstorm = n_elements(storm_times[*,0])
for storm_id=0,nstorm-1 do begin
    storm_time_range = reform(storm_times[storm_id,*])
    print, time_string(storm_time_range)
    
;    if keyword_set(test) then begin
;        sgopen
;        tplot, 'dst', trange=storm_time_range
;        stop
;    endif

    supermag_read_sme, storm_time_range
   
    get_data, 'ae', times, ae
    sme = get_var_data('sm_sme', at=times)
    store_data, 'ae', times, [[ae],[sme]]
    add_setting, 'ae', smart=1, dictionary($
        'display_type', 'stack', $
        'unit', 'nT', $
        'labels', ['AE','SME'], $
        'colors', sgcolor(['black','red']))
   
    storm_duration = total(storm_time_range*[-1,1])/3600d
    xsize = storm_duration*0.4    ; to inch.
    ysize = 4
    base = 'supermag_adp_'+time_string(storm_time_range[0],tformat='YYYY_MMDD')+'_v01.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsiz=xsize, ysize=ysize
    vars = ['dst','ae','sm_sme2d']
    nvar = n_elements(vars)
    ypans = [1,1,1.5]
    margins = [12,4,12,2]
    poss = sgcalcpos(nvar, margins=margins, ypans=ypans)
    options, vars, 'yticklen', -0.005
    tplot, vars, trange=storm_time_range, position=poss
    if keyword_set(test) then stop
    sgclose
endfor


end