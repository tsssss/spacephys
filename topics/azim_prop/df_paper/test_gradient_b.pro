;+
; Calculates the magnetic gradient of a dipole field at certain location.
;-
;

;---Settings.
    test = 0

    re = 6378d  ; km.
    time = time_double('2014-08-28/10:00')
    ps = geopack_recalc(time)

    tilt = 0.           ; makes gsm coincide with sm, thus only bz is non-zero.
    dx = 0.1    ; in Re.
    rad = !dpi/180d
    deg = 180d/!dpi
    
    foreach dis, [5,10] do begin
        rgsm = [-dis,0,0]
        geopack_dip, rgsm[0],rgsm[1],rgsm[2], bx,by,bz, tilt=tilt
        b0gsm = [bx,by,bz]
        geopack_dip, rgsm[0]+dx*0.5,rgsm[1],rgsm[2], bx,by,bz, tilt=tilt
        b1gsm = [bx,by,bz]
        geopack_dip, rgsm[0]-dx*0.5,rgsm[1],rgsm[2], bx,by,bz, tilt=tilt
        b2gsm = [bx,by,bz]

        grad_b = (b2gsm[2]-b1gsm[2])
        print, 'Distance (Re)', dis
        print, 'Grad B, dipole, nT/Re', grad_b
    endforeach
;    stop
    
    
;---Load data for 2014-08-28 event.
    event_id = '2014_0828_10'
    project = azim_prop_load_project()
    event_info = project.events[event_id]
    time_range = event_info.time_range
    tplot_restore, file=event_info.file
    
    probes = event_info.sorted_probes
    omega_azim = event_info.omega_azim/60*rad  ; from deg/min to rad/sec.
    foreach probe, probes do begin
        get_data, probe+'_db_mag', times, dbmag
        ref_dis = snorm(event_info[probe].ref_rsm)  ; in Re.
        v_azim = omega_azim*ref_dis
        gradb = deriv(times,dbmag)/v_azim
        time_lag = event_info[probe].time_lag
        store_data, probe+'_gradb', times, smooth(gradb,10), limits={$
            ytitle:'(nT/Re)', labels:strupcase(project[probe].short_name)}
    endforeach
    
    if keyword_set(test) then begin
        file = 0
    endif else begin
        file = join_path([shomedir(),'test_gradb_'+event_id+'.pdf'])
    endelse
    sgopen, file, xsize=5, ysize=6
    tplot, probes+'_gradb', trange=time_range, /novtitle
    sgclose
    
    
    if keyword_set(test) then begin
        file = 0
    endif else begin
        file = join_path([shomedir(),'test_gradb_dbmag_'+event_id+'.pdf'])
    endelse
    sgopen, file, xsize=5, ysize=6
    tplot, probes+'_db_mag', trange=time_range, /novtitle
    sgclose
    
    
end