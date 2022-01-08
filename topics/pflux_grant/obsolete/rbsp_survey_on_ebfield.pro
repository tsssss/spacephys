;+
; Survey the whole mission for data quality of the hires E field and 4-sec B field.
;
; The purpose of the plot is to show:
;   1. How well data around perigee are covered by R = 4Re.
;   2. How well the model data work, if good then it's much easier to process the data around perigee.
;   3. Plot Vsc to see if it can be used to remove bad E data.
;-

pro rbsp_survey_on_ebfield, project=project, time_range=time_range

test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    probe_infos = project.probe_infos
    probes = probe_infos.keys()
    probes = 'rbspa'
    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_ebfield_survey'])
    efield_resolution = 'hires'
    efield_time_step = 1d/16
    efield_nrec = secofday/efield_time_step
    efield_min_nrec = 0.5*efield_nrec
    boom_lengths = [100d,100,12]
    spin_period = 12d
    dc_offset_nspin = 1
    max_valid_e = 300.
    max_valid_v = 100.

    bfield_resolution = '4sec'
    bfield_time_step = 4.
    bfield_nrec = secofday/bfield_time_step
    bfield_min_nrec = 0.5*bfield_nrec

    orbit_time_step = 60.
    omega = (2*!dpi)/86400d  ;Earth's rotation angular frequency
    perigee_lshell = 4.

    model = 't89'
    par = 2.
    ndim = 3
    uvw = constant('uvw')
    uv = uvw[0:1]
    xyz = constant('xyz')
    rgb = constant('rgb')
    small_range = [-1,1]*10.    ; mV/m.
    large_range = small_range*30

    fig_xsize = 7.
    fig_ysize = 5.
    ypans = [1.,2]
    nypanel = n_elements(ypans)
    margins = [8,3,5,2]
    panel_ypad = 0.5
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, /inch
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=panel_ypad, margins=margins, xchsz=xchsz, ychsz=ychsz)
    sgclose, /wdelete

    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    label_xshift = 6

    log_file = join_path([project.data_dir,'rbsp_survey_on_ebfield.log'])
    if file_test(log_file) eq 0 then ftouch, log_file

    foreach probe, probes do begin
        the_probe = strmid(probe,0,1,/reverse)
        prefix = 'rbsp'+the_probe+'_'

        the_time_range = (n_elements(time_range) eq 2)? time_range: rbsp_info('efw_l2_data_range', probe=the_probe)
;the_time_range >= time_double('2013-12-07')
        ndate = total(the_time_range*[-1,1])/secofday
        dates = smkarthm(the_time_range[0],secofday,ndate,'x0')

        foreach date, dates do begin
            date_time_range = date+[0,secofday]
            efield_times = make_bins(date_time_range, efield_time_step)
            efield_ntime = n_elements(efield_times)
            bfield_times = make_bins(date_time_range, bfield_time_step)
            bfield_ntime = n_elements(bfield_times)
            orbit_times = make_bins(date_time_range, orbit_time_step)
            norbit_time = n_elements(orbit_times)

            year = time_string(date,tformat='YYYY')
            base_name = probe+'_efield_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
            plot_file = join_path([plot_dir,probe,year,base_name])
            if file_test(plot_file) eq 1 then continue
            del_data, '*'

            ; Read B data.
            b_var = prefix+'b_gsm'
            del_data, b_var
            rbsp_read_bfield, date_time_range, probe=the_probe, resolution=bfield_resolution
            get_data, b_var, times, b_gsm
            if n_elements(b_gsm[*,0]) le bfield_min_nrec then begin
                msg = strupcase(probe)+ ': not enough good B data on '+time_string(date,tformat='YYYY/MM-DD')
                errmsg = handle_error(msg)
                lprmsg, msg, log_file
                continue
            endif

            ; Read E data.
            vsvy_var = prefix+'vsvy'
            if date ge time_double('2016-02-28') then begin
                rbsp_read_efw, date_time_range, id='l2%vsvy-highres2', probe=the_probe
            endif else begin
                rbsp_read_efw, date_time_range, id='l2%vsvy-highres', probe=the_probe
            endelse
            rename_var, 'vsvy', to=vsvy_var
            get_data, vsvy_var, times, vsvy
            if n_elements(vsvy) eq 1 then begin
                msg = strupcase(probe)+ ': no V data on '+time_string(date,tformat='YYYY/MM-DD')
                errmsg = handle_error(msg)
                lprmsg, msg, log_file
                continue
            endif
            vsvy = vsvy[*,0:3]
            index = where(abs(vsvy) ge max_valid_v, count)
            if count ne 0 then vsvy[index] = !values.f_nan
            store_data, vsvy_var, times, vsvy

            data_quality = 1
            if n_elements(vsvy[*,0]) le efield_min_nrec then data_quality = 0
            index = where(finite(snorm(vsvy)), count)
            if count le efield_min_nrec then data_quality = 0
            if data_quality eq 0 then begin
                msg = strupcase(probe)+ ': not enough good E data on '+time_string(date,tformat='YYYY/MM-DD')
                errmsg = handle_error(msg)
                lprmsg, msg, log_file
                continue
            endif


            ; Read orbit and quaternion data.
            r_var = prefix+'r_gsm'
            rbsp_read_orbit, date_time_range, probe=the_probe
            q_var = prefix+'q_uvw2gsm'
            rbsp_read_quaternion, date_time_range, probe=the_probe


            ; Calc Vsc.
            vsvy = get_var_data(vsvy_var, at=efield_times)
            vsc_var = prefix+'vsc'
            vsc = [$
                [mean(vsvy[*,0:1], dimension=2)], $
                [mean(vsvy[*,2:3], dimension=2)]]
            store_data, vsc_var, efield_times, vsc
            add_setting, vsc_var, /smart, {$
                display_type: 'vector', $
                unit: 'V', $
                short_name: 'V!S!USC!N!R', $
                coord: 'UVW', $
                coord_labels: uv, $
                colors: rgb[0:1] }


            ; Calc E[uv] and |E|.
            e_uv_var = prefix+'e_uv'
            eu = (vsvy[*,0]-vsvy[*,1])/boom_lengths[0]*1e3   ; V -> V/m -> mV/m.
            ev = (vsvy[*,2]-vsvy[*,3])/boom_lengths[1]*1e3
            ew = fltarr(efield_ntime)
            index = where(abs(eu) ge max_valid_e, count)
            if count ne 0 then eu[index] = !values.f_nan
            index = where(abs(ev) ge max_valid_e, count)
            if count ne 0 then ev[index] = !values.f_nan

            ; Remove dc-offset.
            width = dc_offset_nspin*spin_period/efield_time_step
            eu = eu-smooth(eu, width, /edge_truncate, /nan)
            ev = ev-smooth(ev, width, /edge_truncate, /nan)
            e_uv = [[eu],[ev]]
            store_data, e_uv_var, efield_times, e_uv
            add_setting, e_uv_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'E', $
                coord: 'UVW', $
                coord_labels: uv, $
                colors: rgb[0:1] }

            emag = snorm(e_uv)
            emag_var = prefix+'emag'
            store_data, emag_var, efield_times, emag
            add_setting, emag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'mV/m', $
                short_name: '|E|'}



            ; Calculate B model and dB.
            r_gsm = get_var_data(r_var, at=bfield_times)
            b0_gsm = fltarr(bfield_ntime,ndim)
            for ii=0, bfield_ntime-1 do begin
                tilt = geopack_recalc(bfield_times[ii])
                ; in-situ position
                rx = r_gsm[ii,0]
                ry = r_gsm[ii,1]
                rz = r_gsm[ii,2]
                ; in-situ B field.
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
                geopack_t89, par, rx,ry,rz, dbx,dby,dbz
                b0_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
            endfor
            b0_var = prefix+'b0_gsm'
            store_data, b0_var, bfield_times, b0_gsm
            add_setting, b0_var, /smart, {$
                display_type: 'vector', $
                unit: 'nT', $
                short_name: 'B0!S!UT89!N!R', $
                coord: 'GSM', $
                coord_labels: xyz}

            b_gsm = get_var_data(b_var, at=bfield_times)
            db_gsm = b_gsm-b0_gsm
            db_var = prefix+'db_gsm'
            store_data, db_var, bfield_times, db_gsm
            add_setting, db_var, /smart, {$
                display_type: 'vector', $
                unit: 'nT', $
                short_name: 'dB', $
                coord: 'GSM', $
                coord_labels: xyz}

            dbmag = snorm(db_gsm)
            dbmag_var = prefix+'dbmag'
            store_data, dbmag_var, bfield_times, dbmag
            add_setting, dbmag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'nT', $
                short_name: '|dB|'}


            ; Calculate E model and dE.
            r_gsm = get_var_data(r_var, at=orbit_times)
            v_gsm = fltarr(norbit_time,ndim)
            for ii=0, ndim-1 do v_gsm[*,ii] = deriv(r_gsm[*,ii])
            v_gsm *= (constant('re')/orbit_time_step)
            v_var = prefix+'v_gsm'
            store_data, v_var, orbit_times, v_gsm
            add_setting, v_var, /smart, {$
                display_type: 'vector', $
                unit: 'km/s', $
                short_name: 'V', $
                coord: 'GSM', $
                coord_labels: xyz, $
                colors: rgb}
            b_gsm = get_var_data(b_var, at=orbit_times)
            evxb_gsm = scross(v_gsm,b_gsm)*1e-3   ; convert to mV/m.
            evxb_var = prefix+'evxb_gsm'
            store_data, evxb_var, orbit_times, evxb_gsm
            add_setting, evxb_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'VxB E', $
                coord: 'GSM', $
                coord_labels: xyz, $
                colors: rgb}

            r_gei = cotran(r_gsm, orbit_times, 'gsm2gei')
            vcoro_gei = fltarr(norbit_time,ndim)
            vcoro_gei[*,0] = -r_gei[*,1]*omega
            vcoro_gei[*,1] =  r_gei[*,0]*omega
            vcoro_gei[*,2] = 0.0
            vcoro_gsm = cotran(vcoro_gei, orbit_times, 'gei2gsm')
            ecoro_gsm = scross(vcoro_gsm, b_gsm)
            ecoro_var = prefix+'ecoro_gsm'
            store_data, ecoro_var, orbit_times, ecoro_gsm
            add_setting, ecoro_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'Coro E', $
                coord: 'GSM', $
                coord_labels: xyz, $
                colors: rgb}

            emod_var = prefix+'e0_gsm'
            emod_gsm = evxb_gsm+ecoro_gsm
            store_data, emod_var, orbit_times, emod_gsm
            add_setting, emod_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'Model E', $
                coord: 'GSM', $
                coord_labels: xyz, $
                colors: rgb}
            interp_time, emod_var, efield_times

            e0_uv_var = prefix+'e0_uv'
            the_var = prefix+'e0'
            rbsp_gsm2uvw, the_var+'_gsm', the_var+'_uvw', quaternion=q_var, probe=the_probe
            e0_uv = get_var_data(the_var+'_uvw', at=efield_times)
            e0_uv = e0_uv[*,0:1]
            store_data, e0_uv_var, efield_times, e0_uv
            add_setting, e0_uv_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'E0', $
                coord: 'UVW', $
                coord_labels: uv, $
                colors: rgb[0:1] }

            e0mag_var = prefix+'e0mag'
            e0mag = snorm(e0_uv)
            e0mag = smooth(e0mag, spin_period/efield_time_step)
            store_data, e0mag_var, efield_times, e0mag
            add_setting, e0mag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'mV/m', $
                short_name: '|E0|'}

            de_uv_var = prefix+'de_uv'
            de_uv = e_uv-e0_uv
            store_data, de_uv_var, efield_times, de_uv
            add_setting, de_uv_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'dE', $
                coord: 'UVW', $
                coord_labels: uv, $
                colors: rgb[0:1] }

            demag_var = prefix+'demag'
            demag = snorm(de_uv)
            store_data, demag_var, efield_times, demag
            add_setting, demag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'mV/m', $
                short_name: '|dE|'}


            dis = snorm(r_gsm)
            index = where(dis le perigee_lshell)
            time_ranges = time_to_range(orbit_times[index], time_step=orbit_time_step)


        ;---Load flags for bad data.
            rbsp_read_eclipse_flag, date_time_range, probe=the_probe
            rbsp_read_sdt_flag, date_time_range, probe=the_probe

            ;tplot, 'rbspa_'+['sdt_flag','eclipse_flag','de_uv','dbmag','demag']
            ;timebar, time_ranges, color=sgcolor('red')


            if keyword_set(test) then plot_file = test
            rbsp_survey_on_ebfield_gen_plot, date_time_range, filename=plot_file, probe=the_probe
            if keyword_set(test) then stop
        endforeach
    endforeach

end

rbsp_survey_on_ebfield, project=project
;time_range = time_double(['2013-06-07','2013-06-08'])
;probe = 'a'
;rbsp_survey_on_ebfield_gen_plot, time_range, filename=plot_file, project=project, probe=probe
end
