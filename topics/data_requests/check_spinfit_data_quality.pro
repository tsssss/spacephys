;+
; Load necessary data to check the data quality of spinfit data.
;-

pro check_spinfit_data_quality, time_range, probe=probe

    rbspx = 'rbsp'+probe
    prefix = 'rbsp'+probe+'_'

;---Load L3 data.
    year = time_string(day,tformat='YYYY')
    base = rbspx+'_efw-l3_'+time_string(day,tformat='YYYYMMDD')+'_v04.cdf'
    file = join_path([root_dir,rbspx,'l3',year,base])
    cdf2tplot, file
    ;rbsp_efw_phasef_read_density, time_range, probe=probe
    ;rbsp_efw_phasef_read_l3, time_range, probe=probe

;---Load E uvw.
    rbsp_efw_phasef_read_e_uvw, time_range, probe=probe

;---Load B MGSE to get the angle of <B,spin-axis>.
    rbsp_efw_phasef_read_b_mgse, time_range, probe=probe
    get_data, prefix+'b_mgse', times, b_mgse
    store_data, 'spinfit_angle', times, acos(b_mgse[*,0]/snorm(b_mgse))*constant('deg'), $
        limits={constant: 90+[-1,1]*15}

    tplot, [prefix+'e_uvw','efield_in_corotation_frame_spinfit_'+['mgse','edotb_mgse'],$
        'spacecraft_potential','density','lshell','spinfit_angle'], trange=time_range


end
