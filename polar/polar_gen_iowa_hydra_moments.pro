;+
; Save iowa hydra moments to a single cdf file.
; The original iowa moments are in electron and ion, and do not contain eflux.
;
; Save to <sdata>/polar/hydra/yyyy/po_hyd_moments_yyyy_mmdd_v01.cdf
;-

pro polar_gen_iowa_hydra_moments, date, file=cdffn


    secofday = 86400d
    ev = 1d/11600d
    min_gyro = 2d       ; bad data for gyrotropy > 2.
    fillval = !values.d_nan

    if size(date,/type) eq 7 then ut0 = time_double(date) else ut0 = date
    ut0 = ut0-(ut0 mod secofday)
    utr = ut0+[0d,86400]

    loc_root = join_path([sdiskdir('Research'),'sdata'])
    local_path = join_path([loc_root,'opt_hydra','moment_data'])
    version = (n_elements(version) eq 0)? 'v[0-9.]*': version

    datdir = sdiskdir('Research')+'/sdata/polar/hydra/'+time_string(utr[0],tformat='YYYY')
    if file_test(datdir,/directory) eq 0 then file_mkdir, datdir
    cdffn = datdir+'/po_hyd_moments_'+time_string(utr[0],tformat='YYYY_MMDD')+'_v01.cdf'

    if file_test(cdffn) eq 1 then file_delete, cdffn
    ginfo = {$
        title: 'Polar Hydra moments, merging from the original Hydra moment from Iowa data'}
    if file_test(cdffn) eq 0 then scdfwrite, cdffn, gattribute=ginfo
    store_data, '*', /delete


    species = ['ion','ele']
    foreach tspecie, species do begin
        ion = (tspecie eq 'ion')? 1: 0
        pre0 = (tspecie eq 'ion')? 'ion_': 'ele_'

        base_pattern = keyword_set(ion)? '%Y%m%d_hyd_momi_'+version+'.cdf': '%Y%m%d_hyd_mom_'+version+'.cdf'
        base_names = apply_time_to_pattern(base_pattern, utr)
        files = local_path+path_sep()+base_names
        file = file_search(files[0])
        if file[0] eq '' then begin
            file_delete, cdffn
            return
        endif

        ; the variables to be loaded from CDF.
        vars = ['time','epoch','sc_potential','n','t','temp','v','q','status','b_status', $
            'tparl','tperp','gyrotropy','skew','angle_tparl_bave','mfe_b_avg','b_avg','mfe_v_avg_sigma','b_proxy']
        dat = scdfread(file[0], vars, skt=skt)
        if size(dat,/type) ne 8 then begin
            file_delete, cdffn
            return
        endif
        date_str = time_string(utr[0], tformat='YYYYMMDD')

        time = *dat[where(vars eq 'time')].value
        rec_start = 0
        rec_count = n_elements(time)
        utname = 'ut_sec'

        ; Define the data structure to return
        data_temp = { hydra_moments_data, $
            moments_status:long(0), $
            time:0.0, $
            time_ssm:0.0, $
            date:'????????', $
            mean_energy:0.0, $
            density:0.0, $
            u_pay:fltarr(4),$
            u_gsm:fltarr(4),$
            u_gse:fltarr(4),$
            u_perp_pay:fltarr(4),$
            u_perp_gsm:fltarr(4),$
            u_parl:0.0,$
            ktparl:0.0,$
            ktperp:0.0,$
            pres_parl:0.0,$
            pres_perp:0.0,$
            gyrotropy:0.0, $
            skew:0.0, $
            anisotropy:0.0, $
            b_status:long(0), $
            b_pay:fltarr(4), $
            b_pay_sig:fltarr(4), $
            b_pay_angle_cone:0.0, $
            b_pay_status:0, $
            b_gsm:fltarr(4), $
            b_gse:fltarr(4), $
            b_proxy:fltarr(4), $
            v_alf:0.0,$
            m_tot:0.0,$
            m_parl:0.0,$
            m_perp:0.0,$
            sc_potential:0.0, $
            q_pay:fltarr(4), $
            q_gsm:fltarr(4), $
            q_gse:fltarr(4), $
            q_perp_pay:fltarr(4), $
            q_perp_gsm:fltarr(4), $
            q_parl:0.0, $
            qparl_over_bmag:0.0, $
            diagonalization_status:0, $
            angle_bavg_bproxy:0.0, $
            angle_bavg_bproxy_sig:0.0, $
            cosangle_q_bavg:0.0, $
            cosangle_q_bavg_sig:0.0, $
            cosangle_u_bavg:0.0, $
            cosangle_u_bavg_sig:0.0, $
            qnu_parl:0.0, $ ; Parl current density contribution in microAmps/meter^2
            eflux_parl:0.0, $ ; ergs/cm^2/s
            eflux_parl_over_b:0.0 $ ; ergs/cm^2/s/nT
            }

        data_values = replicate( data_temp, rec_count )

        sc_potential = *dat[where(vars eq 'sc_potential')].value

        density = *dat[where(vars eq 'n')].value
        t = *dat[where(vars eq 't')].value  ; scalar temperature, K.
        temp_tensor = *dat[where(vars eq 'temp')].value ; temperature tensor, K.
        u = *dat[where(vars eq 'v')].value
        q = *dat[where(vars eq 'q')].value
        moments_status = *dat[where(vars eq 'status')].value
        b_status = *dat[where(vars eq 'b_status')].value

        ; Assume that all of the proper diagonalization is taken
        ; place when the moments are generated.  The quantities below
        ; are sensitive to the diagonalization process...
        tparl = *dat[where(vars eq 'tparl')].value
        tperp = *dat[where(vars eq 'tperp')].value
        gyrotropy = *dat[where(vars eq 'gyrotropy')].value
        skew = *dat[where(vars eq 'skew')].value
        angle_bavg_bproxy = *dat[where(vars eq 'angle_tparl_bave')].value
        b_avg = *dat[where(vars eq 'mfe_b_avg')].value
        b_hk = *dat[where(vars eq 'b_avg')].value
        b_avg_sig = *dat[where(vars eq 'mfe_b_avg_sigma')].value
        b_proxy = *dat[where(vars eq 'b_proxy')].value

        b_used= replicate(0, rec_count)
        idx = where(b_used eq -1, cnt)
        if cnt ne 0 then b_avg[idx,*]= 0
        idx = where(b_used eq 1, cnt)
        if cnt ne 0 then b_avg[idx,*] = b_hk[idx,*]
        message, 'Using "mfe_b_avg" from the moments file for all parl/perp determinations...', /continue


        ; Not available in the moments file---------------------
        angle_bavg_bproxy_sig = fltarr(rec_count)
        diagonalization_status = intarr(rec_count) + 1
        b_avg_angle_cone = fltarr(rec_count)
        b_avg_status = intarr(rec_count)

        diagonalization_status = reform(diagonalization_status)
        tparl = reform(tparl)
        tperp = reform(tperp)
        gyrotropy = reform(gyrotropy)
        skew = reform(skew)
        angle_bavg_bproxy = reform(angle_bavg_bproxy)
        angle_bavg_bproxy_sig = reform(angle_bavg_bproxy_sig)
        b_avg = reform(b_avg)
        b_avg_sig = reform(b_avg_sig)
        b_avg_angle_cone = reform(b_avg_angle_cone)
        b_avg_status = reform(b_avg_status)
        b_proxy = reform(b_proxy)

        ; ********************** DERIVED PARAMETERS *******************************

        moments_species = 1 ; for ion.
        mass= ( [ 9.1e-28,1.67e-24 ] )( moments_species ) ; g
        charge= ( [-4.8d-10, 4.8d-10] )( moments_species ) ; esu

        u = transpose(u)
        q = transpose(q)
        temp_tensor = transpose(temp_tensor)
        b_avg = transpose(b_avg)
        b_avg_sig = transpose(b_avg_sig)
        b_proxy = transpose(b_proxy)

        setenv, 'HYDRA_DDCAL_DATA_PATH='+join_path([loc_root,'opt_hydra','l1_data'])
        defsysv, '!hyd_exception', hyd_get_exception()

        ; Calculate the derived parameters(dp), load the data_values structure
        dp_hydra_moments, date_str, time, sc_potential, density, t, u, q, $
            temp_tensor, b_avg, b_avg_sig, b_avg_status, b_avg_angle_cone, $
            diagonalization_status, tparl, tperp, gyrotropy, skew, $
            angle_bavg_bproxy, angle_bavg_bproxy_sig, b_proxy, $
            data_values, charge, mass

        data_values.b_status = reform(b_status)
        data_values.moments_status = reform(moments_status)


    ;---Save to CDF.
        if keyword_set(ion) then begin
            tval =  time+ut0
            ainfo = {$
                FIELDNAM:'UT for the moments',$
                UNITS:'s (UT)',$
                VAR_TYPE:'support_data'}
            scdfwrite, cdffn, utname, value=tval, cdftype='CDF_DOUBLE', attribute=ainfo

            ; spacecraft potential.
            vname = 'sc_potential'
            tval = sc_potential
            ainfo = {$
                FIELDNAM: 'SC potential', $
                UNITS: 'V', $
                VAR_TYPE: 'data', $
                DEPEND_0: utname}
            scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        endif

        ; gyrotropy. bad data for gyrotropy > 2.
        vname = pre0+'gyrotropy'
        tval = gyrotropy
        ainfo = {$
            FIELDNAM: 'Gyrotropy for '+tspecie, $
            UNITS: '#', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; density.
        vname = pre0+'density'
        tval = data_values.density
        sinfo = {$
            FIELDNAM:'Density '+tspecie, $
            UNITS:'cm!U-3!N', $
            VAR_TYPE:'data', $
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; average temperature.
        vname = pre0+'t_avg'
        tval = t*ev
        sinfo = {$
            FIELDNAM:'Temperature '+tspecie, $
            UNITS:'eV', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; parallel temperature.
        vname = pre0+'tparl'
        tval = tparl*ev
        sinfo = {$
            FIELDNAM:'Parallel temperature '+tspecie, $
            UNITS:'eV', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; perpendicular temperature.
        vname = pre0+'tperp'
        tval = tperp*ev
        sinfo = {$
            FIELDNAM:'Perpendicular temperature '+tspecie, $
            UNITS:'eV', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; number flux.
        vname = pre0+'number_flux'
        coef = keyword_set(ion)? 1d/1.6e-10: -1d/1.6e-10
        tval = data_values.qnu_parl*coef    ; in #/cm^2-s.
        ainfo = {$
            FIELDNAM: 'Number flux for '+tspecie, $
            UNITS: '#/(cm!U2!N-s)', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; energy flux.
        vname = pre0+'energy_flux'
        tval = data_values.eflux_parl
        idx = where(gyrotropy ge min_gyro, cnt)
        if cnt ne 0 then tval[idx] = fillval
        ainfo = {$
            FIELDNAM: 'Energy flux for '+tspecie, $
            UNITS: 'mW/m!U2!N', $
            VAR_TYPE: 'data', $
            FILLVAL: 'NaN', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        ; enthalpy.

        ; bulk velocity.
        vname = pre0+'bulk_velocity'
        tval = data_values.u_parl
        ainfo = {$
            FIELDNAM: 'Parallel bulk velocity for '+tspecie, $
            UNITS: 'km/s', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo

        vname = pre0+'bulk_velocity_gsm'
        tval = transpose(data_values.u_gsm[0:2,*])
        ainfo = {$
            FIELDNAM: 'Bulk velocity in GSM for '+tspecie, $
            UNITS: 'km/s', $
            VAR_TYPE: 'data', $
            DEPEND_0: utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
    endforeach

end

;utr0 = time_double(['1997-05-20','2008-12-31'])
;utr0 = time_double(['1998-03-22','2008-12-31'])
;utr0 = time_double(['2004-01-01','2004-12-31'])
utr0 = time_double(['2001-01-01','2003-12-31'])

utr0 = time_double(['1997-04-29','2004-12-31'])

uts = break_down_times(utr0)
nut = n_elements(uts)
for i=0, nut-1 do begin
    date = time_string(uts[i],tformat='YYYY-MM-DD')
    print, 'Processing '+date+' ...'
    polar_gen_iowa_hydra_moments, date, file=cdffn
    print, '    Saved to '+cdffn+' ...'
endfor

end