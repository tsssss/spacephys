;+
; Load E and B field and check Pflux and E/B ratio.
;-


test = 1
plot_dir = join_path([shomedir(),'check_themis_df_events'])


;---Events.
    event_list = list()


    ; Streamer event.
    event_list.add, dictionary($
        'time_range', time_double(['2008-01-21/07:00','2008-01-21/08:00']), $
        'onset_time_range', time_double(['2008-01-21/07:40','2008-01-21/07:50']), $
        'probes', ['a'] )
    event_list.add, dictionary($
        'time_range', time_double(['2008-01-21/07:00','2008-01-21/08:00']), $
        'onset_time_range', time_double(['2008-01-21/07:10','2008-01-21/07:20']), $
        'probes', ['a'] )
    event_list.add, dictionary($
        'time_range', time_double(['2008-01-21/07:00','2008-01-21/08:00']), $
        'onset_time_range', time_double(['2008-01-21/07:22','2008-01-21/07:28']), $
        'probes', ['a'] )
    event_list.add, dictionary($
        'time_range', time_double(['2008-01-21/07:00','2008-01-21/08:00']), $
        'onset_time_range', time_double(['2008-01-21/07:28','2008-01-21/07:38']), $
        'probes', ['a'] )

;    ; Angeloupolos 2008.
;    event_list.add, dictionary($
;        'time_range', ['2007-03-23/11:17','2007-03-23/11:22'], $
;        'probes', ['a','b','c','d','e'] )

;    ; Ogasawara+ 2011.
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-02-26/03:50','2008-02-26/04:20']), $
;        'onset_time_range', time_double(['2008-02-26/04:00','2008-02-26/04:10']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2009-04-14/08:50','2009-04-14/09:20']), $
;        'onset_time_range', time_double(['2009-04-14/09:07','2009-04-14/09:12']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($    ; Duplicated events?
;        'time_range', time_double(['2009-03-23/05:50','2009-03-23/06:20']), $
;        'onset_time_range', time_double(['2009-03-23/06:03','2009-03-23/06:10']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2009-03-14/04:50','2009-03-14/05:20']), $
;        'onset_time_range', time_double(['2009-03-14/05:05','2009-03-14/05:12']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2009-03-14/03:50','2009-03-14/04:20']), $
;        'onset_time_range', time_double(['2009-03-14/04:03','2009-03-14/04:10']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2009-03-14/01:20','2009-03-14/01:50']), $
;        'onset_time_range', time_double(['2009-03-14/01:35','2009-03-14/01:40']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2009-02-28/05:20','2009-02-28/05:50']), $
;        'onset_time_range', time_double(['2009-02-28/05:34','2009-02-28/05:44']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($    ; E field bad on e.
;        'time_range', time_double(['2009-02-22/08:10','2009-02-22/08:45']), $
;        'onset_time_range', time_double(['2009-02-22/08:25','2009-02-22/08:35']), $
;        'probes', ['d'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-03-30/06:15','2008-03-30/06:45']), $
;        'onset_time_range', time_double(['2008-03-30/06:28','2008-03-30/06:33']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-03-27/07:20','2008-03-27/07:45']), $
;        'onset_time_range', time_double(['2008-03-27/07:28','2008-03-27/07:32']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-03-11/05:40','2008-03-11/06:10']), $
;        'onset_time_range', time_double(['2008-03-11/05:53','2008-03-11/06:00']), $
;        'probes', ['d','e'] )
;;    event_list.add, dictionary($    ; E field bad.
;;        'time_range', time_double(['2008-03-03/07:50','2008-03-03/08:30']), $
;;        'onset_time_range', time_double(['2008-03-03/08:00','2008-03-03/08:10']), $
;;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-03-02/06:40','2008-03-02/07:20']), $
;        'onset_time_range', time_double(['2008-03-02/06:58','2008-03-02/07:05']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-02-27/02:10','2008-02-27/03:30']), $
;        'onset_time_range', time_double(['2008-02-27/02:42','2008-02-27/02:55']), $
;        'probes', ['d','e'] )
;    event_list.add, dictionary($    ; E field bad on e
;        'time_range', time_double(['2008-02-05/11:50','2008-02-05/12:30']), $
;        'onset_time_range', time_double(['2008-02-05/12:07','2008-02-05/12:13']), $
;        'probes', ['d'] )
;    event_list.add, dictionary($
;        'time_range', time_double(['2008-01-26/10:00','2008-01-26/10:50']), $
;        'onset_time_range', time_double(['2008-01-26/10:20','2008-01-26/10:30']), $
;        'probes', ['d','e'] )



;---Settings.
    tplot_options, 'labflag', -1
    rgb = constant('rgb')
    xyz = constant('xyz')

;---Loop through events.
    foreach event, event_list do begin
        del_data, '*'
        probes = event.probes
        time_range = event.time_range
        onset_time_range = event.onset_time_range

    ;---Load basic data: E, B, R, density.
        foreach probe, probes do begin
            prefix = 'th'+probe+'_'
            themis_read_efield, time_range, probe=probe
            b_var = themis_read_bfield(time_range, probe=probe)
            r_var = themis_read_orbit(time_range, probe=probe)
            themis_read_density, time_range, probe=probe

            get_data, prefix+'b_gsm', times
            interp_time, prefix+'r_gsm', times
            time_step = total(times[0:1]*[-1,1])

            vars = prefix+['e','b']+'_gsm'
            foreach var, vars do begin
                get_data, var, times, data
                mag = snorm(data)
                index = where(finite(mag))
                data = sinterpol(data[index,*], times[index], times)
                store_data, var, times, data
            endforeach
        endforeach

        vars = []
        foreach probe, probes do begin
            prefix = 'th'+probe+'_'
            vars = [vars, prefix+['e_gsm','b_gsm','ele_n']]
        endforeach
        ;tplot, vars, trange=time_range

    ;---Calc ExB vel.
        foreach probe, probes do begin
            prefix = 'th'+probe+'_'
            e_gsm = get_var_data(prefix+'e_gsm', times=times)
            b_gsm = get_var_data(prefix+'b_gsm')
            vexb_gsm = vec_cross(e_gsm,b_gsm)
            coef = 1e3/snorm(b_gsm)^2
            for ii=0,2 do vexb_gsm[*,ii] *= coef
            store_data, prefix+'vexb_gsm', times, vexb_gsm
            add_setting, prefix+'vexb_gsm', /smart, dictionary($
                'display_type', 'vector', $
                'short_name', 'ExB V', $
                'unit', 'km/s', $
                'coord', 'GSM', $
                'coord_labels', constant('xyz') )
    
        endforeach

        stop

   ;---Calculate higher level quantities.
       foreach probe, probes do begin
           prefix = 'th'+probe+'_'
;            define_fac, prefix+'b_gsm', prefix+'r_gsm'
;            foreach var, prefix+['e','b']+'_' do to_fac, var+'gsm', to=var+'fac'
            pf_var = prefix+'pf_gsm'
            stplot_calc_pflux_mor, prefix+'e_gsm', prefix+'b_gsm', pf_var
            add_setting, pf_var, /smart, dictionary($
                'display_type', 'vector', $
                'short_name', 'S', $
                'unit', 'mW/m!U2!N', $
                'colors', rgb, $
                'coord', 'GSM', $
                'coord_labels', xyz )


        ;---Calc E/B ratio.
            vars = prefix+['e','b']+'_gsm'

            ; settings for wavelet transform.
            s0 = 4d*time_step
            dj = 1d/8
            s1 = 2000
            j1 = floor(alog(s1/s0)/alog(2)/dj)
            s1 = s0*2d^(dj*j1)
            ns = j1+1
            w0 = 6d
            cdelta = 0.776d
            psi0 = !dpi^(-0.25)

            foreach tvar, vars do begin
                get_data, tvar, uts, dat
                dat = snorm(dat)
                mor = wavelet(dat, time_step, /pad, $
                    s0=s0, dj=dj, j=j1, mother='Morlet', param=w0, period=ps, scale=ss)
                psd = abs(mor)^2
                idx = where(uts ge onset_time_range[0] and uts le onset_time_range[1], tnrec)
                ;idx = where(uts ge time_range[0] and uts le time_range[1], tnrec)
                psd = psd[idx,*]
                gws = total(psd,1)/tnrec^2
                ngws = (gws/ss)*(time_step*dj/cdelta)*tnrec
                store_data, tvar+'_tmp', ps, [[gws],[ngws]]
            endforeach

            get_data, prefix+'e_gsm_tmp', ps, edat
            get_data, prefix+'b_gsm_tmp', ps, bdat
            ebratio = sqrt(edat[*,0]/bdat[*,0])*1e3

            ; Calc V_alfven.
            va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
            bmag = median(snorm(get_var_data(prefix+'b_gsm', in=onset_time_range)))
            density = median(get_var_data(prefix+'ele_n', in=onset_time_range))
            mhd_rho = density*1
            va = va0*bmag/sqrt(mhd_rho)
            store_data, prefix+'ebratio', ps, ebratio, [bmag,density,va]
        endforeach

        foreach probe, probes do begin
            prefix = 'th'+probe+'_'

            data = get_var_data(prefix+'ele_n', in=time_range)
            yrange = minmax([0.05,5,minmax(data)])
            options, prefix+'ele_n', 'yrange', yrange
            get_data, prefix+'b_gsm', times, b_gsm, limits=lim
            store_data, prefix+'b_gsm1', times, [[snorm(b_gsm)],[b_gsm]], limits=lim
            options, prefix+'b_gsm1', 'labels', ['|B|',lim.labels]
            options, prefix+'b_gsm1', 'colors', [sgcolor('black'),lim.colors]

            get_data, prefix+'pf_gsm', times, pf_gsm
            pf_para = sdot(pf_gsm, sunitvec(b_gsm))
            store_data, prefix+'pf_para', times, pf_para
            add_setting, prefix+'pf_para', /smart, dictionary($
                'display_type', 'scalar', $
                'unit', 'mW/m!U2!N', $
                'short_name', 'S!D||!N in-situ' )


            vars = prefix+['e_gsm','b_gsm1','pf_para','ele_n']
            nvar = n_elements(vars)+1
            ypans = fltarr(nvar)+1 & ypans[-1] = 2
            margins = [12,5,10,4]
            plot_file = join_path([plot_dir,'event_'+$
                strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_th'+probe+'.pdf'])
            if keyword_set(test) then plot_file = 0
            sgopen, plot_file, xsize=5, ysize=8, xchsz=xchsz, ychsz=ychsz
            poss = sgcalcpos(nvar, margin=margins,ypans=ypans)
            poss[3,-1] -= ychsz*4
            tplot, vars, trange=time_range, position=poss[*,0:nvar-1]
            timebar, onset_time_range

            tpos = poss[*,-1]
            get_data, prefix+'ebratio', ps, ebratio, val
            xrange = [0.1,100]
            plot, 1d3/ps, ebratio, $
                xstyle=1, xlog=1, xtitle='Freq (mHz)', xrange=xrange, xticklen=-0.02, $
                ystyle=0, ylog=1, ytitle='E/B ratio (km/s)', yticklen=-0.01, $
                position=tpos, /noerase
            va = val[2]
            plots, xrange, va+[0,0], linestyle=1
            tx = tpos[2]+xchsz*1
            ty = (convert_coord(xrange[0],va, /data, /to_normal))[1]-ychsz*0.4
            msg = 'V!DA!N'
            xyouts, tx,ty,/normal, msg

            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = 'V!DA!N (km/s): '+string(va,format='(I0)')
            xyouts, tx,ty,/normal, msg

            tpos = poss[*,1]
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            bmag = val[0]
            msg = '|B| (nT): '+string(bmag, format='(F6.2)')
            xyouts, tx,ty,/normal, msg

            tpos = poss[*,3]
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            bmag = val[1]
            msg = 'n (cc): '+string(bmag, format='(F6.3)')
            xyouts, tx,ty,/normal, msg

            tpos = poss[*,0]
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            msg = 'TH-'+strupcase(probe)+'    '+$
                time_string(time_range[0],tformat='YYYY-MM-DD/hh:mm')+' to '+$
                time_string(time_range[1],tformat='hh:mm')
            xyouts, tx,ty,/normal, msg

            fig_labels = letters(nvar)+'.'
            for ii=0,nvar-1 do begin
                tpos = poss[*,ii]
                tx = xchsz*4
                ty = tpos[3]-ychsz*0.8
                xyouts, tx,ty,/normal, fig_labels[ii]
            endfor
            if keyword_set(test) then stop
            sgclose
        endforeach
    endforeach

end
