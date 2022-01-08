;+
; load bias data at given date for certain probe.
;-

pro bias_param_load_data, date, probes = probes

    if n_elements(date) eq 0 then message, 'no date info ...'
    ut0 = time_double(date)

    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    nboom = 6
    suffix = ['1','2','3','4','5','6']
    pairs = ['12','34','56']
    spinrate = 12   ; sec.
    maxni = 20      ; max # of steps in ibias.
    dibias = 5      ; nA, error in ibias.
    dusher = 1      ; V.
    dguard = 1      ; V.

    dt0 = 86400d    ; sec.

    ; **** load data.
    ; rbspx_IEFI_[IBIAS,GUARD,USHER].
    ; rbspx_V[123456], rbspx_[ibias,guard,usher][123456].
    ; rbspx_IV[123456], rbspx_efw_[esvy].

    ; load 2 days of data, b/c a sweep can sit on the edge of 2 days.
    timespan, ut0-dt0, 2, /day
    rbsp_efw_init

    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        
        ; load bias current.
        ; load the variables by ourself, to get cleaner tplot var.
        ; rbspx_IEFI_[GUARD,IBIAS,USHER][123456].
        type = 'hsk_beb_analog'
        files = 'rbsp'+probes[i]+'/l1/'+type+'/YYYY/'+pre0+'l1_'+type+'_YYYYMMDD_v*.cdf'
        files = file_dailynames(file_format=files,trange=tr)
        files = file_retrieve(files, /last_version, _extra=!rbsp_efw)
        cdf2tplot, file = files, /get_support_data, /convert_int1_to_int2, prefix = pre0
    
        vars = pre0+['CONFIG0','IBEB_TEMP','SPARE','gdh_*','ga_*','ccsds_*','packet_*']
        store_data, vars, /delete   ; delete vars of no use.
        
        ; convert from rbspx_IEFI_XXX[123456] to rbspx_xxx[123456].
        ; e.g., rbspx_IEFI_GUARD1 to rbspx_guard1.
        for j = 0, nboom-1 do begin
            vars = tnames(pre0+'IEFI_*'+suffix[j])
            for k = 0, n_elements(vars)-1 do begin
                tvar = strlowcase(pre0+strmid(vars[k],strpos(vars[k],'_',/reverse_search)+1))
                get_data, vars[k], t0, dat
                store_data, tvar, t0, dat
            endfor
        endfor

        ; combine rbspx_IEFI_XXX[123456] into rbspx_IEFI_XXX.
        ; these vars are for backup purpose.
        vars = pre0+'IEFI_'+['IBIAS','GUARD','USHER']
        for j = 0, n_elements(vars)-1 do begin
            tvar = tnames(vars[j]+'?')
            tmp = []
            for k = 0, nboom-1 do begin
;                bias_param_to_step, tvar[k]
                get_data, tvar[k], t0, dat
                tmp = [[tmp],[dat]]
            endfor
            store_data, tvar, /delete
            store_data, vars[j], t0, tmp
        endfor

        ; load probe potential.
        ; rbspx_efw_esvy.
        rbsp_load_efw_waveform, probe = probes[i], type = 'calibrated', datatype = ['esvy','vsvy']
        vars = pre0+['efw_esvy_*','efw_vsvy_*']
        store_data, vars, /delete

        ; split probe potentials into components.
        vars = pre0+'efw_vsvy'
        get_data, vars, t0, vsvy
        for j = 0, 6-1 do begin
            tvar = pre0+'V'+string(j+1,format='(I0)')
            store_data, tvar, t0, vsvy[*,j]
            options, tvar, 'ysubtitle', '[V]'
            options, tvar, 'colors', j
            options, tvar, 'labels', 'V'+string(j+1,format='(I0)')
        endfor
        store_data, 'efw_vsvy', /delete
            
        ; get the i-v relation.
        ; rbsp[ab]_[IV][123456].
        for j = 0, nboom-1 do begin ; each boom.
            ; find bias sweep time.
            sdtutr = [0d,0]
            
            ; get bias times.
            bias_param_get_sdt_times, $
                pre0+['ibias','guard','usher']+suffix[j], sdtutr
            print, time_string(sdtutr)
            if sdtutr[1]-sdtutr[0] le 600 then begin    ; one sweep is ~20 min.
                tvar = pre0+'IV'+suffix[j]
                store_data, tvar, 0, -1
                continue
            endif
            
            ; trim data to sweep time.
            get_data, pre0+'ibias'+suffix[j], t0, ibias
            get_data, pre0+'guard'+suffix[j], t0, guard
            get_data, pre0+'usher'+suffix[j], t0, usher

            idx = where(t0 ge sdtutr[0] and t0 le sdtutr[1], cnt)
            t0 = t0[idx]
            ibias = ibias[idx]
            guard = guard[idx]
            usher = usher[idx]

            ; expand the last record.
            t0 = [t0,t0[-1]+spinrate]
            ibias = [ibias,ibias[-1]]
            guard = [guard,guard[-1]]
            usher = [usher,usher[-1]]
            sdtutr = minmax(t0)
            
            ; round the value to standard values, only for guard and usher.
            store_data, pre0+'ibias'+suffix[j], t0, ibias
            store_data, pre0+'guard'+suffix[j], t0, round(guard)
            store_data, pre0+'usher'+suffix[j], t0, round(usher)
            
            ; average probe potential.
            tvar = pre0+'V'+suffix[j]
            get_data, tvar, t1, vboom
            nrec = n_elements(t0)
            vs = dblarr(nrec)
            for k = 0, nrec-2 do begin
                idx = where(t1 ge t0[k] and t1 lt t0[k+1], cnt)
                vs[k] = (cnt eq 0)? !values.d_nan: mean(vboom[idx])
            endfor
            infos = {utr:sdtutr, i:ibias, v:vs, guard:guard, usher:usher}

            tvar = pre0+'IV'+suffix[j]
            store_data, tvar, 0, infos
        endfor
    endfor
end

date = '2015-11-14'
date = '2013-02-27'
date = '2014-12-17'
bias_param_load_data, date
end
