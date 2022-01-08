pro bias_param_to_step, tvar
    get_data, tvar, t0, dat
    nrec = n_elements(t0)
    t1 = t0[0]
    tmp = dat[0]
    for i = 1, nrec-1 do begin
        t1 = [t1,t0[i],t0[i]]
        tmp = [tmp,dat[i-1],dat[i]]
    endfor
    get_data, tvar, limit = limit
    store_data, tvar, /delete
    store_data, tvar+'_fg', t0, dat
    options, tvar+'_fg', 'psym', 4
    store_data, tvar+'_bg', t1, tmp
    vars = tvar+['_fg','_bg']
    store_data, tvar, data = vars
end

pro bias_param_info, tvar, trange = tr, datarate = dr0
    get_data, tvar, t0, d0
    if n_elements(tr) eq 0 then tr = [min(t0),max(t0)]
    ; analyze data in given time range.
    idx = where(t0 ge tr[0] and t0 le tr[1], nrec)
    t1 = t0[idx]
    d1 = d0[idx]
    ; analyze data at given data rate.
    diff = t1[1:nrec-1]-t1[0:nrec-2]
    if n_elements(dr0) eq 0 then dr0 = mean(diff)
    idx = where(diff le dr0*2, nrec)
    t1 = t1[idx]
    d1 = d1[idx]
    ; round data. val0s are original values, val1s are rounded values.
    val0s = d1[uniq(d1,sort(d1))]
    nval = n_elements(val0s)
    if nval eq 1 then begin
        val1s = val0s[0]
        ts = t1[[0,nrec-1]]
    endif else begin
        diff = val0s[1:nval-1]-val0s[0:nval-2]
        stdv = stddev(diff)
        idx = where(diff ge stdv, nstep)
        if idx[0] eq 0 then nstep-= 1 else idx = [0,idx]
        if idx[nstep] eq nval-1 then nstep-= 1 else idx = [idx,nval-1]
        val1s = val0s
        for i = 0, nstep do val1s[idx[i]:idx[i+1]-1] = mean(val1s[idx[i]:idx[i+1]-1])
        for i = 0, nval-1 do d1[where(d1 eq val0s[i])] = val1s[i]
        ; find where steps are.
        diff = [0,d1[1:nrec-1]-d1[0:nrec-2]]
        idx = where(diff ne 0, nstep)
        if idx[0] eq 0 then nstep-= 1 else idx = [0,idx]
        if idx[nstep] eq nrec-1 then nstep-= 1 else idx = [idx,nrec-1]
        diff = idx[1:nstep]-idx[0:nstep-1]
        tmp = mean(diff)-stddev(diff)
        for i = 0, nstep-1 do if diff[i] lt tmp then idx = [idx[0:i],idx[i+2:*]]
        ts = t1[idx]
    endelse
    nstep = n_elements(ts)-1
    vs = dblarr(nstep)
    for i = 0, nstep-1 do vs[i] = d1[(where(t1 ge (ts[i]+ts[i+1])*0.5))[0]]
    store_data, tvar+'_tstep', ts, vs
end


pro bias_parameter, rbx, probes = probes, trange = tr, erange = er, vrane = vr, $
    reload = reload, date = date, datarate = dr0, guard = guard, usher = usher
        
    ; useful string.
    if n_elements(rbx) eq 0 then rbx = 'a'
    rbspx = 'rbsp'+rbx
    prefix = rbspx+'_'
    if n_elements(probes) ne 2 then probes = [1,2]
    suffix = string(probes,format='(I0)')
    comp = strjoin(suffix,'')
    
    ; time range for plotting.
    if n_elements(tr) ne 2 then tr = timerange()
    
    ; data rate for stepping the parameters.
    if n_elements(dr0) eq 0 then dr0 = 11   ; spin period.
    
    ; info for saving plots.
    if n_elements(date) eq 0 then $
        date = time_string(tr[1], tformat='YYYY_MMDD')
    locdir = '~/bias_param_'+date
    if file_test(locdir,/directory) eq 0 then file_mkdir, locdir  
    dir = locdir+'/RBSP'+strupcase(rbx)+comp
    if file_test(dir,/directory) eq 0 then file_mkdir, dir
    
    ; load data or not.
    load = 0
    vars = prefix+['V'+suffix,['guard','usher','ibias']+comp,'efw_vsvy','efw_esvy']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], data = data
        ; not a combo or structure.
        if n_elements(data) eq 1 and size(data,/type) ne 8 then load = 1
    endfor
    if keyword_set(reload) then load = 1
  
    if load eq 1 then begin
        rbsp_efw_init

        ; load bias current.
        ; rbspx_IEFI_[GUARD,IBIAS,USHER][123456].
        type = 'hsk_beb_analog'
        files = rbspx+'/l1/'+type+'/YYYY/'+rbspx+'_l1_'+type+'_YYYYMMDD_v*.cdf'
        files = file_dailynames(file_format=files,trange=tr)
        files = file_retrieve(files, /last_version, _extra=!rbsp_efw)
        cdf2tplot, file = files, /get_support_data, /convert_int1_to_int2, prefix = prefix
    
        vars = prefix+['CONFIG0','IBEB_TEMP','SPARE','gdh_*','ga_*','ccsds_*','packet_*']
        store_data, vars, /delete
        
        ; convert to step style, and combine the hsk data.
        for i = 0, 3-1 do begin
            tsufx = [string(i*2+1,format='(I0)'),string(i*2+2,format='(I0)')]
            tcomp = strjoin(tsufx,'')
            
            ; find step time, convert to step style.
            vars = tnames(prefix+'IEFI_*'+tsufx)
            for j = 0, n_elements(vars)-1 do begin
                bias_param_info, vars[j], trange = tr, datarate = dr0
                bias_param_to_step, vars[j]
            endfor
            
            ; combine hsk data.
            ; 'rbspx_[guard,usher,ibias]_[12,34,56].
            tvars = [prefix+'IEFI_GUARD'+tsufx[0]+['_fg','_bg'],prefix+'IEFI_GUARD'+tsufx[1]+['_fg','_bg']]
            store_data, prefix+'guard'+tcomp, data = tvars, limit = {labels:'G'+tcomp}
            tvars = [prefix+'IEFI_USHER'+tsufx[0]+['_fg','_bg'],prefix+'IEFI_USHER'+tsufx[1]+['_fg','_bg']]
            store_data, prefix+'usher'+tcomp, data = tvars, limit = {labels:'U'+tcomp}
            tvars = [prefix+'IEFI_IBIAS'+tsufx[0]+['_fg','_bg'],prefix+'IEFI_IBIAS'+tsufx[1]+['_fg','_bg']]
            store_data, prefix+'ibias'+tcomp, data = tvars, limit = {labels:'I'+tcomp}
        endfor

        ; load probe potential.
        ; rbspx_efw_vsvy.
        rbsp_load_efw_waveform, probe = rbx, type = 'calibrated', datatype = ['esvy','vsvy'], trange = tr
        vars = prefix+['efw_esvy_*','efw_vsvy_*']
        store_data, vars, /delete

        ; split probe potentials into components.
        vars = prefix+'efw_vsvy'
        get_data, vars, t0, vsvy
        for i = 0, 6-1 do begin
            tvar = prefix+'V'+string(i+1,format='(I0)')
            store_data, tvar, t0, vsvy[*,i]
            options, tvar, 'ysubtitle', '[V]'
            options, tvar, 'colors', i
            options, tvar, 'labels', 'V'+string(i+1,format='(I0)')
        endfor        
    endif
    
    ; paint color for all components.
    tplot_options, 'labflag', 1  ; center label.
    for i = 0, 1 do options, '*'+suffix[i], 'colors', probes[i]
    
    ; set range for V survey.
    if n_elements(vr) eq 0 then vr = [-30, 5]
    tvar = prefix+'efw_vsvy'
    options, tvar, 'yrange', vr
    options, tvar, 'ystyle', 1
    
    ; set range for E survey.
    if n_elements(er) eq 0 then er = [-10, 10]
    tvar = prefix+'efw_esvy'
    options, tvar, 'yrange', vr
    options, tvar, 'ystyle', 1

    ; final vars to plot.
    vars = prefix+['V'+suffix,['guard','usher','ibias']+comp,'efw_vsvy','efw_esvy']
    options, vars, 'constant', 0
    options, vars, 'yticklen', 1
    options, vars, 'ygridstyle', 1
  
    set_plot, 'x'
    !p.color = 0
    !p.background = 255
    window, xsize = 1000, ysize = 800, xpos = 0, ypos = 600
    title = 'RBSP-'+strupcase(rbx)+' Probe'+comp
    tplot, vars, title = title, trange = tr
    
    ofn = dir+'/bias_param_'+date+ $
        '_RBSP'+strupcase(rbx)+'_Boom'+comp+'.eps'
    set_plot, 'ps'
    device, filename = ofn, /encapsulate, xsize = 10, ysize = 8, /inch, $
        decomposed = 0, /color, bits_per_pixel = 8
    tplot, vars, title = title, trange = tr
    device, /close
    ; ps to pdf.
    spawn, 'cd '+dir+';/opt/local/bin/ps2pdf -dEPSCrop '+ofn
    spawn, 'cd '+dir+';rm -f '+ofn
    set_plot, 'x'
    !p.color = 0
    !p.background = 255
    tplot, vars, title = title, trange = tr
    
    ; subplots on each guard value.
    if keyword_set(guard) then begin
        tdir = dir+'/guard'
        if file_test(tdir,/directory) eq 0 then file_mkdir, tdir
        
        tvar = prefix+'IEFI_GUARD'+suffix[0]
        get_data, tvar+'_tstep', gts, gvals
        
        for gi = 0, n_elements(gts)-2 do begin
            ttr = gts[gi:gi+1]+2*dr0*[-1,1]
            gval = string(round(gvals[(where(gts ge mean(ttr)))[0]-1]),format='(I0)')

            ofn = tdir+'/bias_param_'+date+ $
              '_RBSP'+strupcase(rbx)+'_Boom'+comp+'_Guard='+gval+'.eps'
            set_plot, 'ps'
            device, /encapsulate, xsize = 10, ysize = 8, /inch, $
              decomposed = 0, /color, bits_per_pixel = 8, $
              filename = ofn
            ttitle = title+' Guard = '+gval
            tplot, vars, title = ttitle, trange = ttr
            device, /close
            ; ps to pdf.
            spawn, 'cd '+tdir+';/opt/local/bin/ps2pdf -dEPSCrop '+ofn
            spawn, 'cd '+tdir+';rm -f '+ofn
            set_plot, 'x'
            !p.color = 0
            !p.background = 255
            tplot, vars, title = ttitle, trange = ttr
        endfor
    endif

    if keyword_set(usher) then begin
        tdir = dir+'/usher'
        if file_test(tdir,/directory) eq 0 then file_mkdir, tdir
        
        tvar = prefix+'IEFI_USHER'+suffix[0]
        get_data, tvar+'_tstep', uts, uvals
        
        tvar = prefix+'IEFI_GUARD'+suffix[0]
        get_data, tvar+'_tstep', gts, gvals

        for ui = 0, n_elements(uts)-2 do begin
            ttr = uts[ui:ui+1]+2*dr0*[-1,1]
            uval = string(round(uvals[(where(uts ge mean(ttr)))[0]-1]),format='(I0)')
            gval = string(round(gvals[(where(gts ge mean(ttr)))[0]-1]),format='(I0)')

            ofn = tdir+'/bias_param_'+date+ $
              '_RBSP'+strupcase(rbx)+'_Boom'+comp+'_Guard='+gval+'_Usher='+uval+'.eps'
            set_plot, 'ps'
            device, /encapsulate, xsize = 10, ysize = 8, /inch, $
              decomposed = 0, /color, bits_per_pixel = 8, $
              filename = ofn
            ttitle = title+', Guard = '+gval+', Usher = '+uval
            tplot, vars, title = ttitle, trange = ttr
            device, /close
            ; ps to pdf.
            spawn, 'cd '+tdir+';/opt/local/bin/ps2pdf -dEPSCrop '+ofn
            spawn, 'cd '+tdir+';rm -f '+ofn
            set_plot, 'x'
            !p.color = 0
            !p.background = 255
            tplot, vars, title = ttitle, trange = ttr
        endfor
    endif
    
    tplot, vars, title = title, trange = tr
  
end

;timespan, '2014-12-03', 1
;
;trb12 = time_double(['2014-12-03/01:25','2014-12-03/02:05'])
;e0b12 = 0.1126
;;bias_parameter, 'b', probes = [1,2], trange = trb12, offset = e0b12
;
;trb34 = time_double(['2014-12-03/02:12','2014-12-03/02:47'])
;e0b34 = -4.583
;;bias_parameter, 'b', probes = [3,4], trange = trb34, offset = e0b34


timespan, '2014-12-16', 2
date = '2014_1217'

tra12 = time_double(['2014-12-17/06:12','2014-12-17/06:40'])
bias_parameter, 'a', probes = [1,2], trange = tra12, date = date, /guard;, /usher

tra34 = time_double(['2014-12-17/06:57','2014-12-17/07:26'])
bias_parameter, 'a', probes = [3,4], trange = tra34, date = date, /guard;, /usher

trb12 = time_double(['2014-12-16/23:25','2014-12-16/23:55'])
bias_parameter, 'b', probes = [1,2], trange = trb12, date = date, /guard;, /usher

trb34 = time_double(['2014-12-17/00:10','2014-12-17/00:42'])
bias_parameter, 'b', probes = [3,4], trange = trb34, date = date, /guard;, /usher

trb56 = time_double(['2014-12-17/00:57','2014-12-17/01:25'])
bias_parameter, 'b', probes = [5,6], trange = trb56, date = date, /guard;, /usher

;timespan, '2015-11-14', 1
;
;trb12 = time_double(['2015-11-14/05:25','2015-11-14/06:05'])
;bias_parameter, 'b', probes = [1,2], trange = trb12, /guard, /usher
;
;trb34 = time_double(['2015-11-14/06:12','2015-11-14/06:47'])
;bias_parameter, 'b', probes = [3,4], trange = trb34, /guard, /usher
;
;trb56 = time_double(['2015-11-14/06:55','2015-11-14/07:35'])
;bias_parameter, 'b', probes = [5,6], trange = trb56, /guard, /usher
;
;tra12 = time_double(['2015-11-14/06:27','2015-11-14/07:05'])
;bias_parameter, 'a', probes = [1,2], trange = tra12, /guard, /usher
;
;tra34 = time_double(['2015-11-14/07:10','2015-11-14/07:47'])
;bias_parameter, 'a', probes = [3,4], trange = tra34, /guard, /usher
;
;tra56 = time_double(['2015-11-14/07:57','2015-11-14/08:32'])
;bias_parameter, 'a', probes = [5,6], trange = tra56, /guard, /usher


end
