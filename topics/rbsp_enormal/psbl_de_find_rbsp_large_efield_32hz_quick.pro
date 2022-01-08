; rbspx_de_[mgse,gse,fac,uvw], rbspx_db_[gse,fac], rbspx_[b,bmod]_gse.
; rbspx_[pos_gsm,dis], rbspx_vsc, rbspx_de_flag.
pro psbl_de_load_efield, utr, tprobe, $
    load_spice = load_spice, nodata = nodata, perigeelim = perigeelim

    re = 6374d & re1 = 1d/re
    rgb = sgcolor(['red','green','blue'])


    pre0 = 'rbsp'+tprobe+'_'

    if keyword_set(load_spice) then $
        rbsp_load_spice_kernels, trange = utr, probes = tprobe

    ; load e_despun, 32 Hz.
    ; rbspx_[de_mgse,de].
    tvar = pre0+'de_mgse'
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'esvy')
    if size(dat,/type) ne 8 then begin
        nodata = 1
        return
    endif else nodata = 0
    if size(dat.epoch,/n_elements) le 10 then begin
        nodata = 1
        return
    endif else nodata = 0
    tmp = sfmepoch(dat.epoch, 'unix', /epoch16)
    dr = sdatarate(tmp[sort(tmp)])
    uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
    dat = sinterpol(dat.efield_mgse, tmp, uts)
    dat[*,0] = 0.   ; set x to 0.
    store_data, tvar, uts, dat, limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
    tvar = pre0+'de'
    store_data, tvar, uts, sqrt(total(dat^2,2))



    ; load position.
    ; rbspx_[pos_gsm,dis].
    tvar = pre0+'pos_gsm'
    rbsp_load_spice_state, probe = tprobe, $
        coord = 'gsm', /no_spice_load, trange = utr
    get_data, pre0+'state_pos_gsm', tmp, pos
    pos = sinterpol(pos*re1, tmp, uts)
    store_data, tvar, uts, pos, limits = $
        {ytitle:'R GSM (Re)', colors:rgb, labels:'GSM '+['x','y','z']}
    
    tvar = pre0+'dis'
    store_data, tvar, uts, sqrt(total(pos^2,2)), limits = $
        {ytitle:'Dis (Re)'}
    
    store_data, pre0+'state_*', /delete


    ; load euvw, vsc.
    ; rbspx_[euvw,vsc].
    tvar = pre0+'de_uvw'
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'euvw')
    if size(dat,/type) ne 8 then begin
        nodata = 2
        return
    endif else nodata = 0
    uts = sfmepoch(dat.epoch, 'unix', /epoch16)
    dat.efield_uvw[*,2] = !values.f_nan
    store_data, tvar, uts, dat.efield_uvw, limits = $
        {ytitle:'dE UVW!C(mV/m)', labels:['u','v','w'], colors:rgb, constant:0}

;    tvar = pre0+'vsc'
;    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
;    if size(dat,/type) ne 8 then begin
;        nodata = 3
;        return
;    endif else nodata = 0
;    uts = sfmepoch(dat.epoch, 'unix', /epoch16)
;    tmp = dat.vsvy
;    store_data, tvar, uts, 0.5*(tmp[*,2]+tmp[*,3]), limits = $
;        {ytitle:'Vsc/L 34!C(mV/m)'}


    ; get flags for bad data using various methods.
    padtime = 5*60d ; sec. padd time for bad data flag.
    rbsp_efw_euvw_flag, utr, probes = tprobe, pad = padtime, evar = 'de_uvw'
    rbsp_efw_perigee_flag, utr, probes = tprobe, pad = padtime, /re, mindis = perigeelim
;    rbsp_efw_vsc_flag, utr, probes = tprobe, pad = padtime
        
    ; combine the flags into one flag.
    get_data, pre0+'de_mgse', t0
    tflags = bytarr(n_elements(t0))
;    vars = pre0+[['euvw','perigee','vsc']+'_flag']
    vars = pre0+[['euvw','perigee']+'_flag']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], tmp, dat
        tflags = tflags or interpol(dat, tmp, t0)
    endfor
    store_data, vars, /delete
    tvar = pre0+'de_flag'
    store_data, tvar, t0, tflags, limits = {yrange:[-0.5,1.5], ytitle:'Bad E flag'}
    
    ; mask bad data using the flag.
    idx = where(tflags eq 1, cnt)
    if cnt ne 0 then begin
        vars = pre0+['de_mgse','de']
        for i = 0, n_elements(vars)-1 do begin
            get_data, vars[i], tmp, dat
            dat[idx,*] = !values.f_nan
            store_data, vars[i], tmp, dat
        endfor
    endif

end


pro psbl_de_find_rbsp_large_efield_32hz_quick, utr, tprobe

    if n_elements(utr) eq 0 then message, 'no time range ...'
    if n_elements(tprobe) eq 0 then message, 'no probe ...'

    ; initial settings and constants.
    rbsp_efw_init
    cdf_leap_second_init
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10


    emax = 25  ; mV/m. |E| larger than this threshold will be examined.
    tmax = 60  ; sec. |E| spikes within 1 min will be counted as 1 event.
    tmin = 30  ; sec. min duration of the event.
    perigeelim = 4d ; re. limit for perigee region.
    spc = '    '
    secofday = 86400d
    rgb = [6,4,2]
    rbx = 'rbsp'+tprobe
    pre0 = rbx+'_'
    pre1 = 'rb'+tprobe+'_'  ; prefix for data buffer.
    bufvars = ['de','de_uvw','pos_gsm','dis']
    nbufvar = n_elements(bufvars)
    store_data, bufvars, /delete    ; clear buffer.
    
    rootdir = shomedir()+'/psbl_de/rbsp'+tprobe
    proclog = shomedir()+'/psbl_de/'+pre0+'process.log'     ; record important steps during the loop.
    if file_test(proclog) eq 0 then stouch, proclog

    ; log file to save info of potential events.
    logfn = rootdir+'/list_rbsp'+tprobe+'_large_efield_32hz_quick.log'
    if file_test(logfn) eq 0 then begin
        tmp = file_dirname(logfn)
        if file_test(tmp,/directory) eq 0 then file_mkdir, tmp
        openw, lun, logfn, /get_lun
        tit = '   event id   '+spc+'  start & end time  '+spc+ $
            'dE max'+spc+'      Pos (Re)      '
        printf, lun, tit
        tit = 'YYYY_MMDD_hhmm'+spc+'hh:mm:ss'+spc+'hh:mm:ss'+spc+ $
            '(mV/m)'+spc+'Xgsm'+spc+'Ygsm'+spc+'Zgsm'
        printf, lun, tit
        tit = '--------------'+spc+'--------'+spc+'--------'+spc+ $
            '------'+spc+'----'+spc+'----'+spc+'----'
        printf, lun, tit
        free_lun, lun
    endif

    
    ; important vars related to time.
    etr = stoepoch(utr, 'unix')
    ut1 = utr[0]-(utr[0] mod secofday)  ; start time of a day.

    ; load spice kernel all together.
    rbsp_load_spice_kernels, trange = utr, probes = tprobe


    ; loop through each orbit. 
    go = 1
    load = 1    ; 1 for load a new day into buffer.
    while go do begin
        ; load new data into data buffer if needed.
        if load eq 1 then begin
            openw, loglun, proclog, /get_lun, /append
            printf, loglun, ''
            printf, loglun, '**** loading RBSP-'+strupcase(tprobe)+$
                ', '+time_string(ut1)+' ...'
            free_lun, loglun
            ; load data for a whole day.
            ; get rbspx_[de,de_mgse,pos_gsm,dis,de_uvw,de_flag].
            tutr = ut1+[0,secofday]
            psbl_de_load_efield, tutr, tprobe, nodata = nodata, perigeelim = perigeelim
            
            ; if no data, move to the next day. otherwise add new data to buffer.
            if nodata eq 1 then begin
                openw, loglun, proclog, /get_lun, /append
                printf, loglun, spc+'no E on '+time_string(ut1)+' ...'
                free_lun, loglun
            endif else if nodata eq 2 then begin
                openw, loglun, proclog, /get_lun, /append
                printf, loglun, spc+'no euvw on '+time_string(ut1)+' ...'
                free_lun, loglun
            endif else begin
                for i = 0, nbufvar-1 do begin
                    if tnames(pre1+bufvars[i]) eq '' then begin
                        get_data, pre0+bufvars[i], t0, dat, limits = lim
                        store_data, pre1+bufvars[i], t0, dat, limits = lim
                    endif else begin
                        get_data, pre1+bufvars[i], t1, dat
                        get_data, pre0+bufvars[i], t0, tmp, limits = lim
                        store_data, pre1+bufvars[i], [t1,t0], [dat,tmp], limits = lim
                    endelse
                endfor
            endelse
        endif

        ; get the indices for the first available orbit: idx0, idx1.
        get_data, pre1+'dis', tmp, dis
        idx = where(dis gt perigeelim, cnt)
        if cnt eq 0 then begin  ; all points are around perigee, so load more.
            ut1 = ut1+secofday
            load = 1
            continue
        endif
        idx0 = idx[0]   ; the idx of 1st point out of perigee region.
        idx = where(dis[idx0:*] le perigeelim, cnt)
        if cnt eq 0 then begin
            ut1 = ut1+secofday
            load = 1
            continue
        endif
        idx1 = idx[0]+idx0  ; the idx of 1st point of the nex perigee region.
        orbutr = tmp[[idx0,idx1]]   ; the start and end times of current orbit.
        openw, loglun, proclog, /get_lun, /append
        printf, loglun, spc+$
            'orbit: '+time_string(orbutr[0])+' to '+time_string(orbutr[1])
        free_lun, loglun

        
        ; move the data in current orbit out from the buffer,
        ; then save them to rbspx_xxx.
        for i = 0, nbufvar-1 do begin
            get_data, pre1+bufvars[i], t0, dat, limits = lim
            idx = where(t0 ge orbutr[0] and t0 lt orbutr[1])
            store_data, pre0+bufvars[i], t0[idx], dat[idx,*], limits = lim
            idx = where(t0 ge orbutr[1])
            store_data, pre1+bufvars[i], t0[idx], dat[idx,*], limits = lim
        endfor


        ; check large e field.
        get_data, pre0+'de', tuts, tdemag
        dr = sdatarate(tuts)
        idx = where(tdemag ge emax, cnt)
        if cnt eq 0 then begin
            openw, loglun, proclog, /get_lun, /append
            printf, loglun, spc+spc+'no |dE| > '+sgnum2str(emax)+' mV/m ...'
            printf, loglun, ''
            free_lun, loglun
        endif
        eidx0 = 0
        eidx1 = 0
        while eidx0 lt cnt-1 do begin
            ; the large e field within tmax are grouped into one chunk.
            tmp = eidx0
            for i = eidx0, cnt-2 do if idx[i+1]-idx[i] gt tmax/dr then break
            eidx1 = i
            tidx0 = idx[eidx0]
            tidx1 = idx[eidx1]
            deutr = tuts[[tidx0,tidx1]]
            tmp = deutr[1]-deutr[0]
            openw, loglun, proclog, /get_lun, /append
            printf, loglun, spc+spc+'event: '+time_string(deutr[0])+' to '+ $
                time_string(deutr[1])+', duration: '+sgnum2str(tmp)+' sec ...'
            free_lun, loglun
            if tmp lt tmin then begin   ; too short, check the next loop.
                eidx0 = eidx1+1
                openw, loglun, proclog, /get_lun, /append
                printf, loglun, spc+spc+'shorter than '+sgnum2str(tmin)+' sec ...'
                printf, loglun, ''
                free_lun, loglun
                continue
            endif
            openw, loglun, proclog, /get_lun, /append
            printf, loglun, ''
            free_lun, loglun
            
            ; generate plots to file.
            ofn = rootdir+time_string(tuts[tidx0], tformat='/YYYY/MM/')+ $
                pre0+'large_efield_32Hz_'+ $
                time_string(tuts[tidx0], tformat='YYYY_MMDD_hhmm')+'.pdf'
            tutr = deutr+[-1,1]*60*2    ; 2 min pad.
            
            tvar = pre0+'pos_gsm'
            stplot_split, tvar, newname = tvar+'_'+['x','y','z']
            tvar = pre0+'pos_gsm_x' & options, tvar, 'ytitle', 'X GSM (Re)'
            tvar = pre0+'pos_gsm_y' & options, tvar, 'ytitle', 'Y GSM (Re)'
            tvar = pre0+'pos_gsm_z' & options, tvar, 'ytitle', 'Z GSM (Re)'

            options, pre0+'de', 'constant', emax

            pos1 = [0.2,0.50,0.9,0.90]
            pos2 = [0.2,0.10,0.9,0.30]

;ofn = 0
            sgopen, ofn, xsize = 8, ysize = 11.5, /inch

            vars = pre0+['de_mgse','de']
            nvar = n_elements(vars)
            labs = pre0+['pos_gsm_'+['z','y','x']]
            titl = 'RBSP-'+strupcase(tprobe)+' dE survey '+$
                time_string(tuts[tidx0])
            poss = sgcalcpos(nvar, position = pos1)
            tplot, vars, var_label = labs, trange = tutr, title = titl, $
                /noerase, position = poss
            
            ; zoom in and show Euvw.
            timebar, deutr, color = sgcolor('black'), linestyle = 1
            zoominutr = mean(deutr)+[-1,1]*0.5*tmin
            timebar, zoominutr, color = sgcolor('red')
            
            tvar = pre0+'de_uvw'
            tplot, tvar, trange = zoominutr, $
                title = 'Zoom in dE UVW, dt = '+sgnum2str(tmin)+' sec', /noerase, position = pos2
            
            sgclose
;stop
            
            ; prepare output to log file.
            cmd = time_string(tuts[tidx0], tformat='YYYY_MMDD_hhmm')+spc
            cmd+= time_string(tuts[tidx0], tformat='hh:mm:ss')+spc
            cmd+= time_string(tuts[tidx1], tformat='hh:mm:ss')+spc
            cmd+= string(max(tdemag[tidx0:tidx1],maxdeidx),format='(F6.1)')+spc
            get_data, pre0+'pos_gsm', tmp, tposgsm
            cmd+= string(tposgsm[maxdeidx,0], format='(F4.1)')+spc
            cmd+= string(tposgsm[maxdeidx,1], format='(F4.1)')+spc
            cmd+= string(tposgsm[maxdeidx,2], format='(F4.1)')+spc
            
            openw, lun, logfn, /get_lun, /append
            printf, lun, cmd
            free_lun, lun

            eidx0 = eidx1+1
        endwhile
        
        
        ; update data buffer.
        ; trim the near earth field.
        get_data, pre1+'dis', tmp, dis
        idx = where(dis gt perigeelim, cnt)
        ut2 = tmp[idx[0]]   ; the time when sc exits perigee region.
        if cnt eq 0 then begin
            ut1 = ut1+secofday
            load = 1
            continue
        endif else begin
            for i = 0, nbufvar-1 do begin
                get_data, pre1+bufvars[i], t1, dat
                idx = where(t1 gt ut2)
                store_data, pre1+bufvars[i], t1[idx], dat[idx,*]
            endfor
        endelse
        
        ; add a new day when less than 9 hour.
        if (ut2 mod secofday) gt 15d*3600 then begin
            get_data, pre1+'dis', tmp
            tmp = max(tmp) & tmp = tmp-(tmp mod secofday)
            if tmp lt ut1 then begin    ; load a new day of data.
                ut1 = ut1+secofday  ; move to the next day.
                load = 1
                continue    ; this continue will load data in next day.
            endif else load = 0 ; new day of data already loaded.
        endif else load = 0
        
        if ut1 ge max(utr) then break
    endwhile

end

tr = ['2012-09-25','2016-01-01']
probes = ['b']
utr = time_double(tr)
store_data, '*', /delete
foreach tprobe, probes do $
    psbl_de_find_rbsp_large_efield_32hz_quick, utr, tprobe

end
