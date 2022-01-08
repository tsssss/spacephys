
; check the radiation dose received by RBSP-A and -B
; during 2012 Oct to 2016 Nov. This is the all the current available data.
; load mageis data. pick the energy around 2 MeV.
; By 2017-02-13.


; default settings.
rootdir = shomedir()+'/rbsp_dose'
if file_test(rootdir) eq 0 then file_mkdir, rootdir

probes = ['a','b']
nprobe = n_elements(probes)
secofday = 86400d


; settings.
utr = time_double(['2012-10-01','2016-12-01'])

en0 = 2000      ; in keV.
en0 = 500
den = 0.15*en0   ; dE in keV.
sen0 = sgnum2str(en0)
dt0 = 3600      ; hour resolution.
reload = 0      ; set to reload the data.
dathrfn = rootdir+'/rb_dose_hour_'+sen0+'kev.tplot'


; all the dates.
uts = smkarthm(utr[0], utr[1], secofday, 'dx')
nut = n_elements(uts)


; download data.
download = 0
if download eq 1 then begin  
    logfn = rootdir+'/download_process.log'     ; save download process.
    if file_test(logfn) eq 0 then stouch, logfn
    
    for i = 0, nut-1 do begin
        tutr = uts[i]+[0,secofday]
        
        printf, -1, 'Downloading '+time_string(tutr[0])+' ...'
        
        openw, loglun, logfn, /get_lun, /append
        printf, loglun, 'Downloading '+time_string(tutr[0])+' ...'
        printf, loglun, ''
        free_lun, loglun
        
        for j = 0, nprobe-1 do begin
            tmp = sread_rbsp_mageis_l3(tutr, probe = probes[j])            
        endfor
    endfor
endif


; extract data at the wanted energy from files.
load = 0
if file_test(dathrfn) eq 0 then load = 1
if reload eq 1 then load = 1

if load then begin
    uts = smkarthm(utr[0], utr[1], dt0, 'dx')
    nrec = n_elements(uts)
    flux = dblarr(nrec)
    ebin = dblarr(nrec)

    for j = 0, nprobe-1 do begin
        ; set tuts smaller than utr[0] to force load new data.
        tuts = utr[0]-dt0

        ; loop through each element in flux.
        for i = 0, nrec-1 do begin
            ; load data when the current buffer is done.
            if uts[i] gt tuts[-1] then begin
                tutr = uts[i]-(uts[i] mod secofday)
                mageisl3 = sread_rbsp_mageis_l3(tutr, probe = probes[j])

                ; skip to next day if no data for current day, 
                if size(mageisl3,/type) ne 8 then begin
                    print, 'no data on '+time_string(tutr[0])+' ...'
                    tmp = uts[i]-(uts[i] mod secofday)+secofday
                    idx = where(uts ge tmp, cnt)
                    if cnt eq 0 then break
                    i = idx[0]
                    continue
                endif

                ; load data for current day.
                tuts = sfmepoch(mageisl3.epoch, 'unix')
                tdat = mageisl3.fedu
                tens = mageisl3.fedu_energy

                tnrec = n_elements(tuts)
                dat = dblarr(tnrec)
                ens = dblarr(tnrec)
                for k = 0, tnrec-1 do begin
                    tmp = min(tens[k,*]-en0, enidx, /absolute, /nan)
                    if abs(tmp) gt den then continue
                    if finite(tens[k,enidx]) eq 0 then continue
                    if tens[k,enidx] lt en0-den then continue
                    dat[k] = tdat[k,enidx]
                    ens[k] = tens[k,enidx]
                endfor
            endif


            idx = where(tuts ge uts[i] and tuts lt uts[i]+dt0, cnt)
            if cnt eq 0 then begin
                flux[i] = !values.d_nan
                ebin[i] = !values.d_nan
            endif else begin
                flux[i] = mean(dat[idx],/nan)
                ebin[i] = mean(ens[idx],/nan)
            endelse
        endfor

        ; save data as tplot var.
        idx = where(flux eq 0, cnt)
        if cnt ne 0 then flux[idx] = !values.d_nan

        tvar = 'rbsp'+probes[j]+'_dose_hour_'+sen0+'kev'
        store_data, tvar, uts, flux, ebin
    endfor

    vars = ['rbsp'+probes]+'_dose_hour_'+sen0+'kev'
    tplot_save, vars, filename = dathrfn
endif


tplot_restore, filename = dathrfn
for j = 0, nprobe-1 do begin
    tvar = 'rbsp'+probes[j]+'_dose_hour_'+sen0+'kev'
    var = 'rbsp'+probes[j]+'_int_dose_hour_'+sen0+'kev'
    get_data, tvar, uts, dat, ens
    idx = where(finite(dat,/nan), cnt)
    if cnt ne 0 then dat[idx] = 0
    idx = where(ens lt en0-den or ens gt en0+den, cnt)
    if cnt ne 0 then begin
        ens[idx] = !values.d_nan
        dat[idx] = 0
    endif
    nrec = n_elements(dat)
    for i = 1, nrec-1 do dat[i] = dat[i]+dat[i-1]
    store_data, var, uts, dat*dt0, ens
    var = 'rbsp'+probes[j]+'_dose_energy_bin'
    store_data, var, uts, ens
endfor

vars = ['rbsp'+probes]+'_int_dose_hour_'+sen0+'kev'
tvar = 'rb_int_dose_hour_'+sen0+'kev'
stplot_merge, vars, newname = tvar, limits = $
    {labels:['RBSP-A','RBSP-B'], colors:[4,6], $
    ytitle:'Integrated Flux!C(1/cm!U2!N-sr-keV)'}
    
vars = ['rbsp'+probes]+'_dose_hour_'+sen0+'kev'
tvar = 'rb_dose_hour_'+sen0+'kev'
stplot_merge, vars, newname = tvar, limits = $
    {labels:['RBSP-A','RBSP-B'], colors:[4,6], $
    ytitle:'Flux!C(1/cm!U2!N-s-sr-keV)'}

vars = ['rbsp'+probes]+'_dose_energy_bin'
tvar = 'rb_dose_energy_bin'
stplot_merge, vars, newname = tvar, limits = $
    {labels:['RBSP-A','RBSP-B'], colors:[4,6], $
    ytitle:'Energy!C(keV)', yrange:en0+[-1,1]*den}



; **** load Dst and AE.
omnifn = rootdir+'/omni.tplot'
if file_test(omnifn) eq 1 then begin
    tplot_restore, filename = omnifn
endif else begin
    dat = sread_omni(utr, vars=['Epoch','SYM_H','AE_INDEX'])
    symh0 = -40
    
    uts = sfmepoch(dat.epoch,'unix')
    store_data, 'symh', uts, dat.sym_h, limits = {constant:[0,symh0], ytitle:'Sym/H!C(nT)'}
    store_data, 'ae', uts, dat.ae_index, limits = {ytitle:'AE!C(nT)'}
    tplot_save, ['symh','ae'], filename = omnifn
endelse
    
ofn = rootdir+'/rb_dose_hour_'+sen0+'kev.pdf'
sgopen, ofn, xsize = 8, ysize = 6, /inch

device, decomposed = 0
loadct2, 43
tplot_options, 'labflag', -1

vars = ['rb_int_dose_hour_'+sen0+'kev', $
    'rb_dose_hour_'+sen0+'kev', $
    'rb_dose_energy_bin', $
    'symh','ae']

title = 'RBSP flux around '+sen0+' keV, '+ $
    time_string(utr[0],tformat='YYYY-MM')+' to '+ $
    time_string(utr[1],tformat='YYYY-MM')

tplot, vars, trange = utr, title = title

sgclose


end
