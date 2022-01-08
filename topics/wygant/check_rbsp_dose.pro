
; check the radiation does received by RBSP-A and -B during 2013 to 2017.
; load mageis data. pick the energy around 2 MeV.


rootdir = shomedir()+'/rbsp_dose'
if file_test(rootdir) eq 0 then file_mkdir, rootdir
logfn = rootdir+'/download_process.log'     ; save download process.

probes = ['a','b']
nprobe = n_elements(probes)

; load data.
utr = time_double(['2012-10-01','2016-12-31'])
secofday = 86400d
uts = smkarthm(utr[0], utr[1], secofday, 'dx')
nut = n_elements(uts)



; download data.
download = 0
if download eq 1 then begin  
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



dathrfn = rootdir+'/rb_does_hour.tplot'



; preprocess data from file. extract 2 MeV from file and save monthly data.
en0 = 2275d     ; ~2 MeV.
den = 0.1*en0   ; dE in eV.
preprocess = 0
if preprocess then begin
    for j = 0, nprobe-1 do begin
        ets = []
        dat = []
        ens = []
        for i = 0, nut-1 do begin
            tutr = uts[i]+[0,secofday]
            tmp = sread_rbsp_mageis_l3(tutr, probe = probes[j])
            if size(tmp,/type) ne 8 then continue
            tets = tmp.epoch
            tens = tmp.fedu_energy
            tdat = tmp.fedu
            nrec = n_elements(tets)
            for k = 0, nrec-1 do begin
                tmp = min(tens[k,*]-en0, enidx, /absolute, /nan)
                if abs(tmp) gt den then continue
                if finite(tens[k,enidx]) eq 0 then continue
                if tens[k,enidx] lt en0-den then continue
                ets = [ets,tets[k]]
                ens = [ens,tens[k,enidx]]
                dat = [dat,tdat[k,enidx]]
;                print, k, tens[k,enidx]
            endfor
            
            tet = tets[nrec-1]
            if abs(tet-sepochceil(tet,'mo')) le 6d5 then begin ; split out per month, when end time is within 10 min of a whole month.
                if n_elements(ets) eq 0 then continue
                tvar = 'rbsp'+probes[j]+'_dose'
                tuts = sfmepoch(ets, 'unix')
                store_data, tvar, tuts, dat, ens
                
                ofn = rootdir+'/'+tvar+sfmepoch(tets[0],'_YYYY_MM')+'.tplot'
                tplot_save, tvar, filename = ofn
                ets = []
                dat = []
                ens = []
            endif
        endfor
    endfor
    if file_test(dathrfn) eq 1 then file_delete, dathrfn
endif



; down-sample, by integrating per hour.
dt = 3600
if file_test(dathrfn) eq 1 then begin
    tplot_restore, filename = dathrfn
endif else begin
    
    ; load data from monthly data.    
    uts = smkarthm(utr[0], utr[1]-dt, dt, 'dx')
    nrec = n_elements(uts)
    dat = dblarr(nrec)
    ens = dblarr(nrec)
    
    for j = 0, nprobe-1 do begin
        
        ; set tuts smaller than utr[0] to force load new data.
        tuts = utr[0]-dt
        
        for i = 0, nrec-1 do begin
            ; load data when the current buffer is done.
            if uts[i] gt tuts[-1] then begin
                tfn = rootdir+'/rbsp'+probes[j]+'_dose'+ $
                    time_string(uts[i],tformat='_YYYY_MM')+'.tplot'
                if file_test(tfn) eq 1 then begin
                    print, 'loading '+tfn+' ...'
                    tplot_restore, filename = tfn
                    tvar = 'rbsp'+probes[j]+'_dose'
                    get_data, tvar, tuts, tdat, tens
                endif else begin
                    print, 'skip '+tfn+' ...'
                    ; move to next month.
                    tmp = sepochceil(stoepoch(uts[i],'unix'),'mo')
                    idx = where(uts ge tmp, cnt)
                    if cnt eq 0 then break
                    i = idx[0]
                endelse
            endif
            
            idx = where(tuts ge uts[i] and tuts lt uts[i]+dt, cnt)
            if cnt eq 0 then begin
                dat[i] = !values.d_nan
                ens[i] = !values.d_nan
            endif else begin
                dat[i] = mean(tdat[idx],/nan)
                ens[i] = mean(tens[idx],/nan)
            endelse
        endfor
        store_data, tvar, uts, dat, ens
    endfor
    
    vars = ['rbsp'+probes]+'_dose'
    tplot_save, vars, filename = dathrfn
endelse

for j = 0, nprobe-1 do begin
    tvar = 'rbsp'+probes[j]+'_dose'
    var = 'rbsp'+probes[j]+'_int_dose'
    get_data, tvar, uts, dat, ens
    idx = where(finite(dat,/nan), cnt)
    if cnt ne 0 then dat[idx] = 0
    nrec = n_elements(dat)
    for i = 1, nrec-1 do dat[i] = dat[i]+dat[i-1]
    store_data, var, uts, dat, ens
endfor

vars = ['rbsp'+probes]+'_int_dose'
tvar = 'rb_int_dose'
stplot_merge, vars, newname = tvar, limits = $
    {labels:['RBSP-A','RBSP-B'], colors:[4,6]}
    
vars = ['rbsp'+probes]+'_dose'
tvar = 'rb_dose'
stplot_merge, vars, newname = tvar, limits = $
    {labels:['RBSP-A','RBSP-B'], colors:[4,6], yrange:[1,1e3]}


    
ofn = 0
sgopen, ofn, xsize = 8, ysize = 6, /inch

device, decomposed = 0
loadct2, 43
tplot_options, 'labflag', -1

tplot, 'rb_'+['int_dose','dose'], trange = utr

sgclose


end