
; poynting flux vs photon count rate.
; etr must be in the same day.

pro pfvscr, etr, sites, probes, del = del

    pfvar = ''
    posvar = ''
    type = 'asf'
    asiloc = spreproot('themis')
    asirem = 'http://themis.ssl.berkeley.edu/data/themis'
    utr = sfmepoch(etr,'unix')
    
    ; how fussy, deg x deg square.
    if n_elements(del) eq 0 then del = 1.5  ; deg.
    
    if n_elements(sites) eq 0 then $
        sites = ['atha','chbg','ekat','fsmi','fsim','fykn',$
            'gako','gbay','gill','inuv','kapu','kian',$
            'kuuj','mcgr','pgeo','pina','rank','snkq',$
            'tpas','whit','yknf','nrsq','snap','talo']
    nsite = n_elements(sites)
    
    if n_elements(probes) eq 0 then probes = ['a','b']
    probes = strlowcase(probes)
    
    ; collect times for aurora image.
    flags = bytarr(nsite)+1b    ; 1: have data in etr.
    asiuts = []
    for i = 0, nsite-1 do begin
        tsite = sites[i]
        asifn = sptn2fn_thm_asi(etr[0], tsite, asiloc, type = type)
        asifn = sprepfile(asifn, asiloc, asirem)
        if asifn[0] eq '' then begin
            flags[i] = 0b & continue
        endif
        tmp = scdfread(asifn, strjoin(['thg',type,tsite,'time'],'_'))
        tuts = *(tmp[0].value)
        idx = where(tuts ge utr[0] and tuts le utr[1], cnt)
        if cnt eq 0 then begin
            flags[i] = 0b & continue
        endif
        asiuts = [asiuts,tuts]
        asiuts = asiuts[sort(asiuts)]
        asiuts = asiuts[uniq(asiuts)]
    endfor
    
    nrec = n_elements(asiuts)
    if nrec eq 0 then message, 'no site has data ...'
    
    sites = sites[where(flags eq 1b)]
    nsite = n_elements(sites)
    print, 'asi available at ', sites
    
    foreach tprobe, probes do begin
        pre1 = 'rbsp'+tprobe+'_'
        
        ; read s/c footpoint.
        tvar = pre1+'fpt_lonlat'     ; s/c footpoint in mlon/mlat in deg.
        get_data, tvar, t0, dat
        scmlon = interpol(dat[*,0], t0, asiuts, /spline)
        scmlat = interpol(dat[*,1], t0, asiuts, /spline)
        
        ; track each site.
        foreach tsite, sites do begin
            ; read pixel mlat, mlon in deg.
            ascfn = sptn2fn_thm_asi(etr[0], tsite, asiloc, type = 'asc')
            ascfn = sprepfile(ascfn, asiloc, asirem)
            tvar = strjoin(['thg',type,tsite,''],'_')+['mlat','mlon']
            tmp = scdfread(ascfn, tvar, 1)
            asimlon = *(tmp[1].value)
            asimlat = *(tmp[0].value)
            asimlon = reform(asimlon[1,*,*])
            asimlat = reform(asimlat[1,*,*])
            asimlon = asimlon[*]        ; make 1-d.
            asimlat = asimlat[*]

            ; read all epoch, find which are within trange.
            asifn = sptn2fn_thm_asi(etr[0], tsite, asiloc, type = type)
            asifn = sprepfile(asifn, asiloc, asirem)
            tmp = scdfread(asifn, strjoin(['thg',type,tsite,'time'],'_'))
            tuts = *(tmp[0].value)
            idx = where(tuts ge utr[0] and tuts le utr[1], cnt)
            rec0 = idx[0]
            rec1 = idx[cnt-1]+1
            tvar = strjoin(['thg',type,tsite],'_')
            tmp = scdfread(asifn, tvar, [rec0,rec1])
            asiimgs = *(tmp[0].value)
            
            ; count rate.
            crs = dblarr(nrec)
            
            ; loop through each time.
            for i = 0, n_elements(tuts)-1 do begin
                tmlon = scmlon[i]
                tmlat = scmlat[i]
                timg = (asiimgs[i,*,*])[*]
                idx = where(asimlon le tmlon+del and asimlon ge tmlon-del and $
                    asimlat le tmlat+del and asimlat ge tmlat-del, cnt)
                if cnt eq 0 then continue
                crs[where(asiuts eq tuts[i])] = mean(timg[idx])
            endfor
            idx = where(crs eq 0, cnt)
            if cnt ne 0 then crs[idx] = !values.d_nan
            if cnt eq nrec then continue
            tvar = 'thg_'+tsite+'_cr_rbsp'+tprobe
            tmp = snum2str(del)+'!9'+string("260b)+'!X'
            store_data, tvar, asiuts, crs/1000d, limits = $
                {ytitle:'thg '+tsite+'!Ccount k#', labels:tmp+'x'+tmp+' avg'}
        endforeach
    endforeach
end

; load e/b field, pos, fpoint, etc.
fn = sdiskdir('Works')+'/data/rbsp_thg/rbsp_efw_fld_2013_0501_04.tplot'
tplot_restore, filename = fn

;; load pflux.
;fn = sdiskdir('Works')+'/works/wygant/pflux_thm_asi/'+ $
;    '/rbspb_2013_0501_pflux_scott/RBSPb_mapped_pflux_2013-05-01_11-to-55-sec.dat'
;restore, filename = fn
;store_data, 'rbspb_pf1', etimes, s_para_mapped
;tplot, ['rbspb_pf1','thg_'+['atha','fsmi','gill','tpas']+'_cr']

etr = stoepoch(['2013-05-01/05:00','2013-05-01/09:30'])
probes = ['b']
sites = ['atha','tpas']
pfvscr, etr, sites, probes
end