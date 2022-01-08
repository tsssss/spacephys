
; calc count around given footpoint track(s) within nxn deg box.
; etr must be in the same day.
; fpvars is a tplot var in [[mlon],[mlat]].
; sc: whose foopoint, a string that enters the label.

pro thm_asi_countrate, etr, sites, fpvars, del = del, sc = sc

    pfvar = ''
    posvar = ''
    type = 'asf'
    asiloc = spreproot('themis')
    asirem = 'http://themis.ssl.berkeley.edu/data/themis'
    utr = sfmepoch(etr,'unix')
    
    ; breakdown etr epoch.
    ethr = 3600000d
    etrs = [etr[0]-(etr[0] mod ethr),etr[1]-(etr[1] mod ethr)]
    if (etr[1] mod ethr) eq 0 then etrs[1]-= ethr
    etrs = smkarthm(etrs[0],etrs[1],ethr,'dx')
    
    ; how fussy, deg x deg square.
    if n_elements(del) eq 0 then del = 1d   ; deg.
    if n_elements(sc) eq 0 then sc = ''
    
    if n_elements(sites) eq 0 then $
        sites = ['atha','chbg','ekat','fsmi','fsim','fykn',$
            'gako','gbay','gill','inuv','kapu','kian',$
            'kuuj','mcgr','pgeo','pina','rank','snkq',$
            'tpas','whit','yknf','nrsq','snap','talo']
    nsite = n_elements(sites)
    
    if n_elements(fpvars) eq 0 then message, 'no footpoint variable ...'
    
    ; collect times for aurora image.
    asiuts = sfmepoch(smkarthm(etr[0]-(etr[0] mod 3000),etr[1]-(etr[1] mod 3000), 3000, 'dx'), 'unix')
;    flags = bytarr(nsite)+1b    ; 1: have data in etr.
;    asiuts = []
;    for i = 0, nsite-1 do begin
;        tsite = sites[i]
;        asifn = sptn2fn_thm_asi(etr[0], tsite, asiloc, type = type)
;        asifn = sprepfile(asifn, asiloc, asirem)
;        if asifn[0] eq '' then begin
;            flags[i] = 0b & continue
;        endif
;        tmp = scdfread(asifn, strjoin(['thg',type,tsite,'time'],'_'))
;        tuts = *(tmp[0].value)
;        idx = where(tuts ge utr[0] and tuts le utr[1], cnt)
;        if cnt eq 0 then begin
;            flags[i] = 0b & continue
;        endif
;        asiuts = [asiuts,tuts]
;        asiuts = asiuts[sort(asiuts)]
;        asiuts = asiuts[uniq(asiuts)]
;    endfor
    
    nrec = n_elements(asiuts)
    if nrec eq 0 then message, 'no site has data ...'
    
;    sites = sites[where(flags eq 1b)]
    nsite = n_elements(sites)
    print, 'asi available at ', sites
    
    foreach tfpvar, fpvars do begin
        ; track each site.
        foreach tsite, sites do begin
            ; read pixel mlat, mlon in deg.
            ascfn = sptn2fn_thm_asi(etrs[0], tsite, asiloc, type = 'asc')
            ascfn = sprepfile(paths = ascfn)
            tvar = strjoin(['thg',type,tsite,''],'_')+['mlat','mlon']
            tmp = scdfread(ascfn, tvar, 1)
            asimlon = *(tmp[1].value)
            asimlat = *(tmp[0].value)
            asimlon = reform(asimlon[1,*,*])
            asimlat = reform(asimlat[1,*,*])
            asimlon = asimlon[*]        ; make 1-d.
            asimlat = asimlat[*]

            ; read all epoch, find which are within trange.
            nfn = n_elements(etrs)
            asifns = strarr(nfn)
            for i = 0, nfn-1 do begin
                asifns[i] = sptn2fn_thm_asi(etrs[i], tsite, asiloc, type = type)
                asifns[i] = sprepfile(paths = asifns[i])
            endfor
            asifns = asifns[where(asifns ne '')]
            
            cruts = []
            crs = []
            foreach asifn, asifns do begin
                tmp = scdfread(asifn, strjoin(['thg',type,tsite],'_')+['','_time'])
                tuts = *(tmp[1].value)
                idx = where(tuts ge utr[0] and tuts le utr[1], cnt)
                tuts = tuts[idx]
                asiimgs = (*(tmp[0].value))[idx,*,*]
                
                ; read s/c footpoint.
                get_data, tfpvar, t0, dat
                scmlon = interpol(dat[*,0], t0, tuts, /quadratic)
                scmlat = interpol(dat[*,1], t0, tuts, /quadratic)
                
                ; count rate.
                tcrs = dblarr(n_elements(tuts))
                
                ; loop through each time.
                for i = 0, n_elements(tuts)-1 do begin
                    tmlon = scmlon[i]
                    tmlat = scmlat[i]
                    timg = (asiimgs[i,*,*])[*]
                    idx = where(asimlon le tmlon+del and asimlon ge tmlon-del and $
                        asimlat le tmlat+del and asimlat ge tmlat-del, cnt, complement=idx2)
                    if cnt eq 0 then continue
                    tcrs[i] = mean(timg[idx])
                    tcrs[i] = max(timg[idx])
                endfor
                cruts = [cruts, tuts]
                crs = [crs, tcrs]
            endforeach
            
            idx = where(crs eq 0, cnt)
            if cnt ne 0 then crs[idx] = !values.d_nan
            if cnt eq nrec then continue
            tvar = tfpvar+'_thg_'+tsite+'_cr'
            tmp = snum2str(del)+'!9'+string("260b)+'!X'
            store_data, tvar, cruts, crs/1000d, limits = $
                {ytitle:'thg '+tsite+'!Ccount k#', $
                labels:tmp+'x'+tmp+' max!C  site: '+strupcase(tsite)+'!C  '+sc}
        endforeach
    endforeach
end

;; load e/b field, pos, fpoint, etc.
;fn = sdiskdir('Research')+'/data/rbsp_thg/rbsp_efw_fld_2013_0501_04.tplot'
;tplot_restore, filename = fn
;
;; load pflux.
;fn = sdiskdir('Research')+'/works/wygant/pflux_thm_asi/'+ $
;    '/rbspb_2013_0501_pflux_scott/RBSPb_mapped_pflux_2013-05-01_11-to-55-sec.dat'
;restore, filename = fn
;store_data, 'rbspb_pf1', etimes, s_para_mapped

etr = stoepoch(['2013-05-01/06:00','2013-05-01/08:00'])
fpvars = 'rbsp'+['b']+'_fpt_lonlat'
thm_asi_countrate, etr, sites, fpvars

ylim, 'rbspb_fpt_lonlat_thg_*', 0,2e1, 0
tplot, ['rbspb_pf1','rbspb_fpt_lonlat_thg_'+['atha','fsmi','gill','tpas']+'_cr']
end