
pro cusp_check_pflux_limit, eventid, reload = reload, save_data = save_data, no_plot = no_plot, test = test
    
    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = eventid)
    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif
    
    
    ; directory to save output data and pdfs.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'

    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif
    
    if ~file_test(figdir,/directory) then file_mkdir, figdir
        
    

; **** basic info for both loading data and plotting.
    potr = loginfo.polar.plot_time
    fatr = loginfo.fast.plot_time
    potrcusp = loginfo.polar.cusp_time
    fatrcusp = loginfo.fast.cusp_time
    if potr[1] lt potr[0] then potr[1]+= 86400d
    if fatr[1] lt fatr[0] then fatr[1]+= 86400d
    if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
    if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
    id = loginfo.id

    ; filters, min to max, uniq.
    ; Will take the minimum and maximum filters as noise delimiters,
    ; if delimiters are excplitly set, then use them.
    ; Will take the second large filter as delimiter for fac and wave,
    ; if fac delimiter is set, then use it.

    ; filt0, original filter: f1,f2,f3,f4.
    ; We also want to include 0 and infinity (maxf).
    ;
    ; delims, rule out noise: n1, n2.
    ;
    ; Say n1=f2, n2=f4, the wanted filters are f2, f3, f4, ie, f1',f2',f3'.
    ; The filters used to filter mat spectrogram are: 0,f1',f2',f3',maxf.
    ; There are 5 filters, 4 bands, among the bands, 3 are wanted.
    ; High freq noise band is 0-f1', low freq noise band is f3'-maxf.
    ;
    ; faclim, sets where fac bands are.
    ; Say faclim=f2', then f1'-faclim are wave bands, faclim-f3' are fac bands.

    maxfilt = 1e31

    ; polar.
    pofilt0 = loginfo.polar.filters     ; original filters.
    pofaclim = loginfo.polar.faclim     ; above is fac.
    podelims = loginfo.polar.noisedelim ; below and above are noise.
    if (where(podelims lt 0))[0] eq -1 then podelims = minmax(podelims) ; sort if no -1.

    ; sort from small to large.
    pofilts = pofilt0
    pofilts = pofilts[sort(pofilts)]
    pofilts = pofilts[uniq(pofilts)]
    if pofilts[0] ne 0 then pofilts = [0,pofilts]   ; auto fill 0.

    ; default setting: [0,max original filters].
    if podelims[0] lt 0 then podelims[0] = min(pofilts)
    if podelims[1] lt 0 then podelims[1] = maxfilt

    ; get the wanted filters, add noise delimiters, and fac limit.
    idx = where(pofilts ge podelims[0] and pofilts le podelims[1], ponfilt)
    pofilts = pofilts[idx]
    if min(pofilts) gt podelims[0] then pofilts = [podelims[0],pofilts]
    if max(pofilts) lt podelims[1] then pofilts = [pofilts,podelims[1]]
    
    ; use faclim to omit unecessary filters.
    if pofaclim lt 0 then pofaclim = pofilts[n_elements(pofilts)-2]
    if pofaclim gt podelims[1] then pohasfac = 0 else begin
        pohasfac = 1
        pofacidx = where(pofilts lt pofaclim)
        pofilts = [pofilts[pofacidx],pofaclim,max(pofilts)]
    endelse

    ; includes wave bands and fac (if pohasfac=1).
    ponfilt = n_elements(pofilts)

    ; the filters used in filtering in mat spectrogram.
    pomatfilts = pofilts
    if min(pomatfilts) gt 0 then begin
        pomatfilts = [0,pomatfilts]
        pohashighnoise = 1
    endif else pohashighnoise = 0
    if max(pomatfilts) lt maxfilt then begin
        pomatfilts = [pomatfilts,maxfilt]
        pohaslownoise = 1
    endif else pohaslownoise = 0
    ponmatband = n_elements(pomatfilts)-1
    pomatbandids = string(indgen(ponmatband),format='(I0)')

    ; the wanted bands are within the wanted filters.
    ponband = ponfilt-1     ; includes wave bands and fac, separate later.
    pobandids = 'b'+string(indgen(ponband),format='(I0)')    ; in var name.
    pobandlabels = strarr(ponband)
    for i = 0, ponband-1 do pobandlabels[i] = $
        sgnum2str(pofilts[i],msgn=3)+'-'+sgnum2str(pofilts[i+1],msgn=3)+'s'


    ; fast.
    fafilt0 = loginfo.fast.filters     ; original filters.
    fafaclim = loginfo.fast.faclim     ; above is fac.
    fadelims = loginfo.fast.noisedelim ; below and above are noise.
    if (where(fadelims lt 0))[0] eq -1 then fadelims = minmax(fadelims) ; sort if no -1.

    ; sort from small to large.
    fafilts = fafilt0
    fafilts = fafilts[sort(fafilts)]
    fafilts = fafilts[uniq(fafilts)]
    if fafilts[0] ne 0 then fafilts = [0,fafilts]   ; auto fill 0.

    ; default setting: [0,max original filters].
    if fadelims[0] lt 0 then fadelims[0] = min(fafilts)
    if fadelims[1] lt 0 then fadelims[1] = maxfilt

    ; get the wanted filters, add noise delimiters, and fac limit.
    idx = where(fafilts ge fadelims[0] and fafilts le fadelims[1], fanfilt)
    fafilts = fafilts[idx]
    if min(fafilts) gt fadelims[0] then fafilts = [fadelims[0],fafilts]
    if max(fafilts) lt fadelims[1] then fafilts = [fafilts,fadelims[1]]
    
    ; use faclim to omit unecessary filters.
    if fafaclim lt 0 then fafaclim = fafilts[n_elements(fafilts)-2]
    if fafaclim gt fadelims[1] then fahasfac = 0 else begin
        fahasfac = 1
        fafacidx = where(fafilts lt fafaclim)
        fafilts = [fafilts[fafacidx],fafaclim,max(fafilts)]
    endelse

    ; includes wave bands and fac (if fahasfac=1).
    fanfilt = n_elements(fafilts)

    ; the filters used in filtering in mat spectrogram.
    famatfilts = fafilts
    if min(famatfilts) gt 0 then begin
        famatfilts = [0,famatfilts]
        fahashighnoise = 1
    endif else fahashighnoise = 0
    if max(famatfilts) lt maxfilt then begin
        famatfilts = [famatfilts,maxfilt]
        fahaslownoise = 1
    endif else fahaslownoise = 0
    fanmatband = n_elements(famatfilts)-1
    famatbandids = string(indgen(fanmatband),format='(I0)')

    ; the wanted bands are within the wanted filters.
    fanband = fanfilt-1     ; includes wave bands and fac, separate later.
    fabandids = 'b'+string(indgen(fanband),format='(I0)')    ; in var name.
    fabandlabels = strarr(fanband)
    for i = 0, fanband-1 do fabandlabels[i] = $
        sgnum2str(fafilts[i],msgn=3)+'-'+sgnum2str(fafilts[i+1],msgn=3)+'s'


    
    pot1 = potrcusp[1]-potrcusp[0]
    pot2 = max(loginfo.polar.filters)
    fat1 = fatrcusp[1]-fatrcusp[0]
    fat2 = max(loginfo.fast.filters)
    
    print, ''
    print, '**** event id: '+eventid
    
    if pot1 lt pot2 or fat1 lt fat2 then begin
        print, 'Polar'
        print, 'cusp time (sec): ', sgnum2str(pot1)
        print, 'filters (sec):   ', sgnum2str(pot2)
        
        print, 'Fast'
        print, 'cusp time (sec): ', sgnum2str(fat1)
        print, 'filters (sec):   ', sgnum2str(fat2)
    endif else print, 'No problem ...'
    
end


ids = cusp_id('south_imf')
idx = where(stregex(ids, '1999') eq 0, nid)
if nid eq 0 then message , 'no event found ...'
ids = ids[idx]

for i = 0, nid-1 do cusp_check_pflux_limit, ids[i]


end