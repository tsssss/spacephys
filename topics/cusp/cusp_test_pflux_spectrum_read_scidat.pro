;+
; read key info out of xxx_all_data.tplot.
; only for 1 s/c, either polar or fast.
;-
function cusp_test_pflux_spectrum_read_scidat, id, forpolar = forpolar


    rad = !dpi/180
    deg = 180/!dpi
    re = 6378d
    maxnfilt = 16
    
    info0 = { $
        name:'', $      ; 'polar' or 'fast'.
        id:'', $        ; 'yyyy_mmdd_hh'.
        dis:0d, $       ; average distance, in Re.
        vsc:0d, $       ; average s/c velocity, in km/s.
        vcv:0d, $       ; convection velocity.
        dilat:0d, $     ; cusp width in ilat in deg.
        pflux:0d, $     ; poynting flux calculated from original filters.
        pflux1:0d, $    ; poynting flux calculated from fac scale.
        pflux_hf:0d, $  ; poynting flux due to high freq noise.
        pflux_fac0:0d, $; poynting flux due to fac from original filters.
        hfreqlim:0d, $  ; high freq limit of the original filters, in deg.
        lfreqlim:0d, $  ; low freq limit of the original filters, in deg.
        faclim:0d, $    ; the scale of fac, in deg.
        ebratio:dblarr(maxnfilt), $
        espec:dblarr(maxnfilt), $
        bspec:dblarr(maxnfilt), $
        pfluxes:dblarr(maxnfilt), $
        filters:dblarr(maxnfilt), $
        nfilter:0d, $
        va:0d, $        ; local alfven speed, in km/s.
        kee:0d, $
        kei:0d, $
        tmp:0d}


; **** read basic info.
    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

    ; log info, search on conjunction list, then on polar list.
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = id)
    if size(loginfo,/type) ne 8 then begin
        logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
        loginfo = cusp_read_conjun_list(logfile, event = id, /nofast)
    endif

    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return, info0
    endif

    print, 'processing event: '+id
    store_data, '*', /delete

    ; info, id, pre0, [de,db,pf]idx, plottr, cusptr.
    info = forpolar? loginfo.polar: loginfo.fast
    pre0 = forpolar? 'po_': 'fa_'

    plottr = info.plot_time
    if plottr[1] lt plottr[0] then plottr[1]+= 86400d

    cusptr = info.cusp_time
    if cusptr[1] lt cusptr[0] then cusptr[1]+= 86400d

    ; filt0, filts, delims, facdelim, tscls.
    filt0 = info.filters    ; original filters.
    filts = filt0           ; filters, small > large.
    filts = filts[uniq(filts,sort(filts))]

    delims = info.noisedelim    ; below and above are noise.
    delims = minmax(delims)     ; sort.
    if delims[0] lt 0 then delims[0] = 0
    if delims[1] lt 0 then delims[1] = plottr[1]-plottr[0]

    facdelim = info.faclim      ; above is fac.
    if facdelim eq -1 then facdelim = max(filts)

    tmp = info.scaleinfo        ; preliminary, need datarate to modify.
    tscls = smkgmtrc(tmp[0],tmp[1],tmp[2],'n')


; **** read poynting flux info.
    datdir = rootdir+'/data/cusp'
    infofn = datdir+'/'+id+'_all_data.tplot'
    infovar = 'scidat'
    tplot_restore, filename = infofn
    
    get_data, infovar, scidatt0, scidat

    tmp = forpolar? 'POLAR': 'FAST'
    idx = where(tag_names(scidat) eq tmp)
    scinfo = scidat.(idx)

    dis = (scinfo.cusp.entry.dis+scinfo.cusp.exit.dis)*0.5
    dilat = abs(scinfo.cusp.entry.ilat)-abs(scinfo.cusp.exit.ilat)
    cuspdt = abs(scinfo.cusp.entry.ut-scinfo.cusp.exit.ut)
    vsc = -dis*re*dilat*rad/cuspdt ; positive if along anti-sunward.
    vcv = dis^1.5*0.625d           ; positive if along anti-sunward.

    nband = scinfo.nfilter
    pf0s = scinfo.sb.fs[0:nband-1]
    star = scinfo.pfstar[0:nband-1]     ; S -> S*.
    idx = where(star gt 1, cnt)
    star0 = (cnt eq 0)? 1: mean(star[idx])
    pf0s = pf0s*star0
    filts = scinfo.filters[0:nband-1]
    
    pf_lfreq = pf0s[nband-1]
    pf_hfreq = pf0s[0]
    pf_fac = pf0s[nband-2]
    pfs = pf0s[1:nband-3]
    freqs = filts[2:nband-2]


    tinfo = info0
    tinfo.name = forpolar? 'polar': 'fast'
    tinfo.id = id
    tinfo.dis = dis
    tinfo.vsc = vsc
    tinfo.vcv = vcv
    tinfo.dilat = abs(dilat)
    tinfo.pflux = total(pfs)
    tinfo.pflux_fac0 = pf_fac
    tinfo.pflux_hf = pf_hfreq
    tinfo.hfreqlim = filts[1]*abs(vsc)/dis/re*deg
    tinfo.lfreqlim = filts[nband-2]*abs(vsc)/dis/re*deg
    tinfo.faclim = abs(cuspdt*vsc/dis/re*deg)
    tinfo.va = scinfo.va
    tinfo.kei = scinfo.kei
    tinfo.kee = scinfo.kee

    ; apply common fac lim.
    pflux0 = pf0s[1:nband-3]
    filt0s = filts[2:nband-2]
    filt1s = [filt0s,cuspdt]
    filt1s = filt1s[where(filt1s le cuspdt)]
    if max(filt0s) gt cuspdt then begin
        pflux1 = interpol(pflux0, filt0s, filt1s)
    endif else begin
        filt1s = [filt0s[0:-2],cuspdt]
        pflux1 = interpol(pflux0, filt0s, filt1s)
    endelse
    tinfo.pflux1 = total(pflux1)
    
    tinfo.ebratio = scinfo.ebratio[1:nband-2,0]
    tinfo.espec = scinfo.ebratio[1:nband-2,1]
    tinfo.bspec = scinfo.ebratio[1:nband-2,2]
    tinfo.pfluxes = pflux1
    tinfo.filters = filt1s
    tinfo.nfilter = n_elements(filt1s)

    return, tinfo

end
