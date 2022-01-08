;+
; Combine Polar pre_or and def_or to sheng_or, use def_or if available
; use sdt data to fill in days of missing data.
; 
; This is a replacer for polar_gen_orb_cdf and polar_fill_orb_gap.
;-

pro polar_combine_pre_def_or_to_sheng_or


;---Constants.
    secofday = 86400d   ; sec.
    secofhour = 3600d   ; sec.

    
;---Settings.
    idir = sdiskdir('Research')+'/data/polar/orbit'
    odir = sdiskdir('Research')+'/sdata/polar/orb'
    logfn = shomedir()+'/polar_combine_pre_def_or_to_sheng_or.log'
    
    utr0 = time_double(['1996-02-27','2008-06-14'])
    
    dr0 = 60d   ; sec.
    ddr = 1d    ; error allowed.
    net0 = 1440d
    dnet = 0d

;---Preparation.
    nday = (utr0[1]-utr0[0])/secofday+1
    ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
    ut2s = ut1s+secofday

    for i=0, nday-1 do begin
        stryr = time_string(ut1s[i],tformat='YYYY')
        strut = time_string(ut1s[i],tformat='YYYYMMDD')
        ipath = idir+'/'+stryr
        ifn = 'po_or_[a-z]{3}_'+strut+'_v[0-9]{2}.cdf'
        ifn = sgetfile(ifn, ipath, /local_only)
        
        opath = odir+'/'+stryr
        ofn = 'po_or_sheng_'+strut+'_v01.cdf'
        
        cmd = ''
        if ifn eq '' then begin
            cmd += strut+': no data ...'
        endif else begin
            cdf = scdfread(ifn, 'Epoch')
            tets = *cdf[0].value
            tdr = sdatarate(tets)*1e-3    ; sec.
            tnet = n_elements(tets)
            if abs(tdr-dr0) gt ddr then $
                cmd += strut+': data rate '+sgnum2str(tdr)+' ...'
            if abs(tnet-net0) gt dnet then $
                cmd += strut+': # of record '+sgnum2str(tnet)+' ...'
        endelse
        if cmd ne '' then begin
            openw, lun, logfn, /append, /get_lun
            printf, lun, cmd
            free_lun, lun
            continue
        endif else printf, -1, strut+': copy from '+ifn+' to '+opath+'/'+ofn+' ...'
        
        if file_test(opath,/directory) eq 0 then file_mkdir, opath
        file_copy, ifn, opath+'/'+ofn, /overwrite
    endfor

end