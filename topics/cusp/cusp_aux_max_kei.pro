;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro cusp_aux_max_kei

    lun = -1
    
    ids = cusp_id('default')

    rootdir = shomedir()+'/Google Drive/works/data/cusp'
    
    fns = rootdir+'/'+ids+'_all_data.tplot'
    fns = file_search(fns, count = nfn)
    if nfn eq 0 then message, 'no data ...'

    theids = []
    maxkeis = []
    
    for i = 0, nfn-1 do begin
        print, fns[i]
        tplot_restore, filename = fns[i]
        get_data, 'scidat', tmp, info
        
        potr = [info.polar.cusp.entry.ut,info.polar.cusp.exit.ut]
        get_data, 'po_ion_keflux_map', t0, dat
        idx = where(t0 ge potr[0] and t0 le potr[1])
        
        if info.polar.cusp.entry.ilat lt 0 then dat*= -1
        maxkei = abs(min(dat[idx])) ; negative is upward.
        if maxkei lt 1 then begin
            theids = [theids, info.id]
            maxkeis = [maxkeis, string(maxkei)]
        endif
    endfor
    
    for i = 0, n_elements(theids)-1 do printf, lun, theids[i], string(maxkeis[i])

end
