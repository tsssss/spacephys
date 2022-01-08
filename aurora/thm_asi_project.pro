
; mlat/mlon are already none uniform, i.e., 
; they are the mlat/mlon grid, not pixel grid.


pro thm_asi_project, img0, lat0, lon0

    ; forward projection.
    
    ; interpolation.



end


; explore 1 site.

rootdir = spreproot()
et0 = stoepoch('2013-04-14 07:00:00')       ; wygant time.
et0 = stoepoch('2008-02-10 05:20:00')       ; thm crib time.
sites = 'fsmi'
type = 'ast'
type = 'asf'
asifn = sptn2fn_thm_asi(et0, sites, rootdir, type = type)
ascfn = sptn2fn_thm_asi(et0, sites, rootdir, type = 'asc')

cdfid = cdf_open(asifn)
vroot = strjoin(['thg',type,sites],'_')

; find record.
vname = vroot+'_epoch'
cdf_control, cdfid, variable = vname, get_var_info = tmp, /zvariable
cdf_varget, cdfid, vname, et, rec_count = tmp.maxrec+1, /zvariable
et = reform(et)
tmp = min(et-et0, rec, /absolute)

; raw image and epoch.
et0 = et[rec]
cdf_varget, cdfid, vroot, img, rec_start = rec, /zvariable       
cdf_close, cdfid
npx = (size(img, /dimensions))[0]   ; image size.

; elevation, mlat and mlon, glat and glon.
cdfid = cdf_open(ascfn)
cdf_varget, cdfid, vroot+'_elev', elev, /zvariable
cdf_varget, cdfid, vroot+'_mlat', mlat, /zvariable
cdf_varget, cdfid, vroot+'_mlon', mlon, /zvariable
cdf_varget, cdfid, vroot+'_glat', glat, /zvariable
cdf_varget, cdfid, vroot+'_glon', glon, /zvariable

stop


end