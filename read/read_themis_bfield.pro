;+
; Read themis B field into tplot.
; Read FGS, which is at 3 sec cadence.
; deal with single probe: a-e.
; set newname to change the name of the tplot var.
; set varlist to add the produced tplot vars.
; set times to interpolate all vars to that time.
;-

pro read_themis_bfield, utr0, probe=probe, newname=newname, $
    addto=varlist, times=times

    compile_opt idl2

;---local and remote directory.
    sep = path_sep()
    if n_elements(locroot) eq 0 then locroot = spreproot('themis')
    if n_elements(remroot) eq 0 then $
        remroot = 'ftp://cdaweb.gsfc.nasa.gov/pub/data/themis'
    if n_elements(probe) eq 0 then message, 'Please provide a probe ...'

;---prepare file names.
    type1 = 'fgm'   ; use in filename.
    type2 = 'fgm'   ; use in path.
    vsn = (n_elements(version))? version: 'v[0-9]{2}'
    ext = 'cdf'
    level = 'l2'
    prb = probe[0]

    baseptns = ['th'+prb+'_'+level+'_'+type1,'_YYYYMMDD_',vsn+'.'+ext]
    nbaseptn = n_elements(baseptns)
    ptnflags = [0,1,0]
    rempaths = [remroot+'/th'+prb+'/'+level+'/'+type1,'YYYY',baseptns]
    locpaths = [locroot+'/th'+prb+'/'+level+'/'+type1,'YYYY',baseptns]
    ptnflags = [0,1,ptnflags]
    

;---load files.
    remfns = sprepfile(utr0, paths = rempaths, flags=ptnflags, nbase=nbaseptn)
    locfns = sprepfile(utr0, paths = locpaths, flags=ptnflags, nbase=nbaseptn)
    nfn = n_elements(locfns)
    for i = 0, nfn-1 do begin
        basefn = file_basename(locfns[i])
        locpath = file_dirname(locfns[i])
        rempath = file_dirname(remfns[i])
        locfns[i] = sgetfile(basefn, locpath, rempath)
    endfor
    idx = where(locfns ne '', nfn)    
    if nfn ne 0 then locfns = locfns[idx] else message, 'no file found ...'


;---prepare var names.
    pre0 = 'th'+prb+'_'
    var0s = pre0+['fgs_time','fgs_gsm']   ; 3 sec resolution.
    var1s = idl_validname(var0s)

;---module for variable loading.
    nvar = n_elements(var0s)
    if nvar ne n_elements(var1s) then message, 'mismatch var names ...'
    ptrs = ptrarr(nvar)
    ; first file.
    tmp = scdfread(locfns[0],var0s, skt=skt)
    for j = 0, nvar-1 do ptrs[j] = (tmp[j].value)
    ; rest files.
    for i = 1, nfn-1 do begin
        tmp = scdfread(locfns[i],var0s)
        for j = 0, nvar-1 do begin
            ; works for all dims b/c cdf records on the 1st dim of array.
            *ptrs[j] = [*ptrs[j],*(tmp[j].value)]
            ptr_free, tmp[j].value  ; release pointers.
        endfor
    endfor

    ; fill value.
    fillval = -1e31
    for j = 0, nvar-1 do begin
        idx = where((*ptrs[j]) eq fillval, cnt)
        if cnt eq 0 then continue
        (*ptrs[j])[idx] = !values.d_nan
    endfor



;---Move data to a structure then to tplot.
    ; move data to structure.
    dat = create_struct(var1s[0],*ptrs[0])
    for j = 1, nvar-1 do dat = create_struct(dat, var1s[j],*ptrs[j])
    for j = 0, nvar-1 do ptr_free, ptrs[j]
    
    ; make a uniform time.
    dr0 = 3d    ; sec.
    utr1 = utr0-(utr0 mod dr0)+[1,0]*dr0
    ut0s = smkarthm(utr1[0], utr1[1], dr0, 'dx')
    if keyword_set(times) then ut0s = times
    ; interpol to that uniform time.
    bgsm = sinterpol(dat.(1), dat.(0), ut0s)

    myvars = pre0+['b_gsm']
    if not keyword_set(newname) then newname = 'b_gsm'
    store_data, pre0+newname, ut0s, bgsm, $
        limits = {ytitle:'(mV/m)', colors:[6,4,2], labels:'GSM B'+['x','y','z'], labflag:-1}
    if not keyword_set(varlist) then varlist=[]
    varlist = [varlist, myvars]
end


vars = []
utr = time_double(['2014-08-28/09:00','2014-08-28/12:00'])
tprobe = 'e'
read_themis_bfield, utr, probe=tprobe, addto=vars

end

