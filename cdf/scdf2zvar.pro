;+
; for a given cdf, convert the rvars to zvars to save disk.
;-
pro scdf2zvar, ifn, ofn

    ; check ifn.
    if file_test(ifn) eq 0 then message, 'input file does not exist...'
    icdfid = cdf_open(ifn)
    cdfinq = cdf_inquire(icdfid)
    decoding = cdfinq.decoding
    encoding = cdfinq.encoding
    majority = cdfinq.majority
    extra = create_struct('clobber',1, majority,1, decoding,1, encoding, 1)
    
    ; create ofn.
    cdf_set_cdf27_backward_compatible, /yes
    ocdfid = cdf_create(ofn, _extra = extra)
    cdf_compression, ocdfid, set_compression = 5, set_gzip_level = 9
    cdfinq = cdf_inquire(ocdfid)
    nrvar = cdfinq.nvars
    nzvar = cdfinq.nzvars
    nvar = nrvar+nzvar
    if nvar eq 0 then begin
        ovarnames = ''
        ovarflags = -1
    endif else begin
        ovarnames = strarr(nvar)
        ovarflags = bytarr(nvar)
        for i = 0, nrvar-1 do begin
            varinq = cdf_varinq(ocdfid, i)
            ovarnames[i] = varinq.name
            ovarflags[i] = 0
        endfor
        for i = nrvar, nvar-1 do begin
            varinq = cdf_varinq(ocdfid, i-nrvar, /zvariable)
            ovarnames[i] = varinq.name
            ovarflags[i] = 1
        endfor
    endelse
    
    cdfinq = cdf_inquire(icdfid)
    cdf_control, icdfid, get_numattrs = natt
    
    ngatt = natt[0]
    nvatt = natt[1]
    nrvar = cdfinq.nvars
    nzvar = cdfinq.nzvars
    
    ; create all attributes.
    natt = cdfinq.natts
    attnames = strarr(natt)
    attflags = bytarr(natt)
    for i = 0, natt-1 do begin
        cdf_attinq, icdfid, i, attname, scope, dummy
        attnames[i] = attname
        attflags[i] = strmid(scope,0,1) eq 'G'
        gscope = attflags[i] & vscope = 1-gscope
        if cdf_attexists(ocdfid, attname) then cdf_attdelete, ocdfid, attname
        attid = cdf_attcreate(ocdfid, attname, $
            global_scope = gscope, variable_scope = vscope)
    endfor
    
    ; copy gatt.
    if ngatt gt 0 then begin
        idx = where(attflags eq 1b)     ; 1 for gatt.
        for i = 0, ngatt-1 do begin
            attname = attnames[idx[i]]
            cdf_control, icdfid, attribute = attname, get_attr_info = attinfo
            nentry = attinfo.numgentries
            maxentry = attinfo.maxgentry
            for j = 0, maxentry do begin
                if ~cdf_attexists(icdfid, attname, j) then continue
                cdf_attget, icdfid, attname, j, tval
                cdf_attput, ocdfid, attname, j, tval
            endfor
        endfor
    endif
    
    if nvatt eq 0 then vattnames = '' $
    else vattnames = attnames[where(attflags eq 0)]
    
    ; vars.
    nvar = nrvar+nzvar
    varnames = strarr(nvar)
    varflags = bytarr(nvar)
    for i = 0, nrvar-1 do begin
        varinq = cdf_varinq(icdfid, i)
        varnames[i] = varinq.name
        varflags[i] = 0
    endfor
    for i = nrvar, nvar-1 do begin
        varinq = cdf_varinq(icdfid, i-nrvar, /zvariable)
        varnames[i] = varinq.name
        varflags[i] = 1
    endfor
    
    ; copy all var.
    for j = 0, nvar-1 do begin
        i = varflags[j]? j-nrvar: j
        varinq = cdf_varinq(icdfid, i, zvariable = varflags[j])
        cdf_control, icdfid, variable = i, get_var_info = varinfo, $
            zvariable = varflags[j], get_zmode = zmode, get_cachesize = cache
        varname = varinq.name
        cdftype = varinq.datatype
        nelem = varinq.numelem
        recvary = 'rec_'+varinq.recvar
        maxrec = varinfo.maxrec+1
        dims = varflags[j]? varinq.dim: cdfinq.dim
        dimvary = varinq.dimvar
        extra = create_struct(cdftype,1, recvary,1, $
            'zvariable',1, 'numelem',nelem); treat dim shrink.
        idx = where(dimvary eq 1, cnt)
        if cnt eq 0 then begin      ; scalar.
            dims = 0 & dimvary1 = 0
        endif else begin
            dims = dims[idx] & dimvary1 = dimvary[idx]
        endelse
        if dims[0] ne 0 then extra = create_struct(extra, 'dimensions',dims)
        idx = where(ovarnames eq varname, cnt)
        if cnt ne 0 then cdf_vardelete, ocdfid, varname
        varid = cdf_varcreate(ocdfid, varname, dimvary1, _extra = extra)
;        ; settings.
        cdf_compression, icdfid, variable = varname, $
            get_var_compression = varcomp, get_var_gzip_level = vargzip
        cdf_compression, ocdfid, variable = varname, $
            set_var_compression = varcomp, set_var_gzip_level = vargzip
        ; copy data.
        cdf_varget, icdfid, varname, tval, /string, rec_start = 0
        dims = (dims[0] ne 0)? [maxrec,dims]: [maxrec]
        vals = make_array(type = size(tval,/type), dims)
        for i = 0, maxrec-1 do begin
            cdf_varget, icdfid, varname, tval, /string, rec_start = i
            vals[i,*,*,*,*,*,*,*] = srmdim(tval, dimvary)
        endfor
        vals = transpose(vals, shift(indgen(n_elements(dims)),-1))
        cdf_varput, ocdfid, varname, vals
        cdf_varrename, ocdfid, varname, strtrim(varname,2)
        ; copy vatt.
        for i = 0, nvatt-1 do begin
            if ~cdf_attexists(icdfid, vattnames[i], varname) then continue
            cdf_attget, icdfid, vattnames[i], varname, tval
            if size(tval,/type) eq 7 then $
                if strtrim(tval[0],2) eq '' then continue   ; empty tval.
            cdf_attput, ocdfid, vattnames[i], varname, tval ; no need zvar.
        endfor
    endfor
    
;    for j = 0, nrvar-1 do begin
;        varinq = cdf_varinq(icdfid, j)
;        cdf_control, icdfid, variable = j, get_var_info = varinfo
;        varname = varinq.name
;        cdftype = varinq.datatype
;        nelem = varinq.numelem
;        recvary = 'rec_'+varinq.recvar
;        maxrec = varinfo.maxrec+1
;        dims = cdfinq.dim
;        dimvary = varinq.dimvar
;;if strtrim(varname,2) ne 'INT_IMAGE' then continue
;        extra = create_struct(cdftype,1, recvary,1, $
;            'zvariable',0, 'numelem',nelem); treat dim shrink.
;        dimvary1 = dimvary
;        if dims[0] ne 0 then extra = create_struct(extra, 'dimensions',dims)
;;        cdf_vardelete, ocdfid, varname
;cdf_compression, icdfid, variable = varname, get_var_compression = varcomp, get_var_gzip_level = vargzip
;print, varcomp, vargzip
;        varid = cdf_varcreate(ocdfid, varname, dimvary1, _extra = extra)
;        cdf_compression, ocdfid, variable = varname, set_var_compression = varcomp, set_var_gzip_level = vargzip
;        cdf_varget, icdfid, varname, tval, /string, rec_start = 0, rec_count = maxrec
;        cdf_varput, ocdfid, varname, tval, rec_start = 0
;    endfor
    
    cdf_close, icdfid
    cdf_close, ocdfid
    
end

ifn = shomedir()+'/Downloads/po_level1_uvi_19970501_v01.cdf'
ofn = shomedir()+'/Downloads/po_l1_uvi_19970501_v01.cdf'

ifn = shomedir()+'/Downloads/po_k0_uvi_19970509_v01.cdf'
ofn = shomedir()+'/Downloads/po_k1_uvi_19970509_v01.cdf'

;ifn = shomedir()+'/Downloads/fa_k0_ies_19970501_v01.cdf'
;ofn = shomedir()+'/Downloads/fa_k1_ies_19970501_v01.cdf'

;ifn = shomedir()+'/Downloads/po_k0_uvi_19970501_v01.cdf'
;ofn = shomedir()+'/Downloads/po_k1_uvi_19970501_v01.cdf'
scdf2zvar, ifn, ofn
end