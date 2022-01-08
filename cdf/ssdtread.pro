;+
; Type: function.
; Purpose: I read all data in given sdt file.
; Parameters:
;   fn, in, string, req. Full file name of sdt file.
; Keywords:
;   none.
; Return: struct.
;     sdt = { $
;       __name: 'sdt', $
;       header: header, $
;       att: gatt, $
;       var: var $
;     }
;  
;     header = { $
;       __name: 'sdt.header', $
;       filename: fn, $
;       format: format, $
;       version: version }
;       
;     att = { $
;       natts: 1, $
;       comment: comment }
;       
;     var = { $
;       nvars: nvars, $           ; number of variable.
;       varname1: var1, $         ; the 1st variable.
;       varname2: var2, $         ; the 2nd varialbe.
;       ...
;       varnameN: varN }
;       
;     varN = { $
;       name: varname0, $
;       value: value, $
;       depend_0: time, $
;       nrecs: nrecs, $
;       dims: dims, $
;       att: att }
;     
;     att = {natts: 0}
;
; Dependence: none.
; Notes:
;   * scalar data are in [n].
;   * m-dimension vector data are in [m, n].
;   * I cannot deal with spectrum or image data.
; History:
;   2012-09-13, Sheng Tian, create.
;-
function ssdtread, fn
  compile_opt idl2
  
  ; check if file exists.
  if n_elements(fn) eq 0 then $
    message, 'file name is undefined ...'
  if file_test(fn) ne 1 then $
    message, 'file does not exist ...'
  
  ; open file.
  openr, lun, fn, /get_lun
  
  line = ''
  readf, lun, line
  
  ; ===========
  ; header info
  ; ===========
  header = strarr(7)
  while ~eof(lun) do begin
    if stregex(line,'^%*',/boolean) then break
    readf, lun, line
  endwhile
  
  ;  if reach end of file.
  if eof(lun) then $
    message, 'no data in file ...'
  
;  tab = string(9B)
  readf, lun, header
  dummy = ''
    
  ; format line.
  format = ''
  reads, header[1], dummy, format, format = '(A12, A)'
  
  ; version line.
  version = ''
  reads, header[3], dummy, version, format = '(A15, A)'
  
  ; comment line.
  comment = ''
  reads, header[4], dummy, comment, format = '(A12, A)'
  
  header = { $
    __name: 'sdt.header', $
    filename: fn, $
    format: format, $
    version: version }
  
  print, 'reading file ' + fn + ' ...'
  
  ; ==================
  ; global attributes.
  ; ==================
  gatt = { $
    natts: 1, $
    comment: comment }
  
  ; ==========
  ; variables.
  ; ==========
  nvars = 0
  
  ; read all variables.
  varheader = strarr(6)
  vartail = strarr(2)
  var = { $
    nvars: nvars }
  
  readf, lun, line
  while ~eof(lun) do begin
    if ~stregex(line,'^%Start',/boolean) then begin
      readf, lun, line
      continue
    endif
    
    if eof(lun) then $
      message, 'reach end of file ...'
    
    ; read variable header.
    readf, lun, varheader
    
    ; variable name.
    varname0 = ''
    reads, line, dummy, varname0, format = '(A7, A)'
    
    ; number of data records.
    nrecs = 0UL
    reads, varheader[3], dummy, nrecs, format = '(A8, I)'
    
    ; read the data.
    if nrecs gt 0 then begin
        data = fltarr(3, nrecs)
        readf, lun, data
    endif else data = !values.d_nan
    
    ; read variable tail.
    readf, lun, vartail
    
    ; save the variable.
    date = 0L
    reads, varheader[2], dummy, date, format = '(A5, I)'
    caldat, date, mo, dy, yr
    cdf_epoch, epoch, yr, mo, dy, /compute
    epoch0 = 62167219200000D   ; 1970-01-01 00:00:00.000 UTC.
    time0 = (epoch - epoch0) /1000D
    
    if nrecs gt 0 then begin
        time = transpose(data[1,*]) + time0   ; unix time.
        value = transpose(data[2,*])
    endif else begin
        time = time0
        value = !values.d_nan
    endelse
    size = size(value)
    ndims = size[0] - 1
    dims = (ndims gt 0) ? size[1: 1+ndims-1]: 0
    
    vatt = {natts: 0}

    thisvar = { $
      name: varname0, $
      value: value, $
      depend_0: time, $
      nrecs: nrecs, $
      dims: dims, $
      att: vatt }
    
    varname = idl_validname(varname0, /convert_all)
    if ~strcmp(varname, varname0) then $
      message, 'convert ' + varname0 + ' to ' + varname + ' ...', /continue
      
    var = create_struct(var, varname, thisvar)
    var.nvars = var.nvars + 1
    
    ; continue to read following vars.
    if ~eof(lun) then $
      readf, lun, line
      
  endwhile
  
  free_lun, lun
  
  ssdt = { $
    __name: 'sdt', $
    header: header, $
    att: gatt, $
    var: var $
  }
  
  return, ssdt

end
