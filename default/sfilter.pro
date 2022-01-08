;+
; Type:
;	function.
;
; Name:
;   sfilter.
;
; Purpose:
;   I filter data in given range.
;
; Parameters:
;   data0, in, type = array of certain type, required.
;     can be scalar array in [n].
;     can be m-dimension vector array in [m, n].
;     can be m-dimension vector array in [n, m].
;   
;   p10, in, type = double, required.
;     lower limit of period.
;   
;   p20, in, type = double, required.
;     upper limit of period.
;     
;   datarate, in, type = double, optional.
;     If set data rate, then I assume p10 and p20 are in sec,
;     otherwise, p10 and p20 are in # of record.
;   
; Keywords:
;   _extra = extra, in, optional.
;     Keywords for smooth() in IDL library.
;   
; Return:
;   return, out, type = array of certain type.
;     The filtered data, has same dimension as data0.
;   
; Example:
;   e_filter = sfilter(efield, 0.5, 180, datarate).
;   
; Notes:
;   * When data0 is an array of m-dimension vector,
;     I automatically tell it is in [m, n] or [n, m],
;     based on the assumption that n > m (which is true
;     in most case). However, dangers DOES exist!
;     
;   * The filtering is based on smooth(), it is NOT
;     at all a filter in frequency space.
;
;	* The NaN keyword of smooth() is set by default.
;  
; Dependence:
;   none.
;   
; Author:
;   Sheng Tian
; 
; History:
;   2012-09-18, Sheng Tian, create.
;-

function sfilter, data0, p10, p20, datarate, $
  smooth = smooth, detrnd = detrnd, _extra = extra

  ; p10 and p20 are in sec.
  if n_elements(datarate) ne 0 then begin
    p1 = p10 / datarate
    p2 = p20 / datarate
  ; p10 and p20 are in # of records.
  endif else begin
    p1 = p10
    p2 = p20
  endelse
  
  if p1 eq p2 then $
    message, 'p1 equal p2 ...'
  if p1 gt p2 then begin  ; ensure p1 < p2.
    tmp = p1
    p1 = p2
    p2 = tmp
  endif
  
  thesize = size(data0)
  ndims = thesize[0]
  
  if ndims gt 2 then $
    message, '# of dimension greater than 2 ...'
    
  if ndims eq 0 then $
    message, 'data is scalar ...'
    
  if ndims eq 2 then begin
    ; data in [n, m].
    if thesize[1] gt thesize[2] then begin
      p1 = [p1,1]
      p2 = [p2,1]
    ; data in [m, n].
    endif else begin
      p1 = [1,p1]
      p2 = [1,p2]
    endelse
  endif
  
  if n_elements(smooth) ne 0 then $
    return, smooth(data0, p1, /nan, _extra = extra)
  if n_elements(detrnd) ne 0 then $
    return, data-smooth(data0, p1, /nan, _extra = extra)
  ; p1 < p2.  
  return, smooth(data0, p1, /nan, _extra = extra) - $
    smooth(data0, p2, /nan, _extra = extra)
  
end
