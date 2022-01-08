;+
; Type:
;   function.
;
; Name:
;   sptrn2name.
;
; Purpose:
;   I do strput for date and time in given pattern.
;    
; Parameters:
;   pattern, in, type = string, required.
;     The pattern specify date and time info.
;     yyyy: year    mm: month
;     dd: day       hh: hour
;     mi: minute    ss: second
;     ms: millisec
;     
;   epoch, in, type = double, required.
;     Epoch contains the date and time info.
;       
; Keywords:
;   none.
;
; Return:
;   return, out, type = string.
;   
; Example:
;   fullname = spatterntofullname('~/data/fast/yyyy/yyyymmdd.cdf', epoch).
;   
; Notes:
;   * pattern is case insensitive.
;   * I do not deal with arrays for pattern and epoch.
;   * todo: support doy: dayofyear.
;   
; Dependence:
;   none.
;   
; Author:
;   Sheng Tian.
;
; History:
;   2011-07-26, Sheng Tian, create.
;   2012-07-12, Sheng Tian, revise.
;-

function sptrn2name, pattern, epoch
  
  compile_opt idl2
  
  on_error, 2
  
  ; I do not deal with array.
  fullname = pattern[0]
  
  ; break down the epoch.
  cdf_epoch, epoch, yr, mo, dy, hr, mi, sc, ms, /break
  
  ; replace the pattern.
  while (((ii = strpos(fullname, 'yyyy'))) ne -1) do $
    strput, fullname, string(yr, format = '(I4.4)'), ii
  while (((ii = strpos(fullname, 'mm'))) ne -1) do $ 
    strput, fullname, string(mo, format = '(I2.2)'), ii
  while (((ii = strpos(fullname, 'dd'))) ne -1) do $ 
    strput, fullname, string(dy, format = '(I2.2)'), ii
  while (((ii = strpos(fullname, 'hh'))) ne -1) do $ 
    strput, fullname, string(hr, format = '(I2.2)'), ii
  while (((ii = strpos(fullname, 'mi'))) ne -1) do $ 
    strput, fullname, string(mi, format = '(I2.2)'), ii
  while (((ii = strpos(fullname, 'ss'))) ne -1) do $ 
    strput, fullname, string(sc, format = '(I2.2)'), ii
  while (((ii = strpos(fullname, 'msc'))) ne -1) do $ 
    strput, fullname, string(ms, format = '(I3.3)'), ii
    
  return, fullname
    
end
