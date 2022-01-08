;+
; Type: function.
; Purpose: Convert pattern to filename.
; Parameters:
; 	ptn, in, string, req. The pattern specifies date and time info.
;		yyyy: year, MM: month, dd: day, doy: day of year, hh: hour, mm: minute, 
;		ss: second, lll: milli sec, ccc: micro sec, ppp: pico sec.
;	et0, in, double, req. Epoch of the day.
; Keywords: none.
; Return: string. The file name.
; Notes:
;       Pattern is case-sensitive for MM and mm, otherwise case insensitive.
;   The returned full name does NOT necessarily exist.
; Dependence: stodoy.
; History:
; 	2011-07-26, Sheng Tian, create.
;	2012-07-12, Sheng Tian, revise.
;-

function sptn2fn, ptn, et0
	compile_opt idl2

	fn = ptn[0]
	et = et0[0]
	if size(et,/type) eq size(0d,/type) then begin
	   cdf_epoch, et, yr, mo, dy, hr, mi, sc, ms, /breakdown_epoch
	   mc = 0 & pc = 0
	endif else $
	   cdf_epoch16, et, yr, mo, dy, hr, mi, sc, ms, mc, pc, /breakdown_epoch

    ; case sensitive part.
	while (i = strpos(fn,'MM')) ne -1 do $
		strput, fn, string(mo, format = '(I02)'), i
	while (i = strpos(fn,'mm')) ne -1 do $
		    strput, fn, string(mi, format = '(I02)'), i

    ; case insensitive part.
	while (i = stregex(fn,'yy+',length=len,/fold_case)) ne -1 do $
	   strput, fn, strmid(string(yr,format='(I04)'),4-len,len), i
	
	while (i = stregex(fn,'dd',/fold_case)) ne -1 do $
		strput, fn, string(dy, format = '(I02)'), i
	while (i = stregex(fn,'hh',/fold_case)) ne -1 do $
		strput, fn, string(hr, format = '(I02)'), i
	while (i = stregex(fn,'ss',/fold_case)) ne -1 do $
		strput, fn, string(sc, format = '(I02)'), i
	while (i = stregex(fn, 'doy',/fold_case)) ne -1 do $
		strput, fn, string(stodoy(yr,mo,dy), format = '(I03)'), i
	while (i = stregex(fn, 'lll',/fold_case)) ne -1 do $
		strput, fn, string(ms, format = '(I03)'), i
    while (i = stregex(fn, 'ccc',/fold_case)) ne -1 do $
		strput, fn, string(mc, format = '(I03)'), i
	while (i = stregex(fn, 'ppp',/fold_case)) ne -1 do $
		strput, fn, string(pc, format = '(I03)'), i

	return, fn

end
