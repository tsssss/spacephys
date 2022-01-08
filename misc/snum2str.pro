;+
; Type: function.
; Purpose: Convert number to string, no white space on both sides.
;	For float/double, after period, no ending zeros, keep n digits.
; Parameters:
; 	nums, in, int/long/float/double, req. The given number.
;	ndigit, in, int, opt. Specify number of digits after period, works only if
;		the given number is a float/double. Default is 2.
; Keywords:
;   maxdigit, in, opt. Set max # of digits for ineger.
;   shortfloat, in, opt. Set to remove trailing .0s for float.
; Return: string. String of the number.
; Notes: Do not deal with complex(6) and dcomplex(9), maybe later.
;   If scientific notation is closer to the real vale, then use it.
; Dependence: none.
; History:
; 	2013-03-27, Sheng Tian, create.
; 	2014-11-09, Sheng Tian, add scientific notation function.
;-

function snum2str, nums0, ndigit, maxdigit = mdigit, shortfloat = short

	compile_opt idl2
	ex0 = !except
	
	if n_params() lt 1 then message, 'no input ...'
	
	if size(nums0,/type) ne 10 then nums = nums0 else begin
	repeat nums = *nums0 until size(nums,/type) ne 10
	endelse
	if n_elements(ndigit) eq 0 then ndigit = 2
	if n_elements(mdigit) eq 0 then mdigit = 6     ; '-1.2e3', default scientific notation will always work.
	
	tcode = size(nums,/type)           ; type code.
	if tcode eq 0 then return, ''      ; undefined.
	if tcode eq 7 then return, nums    ; string.
	if tcode eq 8 then message, 'invalid input: structure ...'
	if tcode eq 11 then message, 'invalid input: object ...'
	if tcode eq 6 or tcode eq 9 then $
	    message, 'do not deal with complex ...'
	
	; left with byte(1), int(2), long(3), float(4), double(5), 
	; uint(12), ulong(13), long64(14), ulong64(15).
	nnum = n_elements(nums)
	strs = strarr(nnum)
	isint = ~(tcode eq 4 or tcode eq 5)    ; 4: float, 5: double
	
	for i = 0, nnum-1 do begin
	   num0 = floor(nums[i])
	   strs[i] = strtrim(string(num0),2)
	   b = strlen(strs[i])
	   
	   ; scientific notation.
	   e = -ceil(alog10(1d/abs(nums[i])))  ; exponent.
	   tnum = nums[i]*10d^(-e)             ; digit term, abs value in [1,10].
	   
	   num0 = floor(tnum)
	   b = strlen(strtrim(string(num0),2)) ; sign included.
	   format = '(f'+strtrim(string(b+ndigit+1),2)+'.'+ $
	       strtrim(string(ndigit),2)+')'
	   str1 = string(tnum,format=format)
	   tmp = strsplit(str1,'0',/extract)           ; remove trailing 0s.
	   if n_elements(tmp) eq 1 then str1 = tmp[0]
	   if strmid(str1,0,1,/reverse_offset) eq '.' then $
	       str1 = strmid(str1,0,strlen(str1)-1)    ; remove trailing ".".
	   str1 += 'e'+strtrim(string(e),2)
	   if tnum eq 0 then str1 = '0'
	   
	   if keyword_set(exp) then begin
	       strs[i] = str1 & continue
	   endif
	   
	   b = strlen(strs[i])     ; sign included.
	   if isint then begin     ; for integer, use scientific notation >1e6 and <-1e5.
	       if b gt mdigit then strs[i] = str1 & continue
	   endif
	   
	   format = '(F'+strtrim(string(b+ndigit+1), 2)+'.'+$
	       strtrim(string(ndigit), 2)+')'
	   str2 = string(nums[i], format = format)
	   if ndigit eq 0 then str2 = strmid(str2,0,b) $
	   else if keyword_set(short) then begin
	       while strmid(str2,0,1,/reverse_offset) eq '0' do $
	           str2 = strmid(str2,0,strlen(str2)-1)
	       if strmid(str2,0,1,/reverse_offset) eq '.' then $
	           str2 = strmid(str2,0,strlen(str2)-1)
	   endif
	   
	   del1 = abs(double(str1)-nums[i])
	   del2 = abs(double(str2)-nums[i])
	   strs[i] = (del1 lt del2)? str1: str2
	endfor

    if nnum eq 1 then return, strs[0] else return, strs
end

print, snum2str(1200000)
print, snum2str(120000)
print, snum2str(-120000)
print, snum2str(-12000)
print, snum2str(1.23)
print, snum2str(1.23,/short)
print, snum2str(0.001)
print, snum2str(120000.0)
print, snum2str(120000.0,/short)
end
