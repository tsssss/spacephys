;+
; Type:
; 	function.
; 
; Name:
; 	swvmatcorr.
; 
; Purpose:
; 	Do amplitude correction for moving average transform.
; 
; Parameters:
; 	mat, in, type = dblarr[n,m], required.
;		Moving average transform's core function.
;
;	scales, in/out, type = dblarr[m], required.
;		Scales in Fourier period.
; 
;	width0, in, type = double, required.
;		Widths on which the correction is done.
;
; Keywords:
; 	none.
; 
; Return:
; 	return, out, type = dblarr[n,m].
;		The fixed moving average transform core.
; 
; Example:
; 	mat = swvmatcorr(mat, recscl, width).
; 
; Notes:
; 	none.
; 
; Dependence:
; 	none.
; 
; Author:
; 	Sheng Tian.
; 
; History:
; 	2013-04-01, Sheng Tian, create.
;-

function swvmatcorr, mat0, scales, width0

	compile_opt idl2

	del = min(scales-width0, ns, /absolute)
	if del lt 0 then ns--      ; choos the closest scale < width.
	ws = scales[0:ns]
	u = (!dpi/width0)*ws
	a = sin(u)/u
	mat = mat0
	; coefficients related to layer n.
	c = dblarr(ns+1)
	nc = n_elements(c)
	c[0] = 1
	for j = 1, nc-1 do c[j] = c[j-1]*a[j-1]
	wv0 = mat[*,ns]*(1D/c[nc-1])
	for j = ns-2, 0, -1 do mat[*,j] -= mat[*,ns-1]*(c[j]-c[j+1])
	mat[*,ns] = wv0
	scales[ns] = width0

	return, mat

end
