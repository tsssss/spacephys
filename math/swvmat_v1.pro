;+
; Type:
; 	function.
; 
; Name:
; 	swvmat.
; 
; Purpose:
; 	Perform moving average transform, or inverse transform.
;
;	There are 3 behaviours:
;		(1) Do transform, return dblarr[n,m], when
;			f0 is dblarr[n], inverse is false;
;		(2) Do inverse transform, return dblarr[n], when
;			f0 is dblarr[n,m], inverse is true (optional);
;		(3) Do filter, return dblarr[n], when
;			f0 is dblarr[n], inverse is true.
;
;	Default scale is 50 log scales in [4,0.5*nrecs].
;	Default filter is [scale[0], scale[nscale-1]].
;	scale and filter are round to integer.
; 
; Parameters:
; 	f0, in, type = dblarr[n] or dblarr[n,m], required.
;		Input signal in dblarr[n] or 
;		transform func in dblarr[n,m].
;
; Keywords:
; 	scale = scl, in/out, type = dblarr[m], optional.
;		Set scales for transform or return the default scale.
;		Scale refers to number of record.
;
;	filter = ftr, in/out, type = dblarr[m], optional.
;		Set filter for inverse tranaform.
;		Filter refers to number of record.
;
;	inverse = inverse, in, type = boolean, optional.
;		Set inverse to do inverse transform or filter.
;
;	correct = cor, in/out, type = dblarr[m], optional.
;		Set correct to do amplitude correction.
;		cor[i] is the Fourier period of channel wvt[i,*].
;		After correction, cor is in dblarr[m,n], the corrected wvt.
; 
; Return:
; 	return, out, dblarr[n] or dblarr[n,m].
; 
; Example:
; 	t1 = swvmat(f0).
;	f1 = swvmat(f0, filter = recftr, /inverse)
;	f1 = swvmat(t1, filter = recftr, /inverse)
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
; 	2013-03-24, Sheng Tian, create.
;-

function swvmat_v1, f0, scale = scl, filter = ftr, inverse = inverse, $
	correct = cor

	compile_opt idl2

	; check f0.
	sz = size(f0)
	ndim = sz[0]
	nrec = sz[1]

	; check scale.
	if n_elements(scl) eq 0 then begin
		nscl = 50
		q = (0.125D*nrec)^(1D/(nscl-1))
		scl = dblarr(nscl)
		scl[0] = 4
		for i = 0, nscl-2 do scl[i+1] = scl[i]*q
	endif
	scl = round(scl)
	scl = double(scl[uniq(scl)])
	nscl = n_elements(scl)

	if ndim eq 1 then begin
		; do transformation.
		dflt = f0
		wvt = dblarr(nrec, nscl)
		for i = 0, nscl-1 do begin
			dnow = smooth(dflt, scl[i], /nan, /edge_wrap)
			wvt[*,i] = dflt-dnow
			dflt = dnow
            dnow = smooth(wvt[*,i], scl[i], /nan, /edge_wrap)
            wvt[*,i] -= dnow            ; purify.
            dflt += dnow
		endfor
		if ~keyword_set(inverse) then return, wvt
	endif

	; check ftr.
	if n_elements(ftr) eq 0 then ftr = scl[[0,nscl-1]]
	ftr = round(ftr)

	if ndim eq 2 then wvt = f0

	; do correction.
	if n_elements(cor) ne 0 then begin
		ws = cor
		ns = n_elements(ws)
		u = !dpi*(ws##(1D/ws))
		a = sin(u)/u
		cor = wvt
		; n-1 rounds of correction.
		for i = ns-1, 1, -1 do begin
			; coeffients related to the ith layer.
			c = dblarr(i+1)
			c[0] = 1
			for j = 1, i do c[j] = c[j-1]*a[i,j-1]
			cor[*,i-1] *= (1D/c[i-1])
			for j = i-2, 0, -1 do cor[*,j] -= cor[*,i-1]*(c[j]-c[j+1])
		endfor
		return, cor
	endif

;    if n_elements(cor) ne 0 then begin
;        ws = cor
;        ns = n_elements(ws)
;        u = !dpi*(ws##(1D/ws))
;        a = sin(u)/u
;        cor = wvt
;        ; n-1 rounds of correction.
;        for i = ns-1, 1, -1 do begin
;            ; coeffients related to the ith layer.
;            c = dblarr(i+1)
;            c[0] = 1
;            for j = 1, i do c[j] = c[j-1]*a[i,j-1]
;            cor[*,i] *= (1D/c[i])
;            for j = i-1, 0, -1 do cor[*,j] -= cor[*,i]*(c[j]-c[j+1])
;        endfor
;        return, cor
;    endif

	; inverse transform or filter.
	idx = where(scl ge ftr[0] and scl le ftr[1])
	if idx[0] eq -1 then message, 'wrong filter ...'
	if n_elements(idx) eq 1 then return, wvt[*,idx[0]]
	return, total(wvt[*,idx], 2)

end

