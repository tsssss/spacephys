;+
; Type:
; 	function.
; 
; Name:
; 	swvdas.
; 
; Purpose:
; 	Perform detrend and smooth.
;
;	For normal usage, specify suftr (super filter).
;	If f0 is in [n,m], each m component is filtered.
;
;	For test usage, do not define suftr in pace with swvmat, swvcwt,
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
;	suftr, in, type = dblarr[n], required.
;		"Super" filter. Set it for normal usage.
;		"Super" means it overwrites all other keywords.
;		It refers to number of record.
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
; Return:
; 	return, out, dblarr[n] or dblarr[n,m].
; 
; Example:
; 	t1 = swvdas(f0, suftr).
; 
; Notes:
; 	* The spectrogram is WRONG and only for test purpose.
;		Only usage is f1 = swvdas(f0, suftr).
;		Other usages are strongly discouraged.
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

function swvdas, f0, suftr, scale = scl, filter = ftr, inverse = inverse

	compile_opt idl2

	; check f0.
	sz = size(f0)
	ndim = sz[0]
	nrec = sz[1]

	; normal usage.
	if n_elements(suftr) ne 0 then begin
		p1 = suftr[0]
		p2 = suftr[1]
		if ndim eq 2 then begin
			p1 = [p1,1]
			p2 = [p2,1]
		endif
		return, smooth(f0, p1, /nan, /edge_wrap)-$
			smooth(f0, p2, /nan, /edge_wrap)
	endif

	; check scale.
	if n_elements(scl) eq 0 then begin
		nscl = 50
		q = (0.125D*nrec)^(1D/(nscl-1))
		scl = dblarr(nscl)
		scl[0] = 4
		for i = 0, nscl-2 do scl[i+1] = scl[i]*q
	endif
	scl = scl[uniq(round(scl))]
	nscl = n_elements(scl)

	if ndim eq 1 then begin
		; do transformation.
		dpre = smooth(f0, scl[0], /nan, /edge_wrap)
		t0 = dblarr(nrec, nscl)
		t0[*,0] = f0-dpre
		for i = 1, nscl-1 do begin
			dnow = smooth(f0, scl[i], /nan, /edge_wrap)
			t0[*,i] = dpre-dnow
			dpre = dnow
		endfor
		if ~keyword_set(inverse) then return, t0
	endif

	; check ftr.
	if n_elements(ftr) eq 0 then ftr = scl[[0,nscl-1]]
	ftr = round(ftr)

	; inverse transforma or filter.
	idx = where(scl ge ftr[0] and scl le ftr[1])
	if idx[0] eq -1 then message, 'wrong filter ...'
	return, total(f0[*,idx], 2)
		
end
