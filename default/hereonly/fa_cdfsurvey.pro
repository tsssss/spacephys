pro fa_cdfsurvey, orbs

	compile_opt idl2

	on_error, 2

	device, decomposed=0
	loadct3, 45
	
	norbs = n_elements(orbs)
	for ii=0,norbs-1 do begin
		orb = orbs[ii]
		cdf_eisurvey, orb
		@tplot_com
		t0 = tstring(time_offset)
		ob = string(orb,format='(I05)')
		fn = t0+'_'+ob+'_fa_cdfsvy_ele.ps'
		pstplot, filename = fn
		cdf_eisurvey, orb, /ion
		fn = t0+'_'+ob+'_fa_cdfsvy_ion.ps'
		pstplot, filename = fn
	endfor
end
		

functon tstring, t0
	ts = time_string(t0)
	yr = strmid(ts, 0, 4)
	mm = strmid(ts, 5, 2)
	dd = strmid(ts, 8, 2)
	return, yr+mm+dd
end
