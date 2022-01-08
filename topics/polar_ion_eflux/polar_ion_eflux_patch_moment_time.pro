;+
; Patch on polar_ion_eflux_calc_and_save_moments.
; add corrected uts to cdf files.
;-
pro polar_ion_eflux_patch_moment_time, utr0

    fn = sdiskdir('Research')+'/sdata/polar/timas/'+$
        time_string(utr0[0],tformat='YYYY')+'/po_tim_moments_'+$
        time_string(utr0[0],tformat='YYYY_MMDD')+'.cdf'
    if file_test(fn) eq 0 then return
    tms = scdfread(fn, skt=skt)
    vnames = tms.name
    idx = where(vnames eq 'utsec', cnt)
    if cnt ne 0 then return
    
    idx = where(vnames eq 'ut_sec')
    ut0s = *tms[idx].value
    ut1s = polar_ion_eflux_fix_timas_time(ut0s, error=err)
    ainfo = {$
        FIELDNAM:'Center UT for the moments, corrected by polar_ion_eflux_fix_timas_time',$
        UNITS:'s (UT)',$
        VAR_TYPE:'support_data'}
    scdfwrite, fn, 'utsec', value=ut1s, cdftype='CDF_DOUBLE', attribute=ainfo

end

secofday = 86400d

utr0 = time_double(['1996-03-17','1997-12-31']) 
;utr0 = time_double(['1998-01-01','1998-12-08'])
nday = (utr0[1]-utr0[0])/secofday
ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
ut2s = ut1s+secofday

stop
for i=0, nday-1 do begin
    tutr = [ut1s[i],ut2s[i]]
    polar_ion_eflux_patch_moment_time, tutr
endfor

end
