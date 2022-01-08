;+
; Conclusion: use interpolation.
; There are 3 versions in the comparison: v1 is interpolated;
; v2 is smoothed upto 0.5 of Nyquist period of dB; v3 is interpolated
; version of v2, i.e, smooth first then interpolate.
; v2 is of course better than v1, but v3 looks attenuated slightly.
; Choose v1 because v1 is pretty similar to v3, and v1 is simpler to get.
;-
pro test_efield_interpol_or_smooth

    ; Problem:
    ;     We want to lower the sample rate of E field.
    ; For FAST, data rate is 0.125s for B field, is 0.03125s for E field.
    ; In calculating Poynting flux, there is no need to keep the high
    ; frequency signals that E has but B doesn't.
    ; 
    ;     One way is to resample E using interpolation.
    ; The other way is to use smoothing to remove high frequency part in E.
    
    fn = sdiskdir('SHENG')+'/works/polarcap/data/'+$
        'fa_sdt_fld_19980925_08278.tplot'
    if file_test(fn) eq 0 then message, 'no such file ...'
    tplot_restore, filename = fn
    
    ; get t0 from dB_fac_v.
    maxnrec = 100000L   ; prevent memory overflow.
    get_data, 'dB_fac_v', t0, dbfac
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    t0 = smkarthm(t0[0], t0[nrec-1], dr, 'dx')
    nrec = n_elements(t0)
    print, 'dB_fac_v data rate: ', dr, 's'
    
    ; original field is E_ALONG_V.
    get_data, 'E_ALONG_V', tmp, dev0
    options, 'E_ALONG_V', 'labels', 'dEv v0!C  original'
    dr1 = sdatarate(tmp)
    print, 'dEv data rate: ', dr1, 's'
    print, floor(dr/dr1)
    ; version 1: int.
    dev_interp = sinterpol(dev0, tmp, t0)
    store_data, 'dev_interp', data = {x:t0, y:dev_interp}, $
        limits = {ytitle:'dEv!C(mV/m)', labels:'dEv v1!C  direct interp'}
    ; version 2: sm0.
    dev_smooth0 = smooth(dev0, 0.5*floor(dr/dr1), /nan, /edge_wrap)
    store_data, 'dev_smooth0', data = {x:tmp, y:dev_smooth0}, $
        limits = {ytitle:'dEv!C(mV/m)', labels:'dEv!C  smooth only'}
    ; version 3: sm1.
    dev_smooth1 = sinterpol(dev_smooth0, tmp, t0)
    store_data, 'dev_smooth1', data = {x:t0, y:dev_smooth1}, $
        limits = {ytitle:'dEv!C(mV/m)', labels:'dEv v2!C  smooth&interp'}
    dev_del = dev_smooth1-dev_interp
    store_data, 'dev_del', data = {x:t0, y:dev_del}, $
        limits = {ytitle:'dEv!C(mV/m)', labels:'dEv v1-v2'}
    
    vars = ['E_ALONG_V','dev_smooth0','dev_interp','dev_smooth1']
    ylim, vars, -800, 400, 0
    
    device, decomposed = 0 & loadct2, 43
    tplot, ['dB_fac_v','E_ALONG_V','dev_smooth0','dev_interp', $
        'dev_smooth1','dev_del']
end