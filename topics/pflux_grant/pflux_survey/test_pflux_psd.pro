;+
; Test to plot pflux PSD during storm time.
;-

dst_range = [-1e5,-40]
quiet_range = [-1,1]*5
time_range = time_double(['2012-10-01','2015-10-01'])
probe = 'a'
rbspx = 'rbsp'+probe
prefix = rbspx+'_'

; Load Dst.
data_type = 'dst'
pflux_survey_load_data, data_type, probe=probe
dst_var = prefix+data_type

; Load pflux PSD.
pflux_survey_load_psd_data, probe=probe
pf_psd_var = prefix+'pf_fac_psd'

get_data, pf_psd_var, times, pf_psd, freqs
dst = get_var_data(dst_var, at=times)

index = where_pro(dst, '[]', quiet_range, count=cnt2)
pf2 = pf_psd[index,*]

sgopen, 0, xsize=5, ysize=5

txs = 1d/freqs
xrange = [1,2000]
yrange = [1e-10,1e4]

bg_color = sgcolor('silver')

plot, xrange, yrange, position=sgcalcpos(1), /nodata, $
    xlog=1, ylog=1, xtitle='Period (sec)', ytitle='PSD (mW/m!U2!N)/Hz'

index = where_pro(dst, '[]', dst_range, count=cnt)
index = where_pro(dst, '[]', quiet_range, count=cnt)
pf = pf_psd[index,*]
for ii=0,cnt-1 do oplot, txs, pf[ii,*], color=bg_color, psym=1, symsize=0.2
oplot, txs, total(pf,1)/cnt

;oplot, 1d/freqs, total(pf2,1)/cnt2*10

end