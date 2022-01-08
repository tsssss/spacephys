;+
; Test HOPE odd moments, compare dB to current.
;-
;

rgb = [6,4,2]
xyz = ['x','y','z']
qe = 1.6e-19

tplot_options, 'lab_flag', -1


utr0 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = time_double(['2013-04-14/04:45','2013-04-14/05:15'])
tprobe = 'b'

pre0 = 'rbsp'+tprobe+'_'

hopemom = sread_rbsp_hope_moments(utr0, probe=tprobe)
uts = [hopemom.ut_ion, hopemom.ut_ele]
uts = uts[sort(uts)]
nrec = n_elements(uts)
jhope = dblarr(nrec,3)
nflux = sinterpol(hopemom.e_number_flux, hopemom.ut_ele, uts)-$
    sinterpol(hopemom.p_number_flux, hopemom.ut_ion, uts)-$
    sinterpol(hopemom.o_number_flux, hopemom.ut_ion, uts)-$
    sinterpol(hopemom.he_number_flux, hopemom.ut_ion, uts)
store_data, pre0+'jhope', uts, nflux, limits={ytitel:'#/cm!U2!N-s', colors:rgb, labels:'GSM '+xyz}

emfisis = sread_rbsp_emfisis_l3(utr0, probe=tprobe, coord='gsm', type='hires')
tuts = sfmepoch(emfisis.epoch, 'unix')
bgsm = emfisis.mag
b0gsm = bgsm & for i=0,2 do b0gsm[*,i] = scalcbg(bgsm[*,i])
dbgsm = bgsm-b0gsm
store_data, pre0+'db_gsm', tuts, dbgsm, limits={ytitle:'(nT)', colors:rgb, labels:'GSM B'+xyz}
store_data, pre0+'by_gsm', tuts, bgsm[*,1], limits={ytitle:'(nT)'}

tplot, ['?_number_flux','o_vbulk', pre0+'jhope',pre0+'by_gsm'], trange=utr0
stop

end