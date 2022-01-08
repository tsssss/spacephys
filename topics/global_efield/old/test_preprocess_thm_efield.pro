;+
; Use thm l1 data to get the flags, similar to what i did for rbsp
;-


; good day.
fn = '/Users/Sheng/Downloads/thd_l1_eff_20080203_v01.cdf'
utr0 = time_double(['2008-02-03','2008-02-04'])

; bad day.
;fn = '/Users/Sheng/Downloads/thd_l1_eff_20150317_v01.cdf'
;utr0 = time_double(['2015-03-17','2015-03-18'])

cdf = scdfread(fn)
spinrate = 3d
boomlen = [50,40,6.9]*2

ut0s = *cdf[4].value
dr0 = sdatarate(ut0s)
euvw = float(*cdf[0].value)
;for i=0, 2 do euvw[*,i] = smooth(euvw[*,i],spinrate/dr0)
for i=0, 2 do euvw[*,i] = euvw[*,i]-scalcbg(euvw[*,i])
for i=0, 2 do euvw[*,i] = euvw[*,i]/boomlen[i]
euvw[*,2] = 0

store_data, 'euvw', ut0s, euvw, limits={labels:['u','v','w'], colors:[6,4,2]}

euvw[*,2] = 0
emag = sqrt(euvw[*,0]^2+euvw[*,1]^2)
store_data, 'emag', ut0s, emag

tplot, ['euvw','emag'], trange=utr0


stop


fn = '/Users/Sheng/Downloads/thd_l1_efw_20150317_v01.cdf'
utr0 = time_double(['2015-03-17','2015-03-18'])

fn = '/Users/Sheng/Downloads/thd_l1_efw_20080203_v01.cdf'
utr0 = time_double(['2008-02-03','2008-02-04'])


cdf = scdfread(fn)


euvw = *cdf[0].value
ut0s = *cdf[4].value
stop
nrec = n_elements(ut0s)
dr0 = 86400d/nrec
ut0s = smkarthm(utr0[0],utr0[1],dr0, 'dx')
ut0s = ut0s[0:nrec-1]

store_data, 'euvw', ut0s, euvw, limits={labels:['u','v','w'], colors:[6,4,2]}

euvw[*,2] = 0
emag = snorm(euvw)
store_data, 'emag', ut0s, emag

device, decomposed=0
loadct2, 43
tplot, ['euvw','emag'], trange=utr0

end
