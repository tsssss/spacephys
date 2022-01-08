
;; uncomment this block to load infos.
;infos = [] 
;ids = cusp_id('south_imf')
;foreach id, ids do begin
;    tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 1)
;    infos = [infos,tinfo]
;endforeach
;
;ids = cusp_id('south_imf')
;foreach id, ids do begin
;    tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 0)
;    infos = [infos,tinfo]
;endforeach


scales = []
pfluxs = []
ebrats = []
diss = []
ninfo = n_elements(infos)
for i = 0, ninfo-1 do begin
    tmp = infos[i].nfilter
    scales = [scales, infos[i].filters[0:tmp-1]]
    pfluxs = [pfluxs, infos[i].pfluxes[0:tmp-1]]
    ebrats = [ebrats, infos[i].ebratio[0:tmp-1]]
    diss = [diss, dblarr(tmp)+infos[i].dis]
endfor

idx = where(diss le 3)
scales = scales[idx]
pfluxs = pfluxs[idx]
ebrats = ebrats[idx]
diss = diss[idx]

tx = alog(scales)
ty = alog(pfluxs)
idx = where(finite(ty))
tx = tx[idx]
ty = ty[idx]
coef = linfit(tx, ty)
print, coef
plot, tx, ty, psym = 1

txs = smkarthm(0,8,20,'n')
tys = coef[0]+coef[1]*txs
oplot, txs, tys

stop

tx = alog(scales)
ty = alog(ebrats)
idx = where(finite(ty))
tx = tx[idx]
ty = ty[idx]
coef = linfit(tx, ty)
print, coef
plot, tx, ty, psym = 1

txs = smkarthm(0,8,20,'n')
tys = coef[0]+coef[1]*txs
oplot, txs, tys

stop

end
