
; check correlation b/w RBSP probe potentials





utr = time_double(['2013-05-01/07:35','2013-05-01/07:40'])
tprobe = 'b'

utr = time_double(['2012-11-14/04:41','2012-11-14/04:46'])
tprobe = 'a'


pre0 = 'rbsp'+tprobe+'_'

tmp = sread_rbsp_efw_l1(utr, probes = tprobe, type = ['vsvy'], vars = ['epoch','vsvy'])
uts = sfmepoch(tmp.epoch,'unix')
idx = where(uts ge utr[0] and uts le utr[1])
uts = uts[idx]
vsc = tmp.vsvy[idx,*]
store_data, pre0+'v1', uts, vsc[*,0]
store_data, pre0+'v2', uts,-vsc[*,1]
store_data, pre0+'v3', uts, vsc[*,2]
store_data, pre0+'v4', uts,-vsc[*,3]
store_data, pre0+'vsc', uts, vsc



; **** cross correlation b/w v1 and v2.

tvar = pre0+'vsc'
get_data, tvar, uts, bv



nrec = n_elements(uts)

dr = sdatarate(uts)

mtscl = 1    ; in sec.

mscl = ceil(mtscl/dr)

scls = smkarthm(-mscl,mscl,1,'dx')

nscl = n_elements(scls)

corr = dblarr(nrec,nscl)



tscl = 3/dr     ; 5 sec.

v1 = bv[*,0]
v2 = -bv[*,1]

for i = 0, nscl-1 do begin

    for j = tscl, nrec-tscl-1 do begin

        corr[j,i] = c_correlate(v1[j-tscl:j+tscl],v2[j-tscl:j+tscl],scls[i])

    endfor

endfor



store_data, 'corr_v12', uts, corr, scls*dr*1e3, limits = $

    {ytitle:'Time (msec)', spec:1, no_interp:1, zrange:[-1,1], $

    yrange:mtscl*[-1,1]*1e3, ztitle:'Cross Corr'}

    

maxcorr = dblarr(nrec)

for i = 0, nrec-1 do begin

    tmp = max(corr[i,*], idx)

    maxcorr[i] = scls[idx]*dr*1e3

endfor

store_data, 'maxcorr_v12', uts, maxcorr, limits = $

    {ytitle:'Time (msec)', yrange:mtscl*[-1,1]*1e3}

    

vars = ['corr_v12','maxcorr_v12',pre0+['v1','v2','v3','v4']]
tplot, vars
    

stop







end
