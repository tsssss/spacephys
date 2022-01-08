;+
; test on solving k and omega, given dE an dB in time series.
;-

deg = 180d/!dpi
rad = !dpi/180d
ndim = 3
lmn = ['L','M','N']
rgb = sgcolor(['red','green','blue'])

tplot_options, 'labflag', -1

nrec = 1001 ; # of records.
dt = 1d     ; duration, in sec.
uts = findgen(nrec)/(nrec-1)
dr0 = uts[1]-uts[0]  ; 1 kHz sample rate.
utr = minmax(uts)

; in LMN coord.
k0 = [0d,0,1]*1.3   ; wave vector in normal direction.
p0 = dt*0.2     ; period, in sec.
f0 = 1d/p0      ; frequency, in Hz.
w0 = 2*!dpi/p0  ; angular frequency, in rad/Hz.


; E in MN plane, almost linearly polarized along N.
de0 = [10,20,100]        ; real amplitude, in V/m.
ep0 = [10,100,100]*rad   ; initial phase, in rad.
des = dblarr(nrec,ndim)
for i = 0, ndim-1 do des[*,i] = de0[i]*cos(w0*uts+ep0[i])
store_data, 'delmn', uts, des, limits = {ytitle:'dE LMN', labels:lmn, colors:rgb}


; analytical results of E and B in frequency domain.
ie0 = complexarr(ndim)
for i = 0, ndim-1 do ie0[i] = complex(de0[i]*cos(ep0[i]),de0[i]*sin(ep0[i]))
ib0 = complexarr(ndim)
trb = scross(k0,real_part(ie0))/w0
tib = scross(k0,imaginary(ie0))/w0
for i = 0, ndim-1 do ib0[i] = complex(trb[i],tib[i])
print, '|Bx|,|By|,|Bz|: ', sqrt(trb^2+tib^2)
print, 'Phase: x, y, z (deg): ', atan(tib,trb)*deg



; B in LM plane, calculated from Faraday's law.
dbs = dblarr(nrec,ndim)
for i = 0, nrec-1 do dbs[i,*] = scross(k0,reform(des[i,*]))
dbs = dbs/w0
store_data, 'dblmn', uts, dbs, limits = {ytitle:'dB LMN', labels:lmn, colors:rgb}


sgopen, 0

vars = ['de','db']+'lmn'
nvar = n_elements(vars)
poss = sgcalcpos(nvar)
tplot, vars, position = poss, trange = utr




; **** step 1: determine w1.
w1s = dblarr(ndim)  ; angular frequency in rad/Hz.
dat = des
w = 6d
s2p = 4*!dpi/(w+sqrt(2+w^2))
s0 = 5d
s1 = nrec/2d
nj = 60d
dj = alog(exp(alog(s1/s0)/nj))/alog(2)
ntp = 40d       ; # of periods to be tested.

cdelta = 0.776d
gamma = 2.32
psi0 = !dpi^(-0.25)


for i = 0, ndim-1 do begin
    tdat = dat[*,i]
    mor = wv_cwt(tdat, 'Morlet', w, /pad, $
        start_scale = s0, nscale = nj, dscale = dj, scale = recscls)
    timescls = recscls*dr0
    periods = timescls*s2p
    morpow = abs(mor^2)*2*!dpi*(uts[nrec-1]-uts[0])/nrec
    morspec = total(morpow,1)/nrec
    maxspec = max(morspec, widx)
    w1s[i] = 2*!dpi/periods[widx]
endfor
w1 = mean(w1s)

print, 'original w (rad/Hz):', w0
print, 'calculated w:       ', w1


; **** step 2: determine direction of k1.
eigvs = smva(dbs, bvecs)

k1hat = reform(bvecs[2,*])


; **** step 3: relative phase within components.
p1 = 2*!dpi/w1
ntp = 50
tps = smkarthm(0,p1,ntp,'n')
tns = tps/dr0

ep1 = dblarr(ndim)      ; relative phase, in deg.
tdat = des
for i = 1, ndim-1 do begin
    cor = dblarr(ntp)
    for j = 0, ntp-1 do cor[j] = c_correlate(tdat[*,0],tdat[*,i],-tns[j])
    maxcor = max(cor, pidx)
    ep1[i] = tps[pidx]*w1*deg
endfor
print, 'relative phase of dE (deg):', ep1


bp1 = dblarr(ndim)      ; relative phase, in deg.
tdat = dbs
for i = 1, ndim-1 do begin
    cor = dblarr(ntp)
    for j = 0, ntp-1 do cor[j] = c_correlate(tdat[*,0],tdat[*,i],-tns[j])
    maxcor = max(cor, pidx)
    bp1[i] = tps[pidx]*w1*deg
endfor
print, 'relative phase of dB (deg):', bp1


cor = dblarr(ntp)
for j = 0, ntp-1 do cor[j] = c_correlate(dbs[*,0],des[*,0],-tns[j])
maxcor = max(cor, pidx)
peb = tps[pidx]*w1*deg      ; phase dE-dB = peb.
print, 'relative phase b/w dE and dB (deg):', peb



; **** step 4: amplitude.
de1s = dblarr(nrec,ndim)
dat = des
for i = 0, ndim-1 do begin
    tdat = deriv(dat[*,i])/(dr0*w1)
    de1s[*,i] = sqrt(tdat^2+dat[*,i]^2)
    de1s[*,i] = smooth(de1s[*,i], p1/dr0, /edge_truncate)
endfor
print, 'original amplitude of dE:', de0
print, 'calculated amplitude:    ', reform(de1s[0,*])

db1s = dblarr(nrec,ndim)
dat = dbs
for i = 0, ndim-1 do begin
    tdat = deriv(dat[*,i])/(dr0*w1)
    db1s[*,i] = sqrt(tdat^2+dat[*,i]^2)
    db1s[*,i] = smooth(db1s[*,i], p1/dr0, /edge_truncate)
endfor



; **** step 5: determine |k|.
kxe = complexarr(ndim)      ; k cross E.
wdb = complexarr(ndim)      ; w dot B.
tidx = nrec/2
de1 = reform(de1s[tidx,*])
dde1 = complex(de1*cos(ep1*rad),de1*sin(ep1*rad))
db1 = reform(db1s[tidx,*])
ddb1 = complex(db1*cos(bp1*rad),db1*sin(bp1*rad))
kxe = complex(scross(k1hat,real_part(dde1)), $
    scross(k1hat,imaginary(dde1)))
wdb = w1*ddb1


res = linfit(abs(kxe), abs(wdb))
k1 = res[1]*k1hat

plot, abs(kxe), abs(wdb), psym = -1

stop



end
