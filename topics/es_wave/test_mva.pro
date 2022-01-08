
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.35',':02.41'])  ; the largest sub-wave packet.

tprobe='a'
pre0 = 'rbsp'+tprobe+'_'
dr0 = 1d/4096
tpad = 200*dr0 ; sec.
p0 = 0.02   ; sec.
utr = utr0+[-1,1]*tpad
timespan, utr[0], utr[1]-utr[0], /second

deg = 180/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
rgb = sgcolor(['red','green','blue'])
red = sgcolor('red')

uvw = ['U','V','W']
bpv = ['b','p','v']
fac = ['para','west','north']

chsz = 1

tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'yticklen', 1



_2013_0607_0456_load_burst_data


; **** phase velocity along u.
; cross correlation to get dt.
p0 = 1d/70  ; 70 Hz.
get_data, pre0+'vb1', uts, vb1
for i = 0, 5 do vb1[*,i] = vb1[*,i]-smooth(vb1[*,i],p0/dr0)
v1 = vb1[*,0]
v2 = vb1[*,1]

dnrng = [-1,1]*10
dtrng = dnrng*dr0
dns = smkarthm(dnrng[0],dnrng[1],1,'dx')
ndn = n_elements(dns)
cors = dblarr(ndn)
dts = dns*dr0*1e3
for i = 0, ndn-1 do begin
    cors[i] = c_correlate(v1,v2,dns[i])
endfor

plot, dts, cors

dt1s = smkarthm(dts[0],dts[ndn-1],ndn*100,'n')
cor1s = interpol(cors,dts,dt1s,/quadratic)

stop


; **** rotate fields in "FAC".
b0p = 0.067  ; sec, the period in B0.
get_data, pre0+'eb1_f1', uts
get_data, pre0+'b0_uvw', tuts, dat
dat = sinterpol(dat, tuts, uts)
for i = 0, 2 do dat[*,i] = smooth(dat[*,i],b0p/dr0)
bhat = sunitvec(dat)


; ** method 1: rotation.
p = atan(bhat[*,1],bhat[*,0])   ; angle in GSE x-y plane.
t = atan(bhat[*,2],sqrt(bhat[*,1]^2+bhat[*,0]^2))   ; angle out of GSE x-y plane.

srotate, dat, -p, 2
srotate, dat,  t, 1
b0mag = dat[*,0]
store_data, pre0+'b0mag', uts, b0mag, limits = $
    {ynozero:1, ytitle:'(nT)', labels:'|B0|'}
    
get_data, pre0+'eb1_f1', uts, eb1f1
srotate, eb1f1, -p, 2
srotate, eb1f1,  t, 1
store_data, pre0+'eb1_f1_fac1', uts, eb1f1, limits = $
    {ytitle:'(mV/m)', labels:'dE '+bpv, colors:rgb}

get_data, pre0+'mb1_f1', uts, mb1f1    
srotate, mb1f1, -p, 2
srotate, mb1f1,  t, 1
store_data, pre0+'mb1_f1_fac1', uts, mb1f1, limits = $
    {ytitle:'10!U-3!N (nT)', labels:'dB '+bpv, colors:rgb}


; ** method 2: use position.
get_data, pre0+'eb1_f1', uts
get_data, pre0+'pos_gsm', tuts, dat
dat = sinterpol(dat, tuts, uts)
rhat = sunitvec(dat)
phat = sunitvec(scross(rhat,bhat))
vhat = scross(bhat,phat)

get_data, pre0+'eb1_f1', uts, dat
eb1f1 = [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]]
store_data, pre0+'eb1_f1_fac', uts, eb1f1, limits = $
    {ytitle:'(mV/m)', labels:'dE '+fac, colors:rgb}

get_data, pre0+'mb1_f1', uts, dat
mb1f1 = [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]]
store_data, pre0+'mb1_f1_fac', uts, mb1f1, limits = $
    {ytitle:'10!U-3!N (nT)', labels:'dB '+fac, colors:rgb}


ofn = shomedir()+'/b0_de_db_70hz.dat'
if file_test(ofn) eq 1 then file_delete, ofn
ssavebin, ofn, uts
ssavebin, ofn, eb1f1
ssavebin, ofn, mb1f1
ssavebin, ofn, b0mag

ofn = shomedir()+'/fig_b0_de_db.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
tplot, pre0+['eb1_f1_fac','mb1_f1_fac','b0mag'], trange = utr1
sgclose

stop




options, pre0+'mb1_*', 'ystyle', 1
options, pre0+'mb1_*', 'yrange', [-1,1]*70
options, pre0+'eb1_*', 'ystyle', 1
options, pre0+'eb1_*', 'yrange', [-1,1]*70



; **** plot3d.
ofn = shomedir()+'/fig_3d_hodogram.pdf'
;ofn = 0
sgopen, ofn, xsize = 5, ysize = 5, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


dbrng = [-1,1]*80
get_data, pre0+'mb1_f1_fac', uts, dat
idx = where(uts ge utr2[0] and uts le utr2[1])
dat = dat[idx,*]
uts = uts[idx]

surface, dist(5), /nodata, /save, $
    xrange = dbrng, yrange = dbrng, zrange = dbrng, $
    position = [0.2,0.2,0.8,0.8], $
    xstyle = 1, ystyle = 1, zstyle = 1, $
    xtitle = 'dB'+bpv[0]+' 10!U-3!N (nT)', $
    ytitle = 'dB'+bpv[1]+' 10!U-3!N (nT)', $
    ztitle = 'dB'+bpv[2]+' 10!U-3!N (nT)', $
;    az = -70, ax = 20, $    ; view to show the wave vector.
    az = 15, ax = 40, $     ; view to show the polarization plane.
    charsize = chsz*2

plots, dat[*,0], dat[*,1], dat[*,2], /t3d
nrec = n_elements(uts)
for i = 0, nrec-1 do begin
    plots, dat[i,0], dat[i,1], dat[i,2], /t3d, psym = 4, symsize = chsz*0.4, $
        color = sgcolor(254d*(1-float(i)/(nrec-1)), ct = 40)
endfor
plots, dat[0,0], dat[0,1], dat[0,2], /t3d, psym = 6, color = red
;for i = 0, nrec-2 do begin
;    plots, dat[i:i+1,0], dat[i:i+1,1], dat[i:i+1,2], /t3d, $
;        color = sgcolor(254d*(1-float(i)/(nrec-1)), ct = 40)
;endfor

tx0 = 0.7
ty0 = 0.75
plots, tx0, ty0+ychsz*0.2, /normal, psym = 6, color = red, symsize = chsz*0.8
xyouts, tx0+xchsz*1, ty0, /normal, alignment = 0, '  t = 0', charsize = chsz

eigvs = smva(dat, bvecs)

ndir = bvecs[0,*]       ; normal variation.
mdir = bvecs[1,*]       ; median.
ldir = bvecs[2,*]       ; largest.

labs = ['normal (minimum)','median','maximum']
for i = 0, 2 do begin
    tmp = bvecs[i,*]*max(dbrng)
    plots, [-1,1]*tmp[0],[-1,1]*tmp[1], [-1,1]*tmp[2], color = rgb[i], /t3d
    plots, tx0+[-0.4,1]*xchsz, ty0-ychsz*1.2*(i+1)+ychsz*0.2, color = rgb[i], /normal
    xyouts, tx0+xchsz*1, ty0-ychsz*1.2*(i+1), '  '+labs[i], /normal
    
endfor

; red-min,normal
; green-medien
; blue-max.

sgclose


end
