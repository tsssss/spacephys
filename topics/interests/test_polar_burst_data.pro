
fn = shomedir()+'/Google Drive/works/data/misc/po_burst_fld_1996_0318_v01.sdt'
tr = time_double(['1996-03-18/08:37:20','1996-03-18/08:38:20'])

if file_test(fn) eq 0 then message, 'file does not exist ...'
sdt = ssdtread(fn)
pre = 'po_'


if n_elements(orootdir) eq 0 then orootdir = shomedir()
if ~file_test(orootdir,/directory) then file_mkdir, orootdir

; some quantities and settings.
device, decomposed = 0 & loadct2, 43
red = 6
green = 4
blue = 2
black = 0
rgb = [red,green,blue]   ; r,g,b in colortable 43,45.
labfac = ['v','p','b']   ; for north sunward meridian cross
labspc  = ['xy','z','56']       ; v~xy, p(vxb)~56, b~z.
defactor = 1.3d ; correction to E due to shielding.
!p.font = 1
tplot_options, 'ygap', 0.25
tplot_options, 'ynozero', 1
tplot_options, 'version', 2
tplot_options, 'num_lab_min', 8
tplot_options, 'labflag', -1
time_stamp, /off
ct = 43
posd = [0.15,0.10,0.9,0.45]     ; lower half.
posu = [0.15,0.55,0.9,0.90]     ; upper half.
posl = [0.10,0.1,0.30,0.9]      ; left.
posm = [0.40,0.1,0.60,0.9]      ; middle.
posr = [0.70,0.1,0.90,0.9]      ; right.

if n_elements(titpre) eq 0 then titpre = ''

;**** get uniform time.
maxnrec = 100000ul   ; prevent memory overflow.
t0 = sdt.var.polar_b_spc_z.depend_0
dr = sdatarate(t0) & nrec = n_elements(t0)
if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
tmp = minmax(sdt.var.polar_e_spc_z.depend_0)
t0 = smkarthm(max([tmp[0],t0[0]]),min([tmp[1],t0[nrec-1]]), dr, 'dx')
nrec = n_elements(t0)
tstr = time_string(t0[0], tformat='YYYY_MMDD')
if n_elements(eventid) ne 0 then tstr = eventid
print, 'data rate: ', dr

;**** original b field and spike removal.
ft  = sdt.var.polar_b_spc_z.depend_0
fxy = sdt.var.polar_b_spc_x_y.value
f56 = sdt.var.polar_b_spc_56.value
fz  = sdt.var.polar_b_spc_z.value
b_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
b_spc = sinterpol(b_spc, ft, t0)
; raw B field data, only interpolated to uniform time.
store_data, 'bxy', t0, b_spc[*,0], limits = {labels:'B0xy'}
store_data, 'bz' , t0, b_spc[*,1], limits = {labels:'B0z'}
store_data, 'b56', t0, b_spc[*,2], limits = {labels:'B056'}
; Bxy and Bz have discontinuities and spikes, B56 has short spikes.
sdespike, t0, b_spc, width = 100, _extra = extra
store_data, 'b0xy', t0, b_spc[*,0], limits = {labels:'Bxy'}
store_data, 'b0z' , t0, b_spc[*,1], limits = {labels:'Bz'}
store_data, 'b056', t0, b_spc[*,2], limits = {labels:'B56'}

;**** plot1: original field and spike removal.
vars = ['bxy','bz','b56','b0xy','b0z','b056']
options, vars, 'ytitle', '(nT)'
ylim, ['bxy','b0xy'], min(b_spc[*,0])-20, max(b_spc[*,0])+20, 0
ylim, ['bz','b0z'], min(b_spc[*,1])-20, max(b_spc[*,1])+20, 0
ylim, ['b56','b056'], min(b_spc[*,2])-20, max(b_spc[*,2])+20, 0
fn = orootdir+'/'+tstr+'_polar_b_despike.pdf'
sgopen, fn, xsize = 5, ysize = 7, /inch
sgindexcolor
loadct2, ct
vars = ['bxy','bz','b56']
nvar = n_elements(vars)
titl = 'Polar B_SPC, raw data'
poss = sgcalcpos(nvar, position = posu)
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
vars = ['b0xy','b0z','b056']
nvar = n_elements(vars)
titl = 'Polar B_SPC, after despike'
poss = sgcalcpos(nvar, position = posd)
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
xyouts, 0.5, 0.95, titpre+'Polar original and despiked B!DSPC!N', /normal, alignment = 0.5, charsize = 1.25
sgclose
if keyword_set(noplot) then file_delete, strmid(fn,0,strlen(fn)-3)+'pdf'
vars = ['bxy','bz','b56','b0xy','b0z','b056']
store_data, vars, /delete

;**** total b field. 'po_b'
btotal = sqrt(total(b_spc^2,2))
store_data, pre+'b', t0, btotal, limits = {ytitle:'B mag!C(nT)', ynozero:1}
store_data, pre+'b_spc', t0, b_spc, $
    limits = {ytitle:'B SPC!C(nT)', labels:labspc, colors:rgb}

;**** t96 model.
ft  = sdt.var.polar_model_b_t96_spc_z.depend_0
fxy = sdt.var.polar_model_b_t96_spc_x_y.value
f56 = sdt.var.polar_model_b_t96_spc_56.value
fz  = sdt.var.polar_model_b_t96_spc_z.value
bt96_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
bt96_spc = sinterpol(bt96_spc, ft, t0)

;**** model b field. 'po_b0_spc'.
db_spc = b_spc-bt96_spc
bmod_spc = scalcbg(db_spc)+bt96_spc
store_data, pre+'b0_spc', t0, bmod_spc, $
    limits = {ytitle:'B model SPC!C(nT)', labels:labspc, colors:rgb}
    
;**** db field. 'po_db_spc'.
db_spc = b_spc-bmod_spc
 idx = where(abs(db_spc) gt 500, cnt)
 if cnt ne 0 then db_spc[idx] = 0
store_data, pre+'db_spc', t0, db_spc, $
    limits = {ytitle:'dB SPC!C(nT)', labels:labspc, colors:rgb}
store_data, pre+'b0_spc', t0, bmod_spc

;**** de field. 'po_de_spc'.
ft  = sdt.var.polar_e_spc_z.depend_0
fxy = sdt.var.polar_e_spc_x_y.value
f56 = sdt.var.polar_e_spc_56.value
fz  = sdt.var.polar_e_spc_z.value
de_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
de_spc = sinterpol(de_spc, ft, t0)*defactor
de560 = de_spc[*,0]
; e56_dot0.
de56_dot0 = (b_spc[*,0]*de_spc[*,0]+b_spc[*,1]*de_spc[*,1])/b_spc[*,2]
idx = where(abs(b_spc[*,2]/btotal) le 0.2, cnt)
if cnt ne 0 then de56_dot0[idx] = !values.f_nan
; throw e56 by default, set e56 to keep e56, set edot0 to calc e56.
if ~keyword_set(e56) then de_spc[*,2] = 0
if keyword_set(edot0) then de_spc[*,2] = de56_dot0
store_data, pre+'de_spc', t0, de_spc, $
    limits = {ytitle:'dE SPC!C(mV/m)', labels:labspc, colors:rgb}
store_data, 'de56', t0, de_spc[*,2]

;**** plot2: B model on B, dB, original dE and
; B total, dB 3D before correction, E56_dot0.
fn = orootdir+'/'+tstr+'_polar_e56dot0_dbcorr.pdf'
sgopen, fn, xsize = 10, ysize = 7, /inch
sgindexcolor
loadct2, ct
; plot B model on top of B total.
store_data, 'b0xy', t0, [[bmod_spc[*,0]],[b_spc[*,0]],[bt96_spc[*,0]]], $
    limits = {labels:'Bxy'}
store_data, 'b0z', t0, [[bmod_spc[*,1]],[b_spc[*,1]],[bt96_spc[*,1]]], $
    limits = {labels:'Bz'}
store_data, 'b056', t0, [[bmod_spc[*,2]],[b_spc[*,2]],[bt96_spc[*,2]]], $
    limits = {labels:'B56'}
vars = ['b0xy','b0z','b056'] & nvar = n_elements(vars)
titl = 'Polar B model SPC'
poss = posl & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'colors', [red,black,green]
options, vars, 'ytitle', '(nT)'
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; plot dB SPC in components.
store_data, 'dbxy', t0, db_spc[*,0], limits = {labels:'dBxy'}
store_data, 'dbz', t0, db_spc[*,1], limits = {labels:'dBz'}
store_data, 'db56', t0, db_spc[*,2], limits = {labels:'dB56'}
vars = ['dbxy','dbz','db56'] & nvar = n_elements(vars)
titl = 'Polar dB SPC'
poss = posm & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(nT)'
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; plot dE SPC in components.
store_data, 'dexy', t0, de_spc[*,0], limits = {labels:'dE0xy'}
store_data, 'dez' , t0, de_spc[*,1], limits = {labels:'dE0z'}
store_data, 'de56', t0, de560, limits = {labels:'dE056'}
vars = ['dexy','dez','de56'] & nvar = n_elements(vars)
titl = 'Polar dE SPC'
poss = posr & poss[1] = 0.40 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(mV/m)'
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete & de560 = 0
; plot dE56 dot0.
store_data, 'de56_dot0', t0, de56_dot0
vars = ['de56_dot0'] & nvar = n_elements(vars)
titl = 'Polar dE56_dot0'
poss = posr & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(mV/m)'
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; plot dB SPC before correction.
store_data, 'db', t0, db_spc, limits = {ytitle:'(nT)', $
    labels:labspc, colors:rgb}
vars = ['db'] & nvar = n_elements(vars)
titl = 'Polar dB SPC 3D'
poss = posm & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(nT)'
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; plot B magnitude.
tmp = sqrt(total(bmod_spc^2,2))
tmp2 = sqrt(total(bt96_spc^2,2))
store_data, 'btotal', t0, [[tmp],[btotal],[tmp2]], limits = $
    {colors:[red,black,green], labels:['Bsmth','in situ','T96']}
vars = ['btotal'] & nvar = n_elements(vars)
titl = 'Polar B magnitude'
poss = posl & poss[3] = 0.25 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(nT)'
options, vars, 'ystyle', 1
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
xyouts, 0.5, 0.95, titpre+'Polar original B and E, B!Dmodel!N dB separation, dE!Ddot0!N', /normal, alignment = 0.5, charsize = 1.25
sgclose
if keyword_set(noplot) then file_delete, strmid(fn,0,strlen(fn)-3)+'pdf'

;**** ilat, mlt, dis.
t1  = sdt.var.polarinvariantlatitude.depend_0
ilat = sdt.var.polarinvariantlatitude.value
store_data, pre+'ilat', data = {x:t1, y:ilat}, limits = {ytitle:'ILat'}
mlt  = sdt.var.polarmlt.value
store_data, pre+'mlt', data = {x:t1, y:mlt}, limits = {ytitle:'MLT'}
dis  = sdt.var.polarspcraftdist.value
store_data, pre+'dis', data = {x:t1, y:dis}, limits = {ytitle:'Dist (Re)'}

;**** rotate from SPC to FAC.
polar_sdt_spc2fac, bmod_spc, db_spc, de_spc, db_fac, de_fac, b1 = lon, b2 = lat
store_data, pre+'db_fac', t0, db_fac, $
    limits = {ytitle:'dB FAC!C(nT)', labels:labfac, colors:rgb}
store_data, pre+'de_fac', t0, de_fac, $
    limits = {ytitle:'dE FAC!C(mV/m)', labels:labfac, colors:rgb}
store_data, pre+'spc2fac', t0, [[lon],[lat]]*(180/!dpi), limits = $
    {ytitle:'SPC2FAC!C(deg)', labels:['lon','lat'], colors:rgb[0:1]}
    
; plot3: show dB and dE rotated from SPC to FAC.
fn = orootdir+'/'+tstr+'_polar_spc2fac.pdf'
sgopen, fn, xsize = 7, ysize = 7, /inch
sgindexcolor
loadct2, ct
; dE FAC.
vars = 'de_fac'+labfac
nvar = n_elements(vars)
stplot_split, pre+'de_fac', newname = vars
tminmax = [0d,0d]
for i = 0, nvar-1 do begin
    twd = 10d/dr         ; 10 sec smooth.
    get_data, vars[i], t0, tmp
    dat = smooth(tmp, twd, /edge_mirror)
    store_data, vars[i], t0, [[tmp],[dat]], $
        limits = {colors:[0,red], labels:'dE'+labfac[i]}
    tminmax = minmax([tminmax,dat])
endfor
titl = 'Polar dE FAC'
poss = posu & poss[2] = 0.45 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(mV/m)'
options, vars, 'yrange', tminmax
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; dB FAC.
vars = 'db_fac'+labfac
nvar = n_elements(vars)
options, pre+'db_fac', 'labels', 'dB'+labfac
stplot_split, pre+'db_fac', newname = vars
get_data, pre+'db_fac', t0, tmp
tminmax = minmax(tmp)
titl = 'Polar dB FAC'
poss = posu & poss[0] = 0.60 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(nT)'
options, vars, 'yrange', tminmax
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; dE SPC.
vars = 'de_spc'+labspc
nvar = n_elements(vars)
stplot_split, pre+'de_spc', newname = vars
tminmax = [0d,0d]
for i = 0, nvar-1 do begin
    twd = 10d/dr         ; 60 sec smooth.
    get_data, vars[i], t0, tmp
    dat = smooth(tmp, twd, /edge_mirror)
    store_data, vars[i], t0, [[tmp],[dat]], $
        limits = {colors:[0,red], labels:'dE'+labspc[i]}
    tminmax = minmax([tminmax,dat])
endfor
titl = 'Polar dE SPC'
poss = posd & poss[2] = 0.45 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(mV/m)'
options, vars, 'yrange', tminmax
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
; dB SPC.
vars = 'db_spc'+labspc
nvar = n_elements(vars)
options, pre+'db_spc', 'labels', 'dB'+labspc
stplot_split, pre+'db_spc', newname = vars
get_data, pre+'db_spc', t0, tmp
tminmax = minmax(tmp)
titl = 'Polar dB SPC'
poss = posd & poss[0] = 0.60 & poss = sgcalcpos(nvar, position = poss)
options, vars, 'ytitle', '(nT)'
options, vars, 'yrange', tminmax
tplot, vars, title = titl, trange = tr, position = poss, /noerase
if n_elements(cusptr) ne 0 then timebar, cusptr, color = red
store_data, vars, /delete
xyouts, 0.5, 0.95, titpre+'Polar SPC to FAC rotation of dE and dB', /normal, alignment = 0.5, charsize = 1.25
sgclose
if keyword_set(noplot) then file_delete, strmid(fn,0,strlen(fn)-3)+'pdf'


; plot 4: burst E and V.
t0 = sdt.var.polar_e_mbl_spc_z.depend_0
idx = where(t0 gt tr[0] and t0 le tr[1], nrec)
t0 = t0[idx]
ft = t0
dr = sdatarate(t0)
if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
t0 = smkarthm(t0[0],t0[nrec-1], nrec, 'n')

fxy = sdt.var.polar_e_mbl_spc_x_y.value[idx]
f56 = sdt.var.polar_e_mbl_spc_56.value[idx]
fz  = sdt.var.polar_e_mbl_spc_z.value[idx]
e_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
e_spc = sinterpol(e_spc, ft, t0)
; raw B field data, only interpolated to uniform time.
store_data, 'burstexy', t0, e_spc[*,0], limits = {labels:'Burst Exy'}
store_data, 'burstez' , t0, e_spc[*,1], limits = {labels:'Burst Ez'}
store_data, 'burste56', t0, e_spc[*,2], limits = {labels:'Burst E56'}
store_data, pre+'burst_de_spc', t0, e_spc, $
    limits = {ytitle:'Burst dE SPC!C(mV/m)', labels:labspc, colors:rgb}

dat = sdt.var.mb_v1l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv1', t0, dat, limits = {labels:'Burst V1'}
dat = sdt.var.mb_v2l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv2', t0, dat, limits = {labels:'Burst V2'}
dat = sdt.var.mb_v3l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv3', t0, dat, limits = {labels:'Burst V3'}
dat = sdt.var.mb_v4l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv4', t0, dat, limits = {labels:'Burst V4'}
dat = sdt.var.mb_v5l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv5', t0, dat, limits = {labels:'Burst V5'}
dat = sdt.var.mb_v6l.value[idx]
dat = sinterpol(dat, ft, t0)
store_data, 'burstv6', t0, dat, limits = {labels:'Burst V6'}

vars = 'burst'+['v1','v2','v3','v4','v5','v6']
stplot_merge, vars, newname = pre+'burst_v', colors = [1,2,3,4,5,6]

fn = orootdir+'/'+tstr+'_polar_burst_e_v.pdf'
sgopen, fn, xsize = 7, ysize = 7, /inch
sgindexcolor
loadct2, ct

vars = 'burst'+['exy','ez','e56','v1','v2','v3','v4','v5','v6']
tplot, vars, trange = tr, title = 'Polar Burst E and V'

sgclose

vars = 'burst'+['exy','ez','e56','v1','v2','v3','v4','v5','v6']
store_data, vars, /delete


; check MAT spectrogram for E56.
nscale = 100
trange = reverse(1d/[50,800])    ; Hz
trange = [1d/1000,0.1]
tscales = smkarthm(trange[1],trange[0], nscale, 'n') ; 5 to 800 Hz.
tscales = smkgmtrc(trange[1],trange[0], nscale, 'n') ; 5 to 800 Hz.
tvar = 'po_burst_de_spc'
get_data, tvar, uts, dat
tvar = 'po_burst_de56'
store_data, tvar, uts, dat[*,2], limits = {labels:'Burst E56'}
stplot_mat, tvar, scale = tscales, newname = tvar+'_mat'
tvar = 'po_burst_de56_mat'
options, tvar, 'zrange', [-30,30]
options, tvar, 'yrange', [min(trange),max(trange)]
options, tvar, 'ylog', 1

;; change to Hz.
;get_data, tvar, uts, dat, val
;store_data, tvar, uts, reverse(dat,2), reverse(1d/val)
;options, tvar, 'yrange', 1d/trange

;vars = 'po_burst_de56'+['','_mat']
;tplot, vars


; **** rotate burst data into FAC.
get_data, 'po_burst_de_spc', uts, dat
nrec = n_elements(uts)
utr = minmax(uts)

get_data, 'po_b_spc', tuts, bhat
bhat = sunitvec(bhat)
bhat = sinterpol(bhat, tuts, uts)
ponhat = [[-bhat[*,1]],[bhat[*,0]],[dblarr(nrec)]]
poffhat = scross(bhat,ponhat)

dat = [[sdot(dat,bhat)],[sdot(dat,ponhat)],[sdot(dat,poffhat)]]
store_data, 'po_burst_de_fac', uts, dat, $
    limits = {ytitle:'Burst dE FAC!C(mV/m)', labels:['Epar','Eperon','Eperoff'], colors:rgb}


tvar = 'po_burst_de_fac'
vars = 'po_burst_de_fac_'+['b','on','off']
stplot_split, tvar, newname = vars
options, vars, 'yrange', [-400,400]
tplot, vars, trange = utr



; **** power spectrum.
tvar = 'po_burst_de_fac'
get_data, tvar, uts, dat
dr = sdatarate(uts)
nrec = n_elements(uts)

device, decomposed = 0
loadct2, 43
!p.charsize = 2
xr = [6e-2,600]
yr = [1e-4,1e2]
xr = [6e-2,600]
yr = [1e-6,1e5]
plot, xr, yr, /nodata, /xlog, /ylog, $
    xstyle = 1, ystyle = 1, xtitle = 'f (Hz)', ytitle = 'PS (mV/m)!U2!N/Hz'
fs = findgen(nrec)/(nrec*dr)
colors = [6,4,2]
tau = 2*!dpi*(utr[1]-utr[0])
for i = 0, 2 do begin
    tdat = dat[*,i]
    ps = abs(fft(tdat))
    ps = smooth(ps,8)
    ps = ps^2*tau
    oplot, fs, ps, color = colors[i]
endfor


stop


; **** cross correlation b/w v1 and v2.
get_data, 'po_burst_v', uts, bv

nrec = n_elements(uts)
dr = sdatarate(uts)
mtscl = 0.01    ; 10 ms.
mscl = ceil(mtscl/dr)
scls = smkarthm(-mscl,mscl,1,'dx')
nscl = n_elements(scls)
corr = dblarr(nrec,nscl)

tscl = mscl*10
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



stop



end

