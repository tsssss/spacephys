;+
; Type: crib.
; Purpose: Try to reproduce L3 data from L2 data, determine the time and angle
;   of each pixel at each energy bin.
; Parameters:
;   date, in, string, req. 'YYYY-MM-DD'.
;   ut0, in, string, req. 'YYYY-MM-DD/hh:mm', the time to be plotted.
;   probe, in, string, req. 'a','b'.
;   ion, in, string, req. 'proton','oxygen','helium'.
; Keywords: none.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-11, Sheng Tian, create.
;-


; date and time.
date = '2013-04-14'
date = '2014-02-19'
ut0 = '2013-05-01/07:39'
ut0 = '2012-11-14/04:42'
ut0 = time_double(ut0)    ; sample time.
en0 = 1e3                 ; sample energy.
type = 'proton'           ; sample species.
utr = ut0-(ut0 mod 86400)+[0,86400d]      ; time range to load data.
probe = ['b']

case type of
    'electron': type0 = 'FEDU'
    'proton': type0 = 'FPDU'
    'oxygen': type0 = 'FODU'
    'helium': type0 = 'FHEDU'
endcase


; suppress math exception.
!except = 0

; constants.
npxl = 5
nenergy = 72
nsec = 16

tsettle = 1.302e-3  ; 1.302 ms, voltage settle time for each esa sweep.
tmeasure = 9.115e-3 ; 9.115 ms, data acquire time for each esa sweep.
tesastep = tsettle+tmeasure ; 10.417 ms.

rad = !dpi/180d
deg = 180d/!dpi

; settings.
sgtruecolor
!p.background = sgcolor('white')
timespan, date, 1, /day

tplot_options, 'labflag', -1
tplot_options, 'xmargin', [20,15]
tplot_options, 'constant', 0


; L3 pitch angles.
l3pa0s = [4.5,18,36,54,72,90,108,126,144,162,175.5]
dpa0s = [4.5,9,9,9,9,9,9,9,9,9,4.5]
nl3pa0 = n_elements(l3pa0s)


; **** load data.
load = 0
vars = ['hopel2','hopel3']
for i = 0, n_elements(vars)-1 do if tnames(vars[i]) eq '' then load = 1

; b_uvw. B in uvw coord, directly comnparable to pixel direction.
; hopel2. level 2 data.
; hopel3. level 3 data.

if load eq 1 then begin

    ; load hope l2 and l3 data.
    hopel2 = sread_rbsp_hope_l2(utr, probes = probe)
    store_data, 'hopel2', 0, hopel2

    hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
    store_data, 'hopel3', 0, hopel3
    
    ; load spice kernal.
    rbsp_load_spice_kernels, /all
    
    ; load b field at time of interest.
    emfisis = sread_rbsp_emfisis(ut0, probes = probe)
    
    b_uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    b_gse = emfisis.mag
    
    store_data, 'b_gse', b_uts, b_gse
    rbsp_uvw2gse, 'b_gse', newname = 'b_uvw', probe = prob, /inverse, /no_spice_load
endif


; load the l2 and l3 data from tplot var.
get_data, 'hopel2', tmp, hopel2
get_data, 'hopel3', tmp, hopel3

; pick out the data at sample time and sample energy and sample species.
uts = sfmepoch(hopel3.epoch_ion,'unix')
tmp = min(uts-ut0, tidx, /absolute)
en = hopel3.hope_energy_ion[tidx,*]
tmp = min(en-en0, eidx, /absolute)
en = en[eidx]   ; energy, same for L3 and L2.
ut = uts[tidx]  ; time in the middle of 1 spin, same for L3 and L2.
dt = hopel3.epoch_ion_delta[tidx]*1e-3  ; half spin, same for L3 and L2.
datl3 = reform((hopel3.(where(tag_names(hopel3) eq type0)))[tidx,eidx,*])
l3pas = hopel3.pitch_angle

print, '**** HOPE/L3'
print, 'time: '+time_string(ut)
print, 'energy: '+string(en)+' eV'


datl2 = reform((hopel2.(where(tag_names(hopel2) eq type0)))[tidx,eidx,*,*])
idx = where(datl2 eq 0, cnt)
datl2 = double(datl2)
if cnt ne 0 then datl2[cnt] = !values.d_nan
smode = reform(hopel2.sector_collapse_cntr[tidx,*])
nsec = 15   ; according to brian, 15 azimuthal sectors.
nsec = 16

print, '**** HOPE/L2'
print, 'time: '+time_string(ut)
print, 'energy: '+string(en)+' eV'

; it's generally ok to assume b is constant within 1 spin.
get_data, 'b_uvw', tmp, buvw
buvw = sunitvec(reform(buvw[tidx,*]))

print, '**** HOPE/EMFISIS'
print, 'b_uvw: '+strjoin(string(buvw),', ')


; the pixel's look dir in the fixed uvw coord.
; in this coord, background b is fixed, polar pixels rotate around spin (w).
dangb = 360/nsec
angas = 90-(findgen(npxl)*36+18)        ; polar angle.
angbs = (findgen(nsec)-nsec/2)*dangb    ; azimuth angle.
angas*= rad
angbs*= rad

plot, [0,180], [0,6], xstyle = 1, ystyle = 1, /nodata, $
    xtitle = 'Pitch angle (deg)', ytitle = 'Pixel #', $
    title = 'B_hat: '+strjoin(string(buvw),',')+' at '+time_string(ut)
l2pas = dblarr(nsec,npxl)
for i = 0, npxl-1 do begin
    print, 'pixel: '+string(i+1)
    for j = 0, nsec-1 do begin
        puvw = [cos(angas[i])*[cos(angbs[j]),sin(angbs[j])],sin(angas[i])]
        l2pas[j,i] = sang(puvw, buvw, /deg)
        print, l2pas[j,i]
    endfor
    plots, l2pas[*,i], dblarr(nsec)+i+1, psym = 1
endfor
stop

device, decomposed = 0
loadct2, 43
plot, [0,180],[2,2e6], /nodata, xstyle = 1, ystyle = 1, ylog = 1
oplot, l3pas, datl3, psym = 4, color = 6
for i = 0, npxl-1 do begin
    ; combine the collapsed sectors, use the meam pitch angle.
    ty = datl2[0:nsec-1:smode[i],i] > 3
    nval = n_elements(ty)
    tx = fltarr(nval)
    for j = 0, nval-1 do tx[j] = mean(l2pas[j*smode[i]:(j+1)*smode[i]-1,i])
;    print, nval
;    tx = l2pas[*,i]
;    ty = datl2[0:nsec-1,i]
    oplot, tx, ty, psym = 1, color = i
endfor

datl2 = reform(datl2,npxl*nsec)
l2pas = reform(l2pas,npxl*nsec)
; bin the L2 pitch angles and take the mean value.
unidat = dblarr(nl3pa0)
for i = 0, nl3pa0-1 do begin
    idx = where(abs(l2pas-l3pa0s[i]) le dpa0s[i], cnt)
    if cnt eq 0 then unidat[i] = !values.d_nan else unidat[i] = mean(datl2[idx])
endfor
oplot, l3pa0s, unidat, psym = 5, color = 6

end
