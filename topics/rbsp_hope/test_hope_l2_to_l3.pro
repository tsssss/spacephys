;+
; Type: crib.
; Purpose: Try to reproduce L3 data from L2 data, make plot to compare
;   the pitch angle spectrogram and energy spectrogram of L2 and L3 results.
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
date = '2013-05-01'
ut0 = time_double(date+'/07:40')
utr = time_double(date)+[0,86400d]
probe = ['b']
ion = 'proton'

case ion of
    'proton': type0 = 'FPDU'
    'oxygen': type0 = 'FODU'
    'helium': type0 = 'FHEDU'
endcase


; suppress math exception.
!except = 0

; constants.
npixel = 5
nenergy = 72
nsection = 16

rad = !dpi/180d
deg = 180d/!dpi

timespan, date, 1, /day

; L3 pitch angles.
l3pas = [4.5,18,36,54,72,90,108,126,144,162,175.5]
nl3pa = n_elements(l3pas)


; **** load data.
load = 0
vars = ['hopel2','hopel3','hopel2_en','hopel3_en','hopel2_pa','hopel3_pa','b_uvw']
nvar = n_elements(vars)

for i = 0, nvar-1 do begin
    ; no data.
    get_data, vars[i], tmp, dat
    if size(dat,/type) eq 8 then continue
    if n_elements(dat) le 1 then load = 1
    ; no time overlap.
    if n_elements(tmp) gt 2 then begin
        if max(utr) lt min(tmp) or min(utr) gt max(tmp) then load = 1
    endif
endfor


if load eq 1 then begin
    
    ; load spice kernal.
    rbsp_load_spice_kernels


    ; load b field.
    emfisis = sread_rbsp_emfisis(utr)

    b_uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    b_gse = emfisis.mag

    store_data, 'b_gse', b_uts, b_gse
    tvar = 'b_gse'
    options, tvar, 'labels', ['x','y','z']

    rbsp_uvw2gse, 'b_gse', newname = 'b_uvw', probe = prob, /inverse

    tvar = 'b_uvw'
    options, tvar, 'labels', ['u','v','w']
    get_data, tvar, t0, dat
    store_data, tvar+'_mag', t0, snorm(dat), limits = {ynozero:1}
    options, tvar, 'ytitle', 'B!IUVW!N!C(nT)'


    ; load hope l2 data.
    hopel2 = sread_rbsp_hope(utr, probes = probe)
    store_data, 'hopel2', 0, hopel2


    ; load l3 data.
    hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
    store_data, 'hopel3', 0, hopel3

    i0 = where(tag_names(hopel3) eq type0, cnt)
    if cnt eq 0 then message, 'no such ion species ...'
    datl3 = hopel3.(i0)
    uts = sfmepoch(hopel3.epoch_ion,'unix')
    idx = where(datl3 eq -1e31, cnt)
    if cnt ne 0 then datl3[idx] = 0;!values.d_nan

    dat = total(datl3,3,/nan)
    store_data, 'hopel3_en', uts, dat, hopel3.hope_energy_ion

    dat = total(datl3,2,/nan)
    store_data, 'hopel3_pa', uts, dat, hopel3.pitch_angle

endif


tplot_options, 'labflag', -1
tplot_options, 'xmargin', [20,15]
tplot_options, 'constant', 0


; load the l2 and l3 data.
get_data, 'hopel2', tmp, hopel2
get_data, 'hopel3', tmp, hopel3


; calculate pixel looking angle.
t0 = sfmepoch(hopel2.epoch,'unix')
ion_t0 = sfmepoch(hopel2.epoch_ion,'unix')
ion_dt = hopel2.epoch_ion_delta*1d-3
nrec = n_elements(ion_t0)
smode = hopel2.sector_collapse_cntr


; the spectrogram.

i0 = where(tag_names(hopel2) eq type0, cnt)
if cnt eq 0 then message, 'no such ion species ...'
datl2 = hopel2.(i0)                     ; raw l2 data.

i0 = where(tag_names(hopel3) eq type0, cnt)
if cnt eq 0 then message, 'no such ion species ...'
datl3 = hopel3.(i0)                     ; raw l3 data.
datl4 = dblarr(nrec,nenergy,nl3pa)      ; our version of l3 data.

get_data, 'b_uvw', b_t0, buvw

tmp = min(ion_t0-ut0, idx, /absolute)
for i = idx, nrec-1 do begin
    
    ; the values at this record.
    tdt = ion_dt[i]                     ; half integration time.
    tt0 = ion_t0[i]-tdt                 ; the start of a spin.
    tdatl2 = reform(datl2[i,*,*,*])     ; [72,16,5].
    tsmode = reform(smode[where(t0 eq tt0),*])  ; [5], default [4,2,1,2,4].
;    tsmode = fltarr(npixel)+1          ; [1,1,1,1,1].
    nsecs = 16/tsmode                   ; # of section.

    nl2pa = total(nsecs)                ; # of pitch angles.
    l2pas = fltarr(nl2pa)               ; the pitch angle.

    l2spec = fltarr(nenergy,nl2pa)      ; the raw pa spectrogram.
    l4spec = fltarr(nenergy,nl3pa)      ; the output pa spectrogram.
    
    l3pas = hopel3.pitch_angle          ; [11].
    l3spec = reform(hopel3.fpdu[i,*,*]) ; [72,11].

    angas = findgen(npixel)*36+18       ; polar angle.
;   agnbs = 0       ; pixels are at 0 azimuthal angle in uvw coord.

    ; calc the pitch angle for each pixel and sector.
    for j = 0, npixel-1 do begin

        nsec = nsecs[j]
        secidx = fix(smkarthm(0,16/nsec,nsec,'x0'))
        
        tt1 = smkarthm(tt0,tdt*2/nsec,nsec,'x0')
        tdat = tdatl2[*,secidx,j]       ; [72,nsec].
        tang = angas[j]*rad
        tb = sinterpol(buvw, b_t0, tt1, /spline)
        
        ; calc each sector.
        for k = 0, nsec-1 do begin
            vec1 = reform([sin(tang),0,cos(tang)])
            vec2 = reform(sunitvec(tb[k,*]))
            idx = ((j eq 0)?0:total(nsecs[0:j-1]))+k
            l2pas[idx] = sang(vec1,vec2)*deg
            l2spec[*,idx] = tdat[*,k]
        endfor
    endfor



; for each energy bin, compare the l2 and l3 pitch angle data.
    for j = 0, nenergy-1 do begin

        ty0 = reform(l2spec[j,*])
        tx0 = reform(l2pas)

        idx = sort(tx0)
        ty1 = ty0[idx]
        tx1 = tx0[idx]

        idx = uniq(tx1)
        tx2 = []
        ty2 = []

        for k = 0, n_elements(tx1)-1 do begin
            l = k+1
            while (l lt n_elements(tx1)-1) and ((tx1[l]-tx1[k]) le 1) do l++
            if l lt n_elements(tx1)-1 then l = l-1
            if l lt k then l = k
            tx2 = [tx2,mean(tx1[k:l])]
            ty2 = [ty2,mean(ty1[k:l])]
            k = l
        endfor
        
        idx = where(ty2 ne 0, cnt)
        if cnt eq 0 then continue
        
        ty2 = ty2[idx]
        tx2 = tx2[idx]


        ; method 1.
;        l4spec[j,*] = interpol(ty2,tx2,l3pas)

        ; method 2.
;        l4spec[j,*] = interpol(ty2,cos(tx2*rad),cos(l3pas*rad))

        ; method 3.        
        order = 3
        cfit = svdfit(cos(tx2*rad), ty2, order, /legendre, yfit=ty3)
        for k = 0, nl3pa-1 do l4spec[j,k] = total(svdleg(cos(l2pas[k]*rad),order)*cfit)
      
        ; remove extropolation points.
        p0 = 18
        for k = 0, nl3pa-1 do begin
            idx = where(tx2 ge (l3pas[k]-p0) and tx2 le (l3pas[k]+p0), cnt)
            if cnt eq 0 then l4spec[j,k] = 0
        endfor

        sgtruecolor
        plot, tx2, ty2, psym = 4, $
            background = sgcolor('white'), color = sgcolor('black'), $
            yrange = [1e2,1e12], ystyle = 1, ytitle = 'Energy (eV)', ylog = 1, $
            xrange = [0,180], xstyle = 1, xtitle = 'pitch angle (deg)'

        oplot, l3pas, l4spec[j,*], psym = 2, color = sgcolor('blue')
        oplot, l3pas, l3spec[j,*], psym = 2, color = sgcolor('red')

        stop

        
    endfor

    
    idx = where(l4spec le 0, cnt)
    if cnt ne 0 then l4spec[idx] = 0
    datl4[i,*,*] = l4spec
    
endfor


c0 = 1d/1.27
c0 = 1d
datl4 *= c0

store_data, 'hopel2_pa', ion_t0, total(datl4,2,/nan), l3pas
store_data, 'hopel2_en', ion_t0, total(datl4,3,/nan), hopel2.hope_energy_ion


vars = ['hopel2_en','hopel3_en']
options, vars, 'ylog', 1
options, vars, 'yrange', [1,5e4]

vars = ['hopel2_pa','hopel3_pa']
options, vars, 'ylog', 0
options, vars, 'yrange', [0,180]


vars = ['hopel2_en','hopel2_pa','hopel3_en','hopel3_pa']
options, vars, 'spec', 1
options, vars, 'no_interp', 1
options, vars, 'zlog', 1
options, vars, 'zrange', [1e2,1e12]
options, vars, 'ystyle', 1
options, vars, 'ztitle', 's!E-1!Ncm!E-2!Nsr!E-1!NeV!E-1!N'

options, 'hopel2_en', 'ytitle', 'Hope L2!CEnergy!C(eV)'
options, 'hopel2_pa', 'ytitle', 'Hope L2!CPitch Angle!C(Deg)'
options, 'hopel3_en', 'ytitle', 'Hope L3!CEnergy!C(eV)'
options, 'hopel3_pa', 'ytitle', 'Hope L3!CPitch Angle!C(Deg)'


; plot 1.

titl = 'RBSP-'+strupcase(probe)+' Proton'
vars = ['hopel2_en','hopel3_en','hopel2_pa','hopel3_pa']

sgindexcolor, 43
tplot, vars, title = titl;, trange = utr

end
