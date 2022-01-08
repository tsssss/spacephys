;+
; Type: crib.
; Purpose: Plot HOPE pitch angle for L2 & L3 data, at each energy bin
;   at certain time.
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
;date = '2013-05-01'
ut0 = time_double(date+'/07:40')
utr = time_double(date)+[0,86400d]
probe = ['b']
ion = 'proton'


; suppress math exception.
!except = 0

; constants.
npxl = 5

rad = !dpi/180d
deg = 180d/!dpi

; set time span to the day to load data.
timespan, date, 1, /day

; type0 is for loading data.
case ion of
    'electron': type0 = 'FEDU'
    'proton': type0 = 'FPDU'
    'oxygen': type0 = 'FODU'
    'helium': type0 = 'FHEDU'
endcase


; **** load data.
load = 0

vars = ['hopel2','hopel3','b_uvw']
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
    vars = emfisis.name

    b_uts = sfmepoch(*emfisis[where(vars eq 'Epoch')].value,'unix',/tt2000)
    b_gse = *emfisis[where(vars eq 'Mag')].value

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
    hopel2 = sread_rbsp_hope_l2(utr, probes = probe)
    store_data, 'hopel2', 0, hopel2
    
    ; load l3 data.
    hopel3 = sread_rbsp_hope_l3(utr, probes = probe)
    store_data, 'hopel3', 0, hopel3

endif


tplot_options, 'labflag', -1
tplot_options, 'xmargin', [20,15]
tplot_options, 'constant', 0


; load the l2 and l3 data.
get_data, 'hopel2', tmp, hopel2
get_data, 'hopel3', tmp, hopel3

; load b data.
get_data, 'b_uvw', b_uts, b_uvw


; basic info.
ion_uts = sfmepoch(hopel2.epoch_ion,'unix')
tmp = min(ion_uts-ut0,rec,/absolute)
ut0 = ion_uts[rec]
i0 = where(tag_names(hopel2) eq type0)
ion_dt = hopel2.epoch_ion_delta[rec]*1e-3*2


; level 2 data.
datl2 = reform((hopel2.(i0))[rec,*,*,*])            ; in [72,16,5].
tmp = where(datl2 eq -1e31, cnt)
if cnt ne 0 then datl2[tmp] = !values.d_nan

enl2s = reform(hopel2.hope_energy_ion[rec,*])       ; energy bins.
enmode = hopel2.energy_collapsed[rec]               ; 0 for collapsed bins.
if enmode eq 0 then begin
    enl2s = enl2s[0:*:2]
    datl2 = datl2[0:*:2,*,*]
endif
nenl2 = n_elements(enl2s)

scmode = reform(hopel2.sector_collapse_cntr[rec,*]) ; sector collapse mode.
; scmode = intarr(npxl)+1                            ; no collapsing.
; adjust value for collapsing.
for j = 0, npxl-1 do begin
    if scmode[j] eq 1 then continue
;    datl2[*,*,j] = (datl2[*,*,j]-16384)*64+16384
;    datl2[*,*,j] = datl2[*,*,j]/scmode[j]
endfor

; calc the pitch angle for all looing directions.
vec1 = reform(sinterpol(b_uvw,b_uts,ut0-ion_dt, /spline))
nsecs = 16/scmode
npal2 = total(16/scmode)
pal2s = dblarr(npal2)
for j = 0, npxl-1 do begin
    nsec = nsecs[j]
    ang0 = 11.25*scmode[j]  ; offset angle accounts for collapsing.
    anga = (j*36+18)*rad                            ; polar angle.
    angb = smkarthm(ang0,360d/nsec,nsec,'x0')*rad   ; azimuthal angle.
    for k = 0, nsec-1 do begin
        vec2 = [sin(anga)*[cos(angb[k]),sin(angb[k])],cos(anga)]
        pal2s[total(nsecs[0:j])-nsecs[j]+k] = sang(vec1,vec2,/degree)
    endfor
endfor

nenl4 = nenl2
datl4 = dblarr(nenl4,npal2)
for i = 0, nenl4-1 do begin
    tmp = []
    for j = 0, npxl-1 do tmp = [tmp,reform(datl2[i,0:*:scmode[j],j])]
    datl4[i,*] = tmp
endfor
datl2 = datl4


; level 3 data.
i0 = where(tag_names(hopel3) eq type0)
datl3 = reform((hopel3.(i0))[rec,*,*])               ; in [72,11].
tmp = where(datl3 eq -1e31, cnt)
if cnt ne 0 then datl3[tmp] = !values.d_nan

enl3s = reform(hopel3.hope_energy_ion[rec,*])
nenl3 = n_elements(enl3s)
pal3s = reform(hopel3.pitch_angle)
npal3 = n_elements(pal3s)

if enmode eq 0 then begin
    enl3s = enl3s[0:*:2]
    nenl3 = n_elements(enl3s)
    datl3 = datl3[0:*:2,*,*]
endif


; level 4 data.
pal4s = pal3s
npal4 = n_elements(pal4s)



; for each energy bin, compare the l2 and l3 pitch angle data.
for j = 0, nenl2-1 do begin

    ; level 2 data.
    ty2 = reform(datl2[j,*])
    tx2 = reform(pal2s)
    ; sort.
    idx = sort(tx2)
    ty2 = ty2[idx]
    tx2 = tx2[idx]
    ; uniq.
    idx = uniq(tx2)
    ty2 = ty2[idx]
    tx2 = tx2[idx]
    ; remove 0.
    idx = where(ty2 ne 0, cnt)
    if cnt eq 0 then continue
    ty2 = ty2[idx]
    tx2 = tx2[idx]


    ; level 3 data.
    ty3 = reform(datl3[j,*])
    tx3 = reform(pal3s)
    ; sort.
    idx = sort(tx3)
    ty3 = ty3[idx]
    tx3 = tx3[idx]
    ; uniq.
    idx = uniq(tx3)
    ty3 = ty3[idx]
    tx3 = tx3[idx]
    ; remove 0.
    idx = where(ty3 ne 0, cnt)
    if cnt eq 0 then continue
    ty3 = ty3[idx]
    tx3 = tx3[idx]

    ; level 4 data.
    tx4 = reform(pal4s)
    ty4 = dblarr(npal4)

;    ; method 1.
;    ty4 = interpol(ty2,tx2,pal4s)

    ; method 2.
;    ty4 = interpol(ty2,cos(tx2*rad),cos(pal4s*rad))

    ; method 3.        
    order = 4<n_elements(ty2)
    cfit = svdfit(cos(tx2*rad), ty2, order, /legendre)
    for k = 0, npal4-1 do ty4[k] = total(svdleg(cos(pal4s[k]*rad),order)*cfit)

    ; remove extropolation points.
    p0 = 18
    for k = 0, npal4-1 do begin
        idx = where(tx2 ge (pal4s[k]-p0) and tx2 le (pal4s[k]+p0), cnt)
        if cnt eq 0 then ty4[k] = 0
    endfor

    ; plot number flux vs pitch angle for each energy per record.
    sgtruecolor
    
    tpos = [0.1,0.1,0.9,0.9]
    clrs = ['red','lime_green','black','blue','orange']

    ; level 2.
    plot, tx2, ty2, /nodata, $
        position = tpos, $
        title = 'RBSP-'+strupcase(probe)+' '+ion+' Energy = '+ $
            sgnum2str(enl2s[j],nsgn=3)+ 'eV!C'+time_string(ut0), $
        background = sgcolor('white'), color = sgcolor('black'), $
        yrange = [1e2,1e12], ystyle = 1, ytitle = 'Number Flux (s!E-1!Ncm!E-2!Nsr!E-1!NeV!E-1!N)', ylog = 1, $
        xrange = [0,180], xstyle = 1, xtitle = 'pitch angle (deg)'
    ty2 = reform(datl2[j,*])
    tx2 = reform(pal2s)
    for k = 0, npxl-1 do begin
        tmp = total(nsecs[0:k])-nsecs[k]
        oplot, tx2[tmp:tmp+nsecs[k]-1], ty2[tmp:tmp+nsecs[k]-1], psym = 2, color = sgcolor(clrs[k])
    endfor
    ; level 4.
    oplot, tx4, ty4, psym = 5, color = sgcolor('pink')
    ; level 3.
    oplot, tx3, ty3, psym = 5, color = sgcolor('red')

    tpos = [0.12,0.7,0.22,0.88]
    plot, [0,1], [0,1], /nodata, /noerase, $
        position = tpos, color = sgcolor('black'), $
        xstyle = 1, xtickformat='(A1)', xticks = 1, xminor = 0, $
        ystyle = 1, ytickformat='(A1)', yticks = 1, yminor = 0
    plots, 0.15, 0.78, psym = 2, color = sgcolor('black')
    xyouts, 0.3, 0.75, 'L2 raw', color = sgcolor('black')
    plots, 0.15, 0.48, psym = 5, color = sgcolor('pink')
    xyouts, 0.3, 0.45, 'L2 interp', color = sgcolor('black')
    plots, 0.15, 0.18, psym = 5, color = sgcolor('red')
    xyouts, 0.3, 0.15, 'L3', color = sgcolor('black')

;    wait, 0.5
    stop
endfor


end
