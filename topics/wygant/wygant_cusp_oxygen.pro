pro numspec2density, vname, mass, newname = newname
    if n_elements(mass) eq 0 then mass = 1  ; proton.
    get_data, vname, t0, val, en    ; val in #/cm^2-s-sr-keV.
    val*= 4*!dpi                    ; val in #/cm^2-s-keV.
    den = en[1:*]-en[0:n_elements(en)-2] & den = [0,den]
    for j = 0, n_elements(t0)-1 do $; en in keV.
        val[j,*]*= den              ; val in #/cm^2-s.
    cc = sqrt(2*1.6/1.67)*1e6       ; convert sqrt(eV/mp) to cm/s.
    cc1 = 1d/cc
    for j = 0, n_elements(t0)-1 do $; mass in mp.
        val[j,*]*= sqrt(mass/(en*1e3))*cc1 ; val in #/cm^3.
    den = total(val,2)
    if n_elements(newname) eq 0 then newname = vname+'_density'
    store_data, newname, t0, den
end

pro wygant_cusp_oxygen

ut0 = systime(1)
ofn = shomedir()+'/cusp_polar_ion_'+time_string(ut0,tformat='vYYYY_MMDD')+'.pdf'
sgopen, ofn, xsize = 10, ysize = 8, /inch
erase

device, decomposed = 0
loadct2, 43
posl = [0.1,0.15,0.4,0.9]
posr = [0.6,0.15,0.9,0.9]

!p.font = 1 & !y.charsize = 1 & !x.charsize = 0.9 &!z.charsize = 0.8
tplot_options, 'ygap', 0.25
tplot_options, 'ynozero', 1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
tplot_options, 'labflag', 1
time_stamp, /off

tr = time_double(['2004-11-07/19:00','2004-11-07/21:00'])
tr = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
z0 = 1e3 & z1 = 1e7
enz0 = 1e3 & enz1 = 1e7
dpi = !dpi
polabs = 'po_'+['ilat','mlt','dis']
massh = 1d & masso = 16d & masshe = 4d

; **** load polar pos.
tmp = sread_polar_orbit(tr)
t0 = sfmepoch(tmp.epoch,'unix')
tval = sqrt(total(tmp.pos_gse^2,2))*(1d/6378.137d)
store_data, 'po_dis', t0, tval, limits = {ytitle:'Dist (Re)'}
store_data, 'po_mlt', t0, tmp.mlt, limits = {ytitle:'MLT (hr)'}
tval = acos(sqrt(1D/tmp.lshell))*(180d/!dpi)   ; in degree.
;idx = where(tmp.mag_latitude lt 0, cnt)
;if cnt gt 0 then tval[idx] = -tval[idx]
store_data, 'po_ilat', t0, tval, limits = {ytitle:'ILat (deg)'}

; **** load polar timas.
tmp = sread_polar_timas(tr, type = 'k1')
opt = {spec:1,ylog:1,zlog:1, xstyle:1, ystyle:1,no_interp:1}
uts = sfmepoch(tmp.epoch,'unix')
tmp.energy*=1e-3    ; eV->keV.
store_data, 'fh', uts, tmp.flux_h, tmp.energy, limits = opt
store_data, 'fo', uts, tmp.flux_o, tmp.energy, limits = opt
store_data, 'fhe1', uts, tmp.flux_he_1, tmp.energy, limits = opt
store_data, 'fhe2', uts, tmp.flux_he_2, tmp.energy, limits = opt
options, 'fh', 'ytitle', 'number flux!CH+!C(keV)'
options, 'fo', 'ytitle', 'number flux!CO+!C(keV)'
options, 'fhe1', 'ytitle', 'number flux!CHe+!C(keV)'
options, 'fhe2', 'ytitle', 'number flux!CHe++!C(keV)'
options, 'fh', 'ztitle', 'Log Nflux!C(#/cm!U2!N-s-sr-keV)'
options, 'fo', 'ztitle', 'Log Nflux!C(#/cm!U2!N-s-sr-keV)'
options, 'fhe1', 'ztitle', 'Log Nflux!C(#/cm!U2!N-s-sr-keV)'
options, 'fhe2', 'ztitle', 'Log Nflux!C(#/cm!U2!N-s-sr-keV)'
store_data, 'mh', uts, tmp.flux_h*sqrt(massh), tmp.energy, limits = opt
store_data, 'mo', uts, tmp.flux_o*sqrt(masso), tmp.energy, limits = opt
store_data, 'mhe1', uts, tmp.flux_he_1*sqrt(masshe), tmp.energy, limits = opt
store_data, 'mhe2', uts, tmp.flux_he_2*sqrt(masshe), tmp.energy, limits = opt
options, 'mh', 'ytitle', 'mass flux!CH+!C(keV)'
options, 'mo', 'ytitle', 'mass flux!CO+!C(keV)'
options, 'mhe1', 'ytitle', 'mass flux!CHe+!C(keV)'
options, 'mhe2', 'ytitle', 'mass flux!CHe++!C(keV)'
options, 'mh', 'ztitle', 'Log Mflux!C(kg!U1/2!N/cm!U2!N-s-sr-keV)'
options, 'mo', 'ztitle', 'Log Mflux!C(kg!U1/2!N/cm!U2!N-s-sr-keV)'
options, 'mhe1', 'ztitle', 'Log Mflux!C(kg!U1/2!N/cm!U2!N-s-sr-keV)'
options, 'mhe2', 'ztitle', 'Log Mflux!C(kg!U1/2!N/cm!U2!N-s-sr-keV)'

; **** load polar hydra.
tmp = sread_polar_hydra(tr, type = 'k0')
store_data, 'ne', sfmepoch(tmp.epoch,'unix'), tmp.nele, $
    limits = {ytitle:'number density!C(cm!U-3!N)',labels:'Hydra ele'}

;tplot_restore, filename = shomedir()+'/wygant_cusp_oxygen_2004_1107.tplot'

; for flux spectrogram.
vars = ['fh','fo','fhe1','fhe2']
titl = ''
z0 = 1e4 & z1 = 1e8
zlim, vars, z0, z1, 1

; for flux spectrogram*sqrt(mass).
vars = ['mh','mo','mhe1','mhe2']
titl = ''
z0 = 1e4 & z1 = 1e8
zlim, vars, z0, z1, 1

; **** integrate timas spectrogram to get density.
; and load hydra electron density.
; flux = n*v = n*sqrt(2E/m), we want to see n, so we need to multiply flux by sqrt(m).
numspec2density, 'fh', massh, newname = 'nh'
numspec2density, 'fo', masso, newname = 'no'
options, 'nh', 'ytitle', 'density!C(cm!U-3!N)'
options, 'nh', 'labels', 'Timas H+'
options, 'no', 'ytitle', 'density!C(cm!U-3!N)'
options, 'no', 'labels', 'Timas O+'

; **** calc the ratios.
; number flux spec ratio.
get_data, 'fh', t0, fh, en, limits = lims
get_data, 'fo', t0, fo, en
no2hspec = double(fo)/fh
;idx = where(fh le z0)
;no2hspec[idx] = 0
nh2ospec = double(fh)/fo
store_data, 'nh2ospec', t0, nh2ospec, en, limits = lims
options, 'nh2ospec', 'ytitle', 'number flux!Cratio H+/O+'
zlim, 'nh2ospec', 1, 1e3, 1

; mass flux spec ratio.
get_data, 'mh', t0, mh, en, limits = lims
get_data, 'mo', t0, mo, en
mh2ospec = double(mh)/mo
store_data, 'mh2ospec', t0, mh2ospec, en, limits = lims
options, 'mh2ospec', 'ytitle', 'mass flux!Cratio H+/O+'
zlim, 'mh2ospec', 1, 1e3, 1

get_data, 'nh', t0, nh
get_data, 'no', t0, no
store_data, 'nh2o', t0, nh/no, limits = {ytitle:'number density!Cratio', labels:'Timas H+/O+'}
store_data, 'rhoh2o', t0, nh*massh/(no*masso), limits = {ytitle:'mass density!Cratio', labels:'Timas H+/O+'}

envars = vars+'_en'
vars = ['fh','fo','mh','mo']
nvar = n_elements(vars)
for i = 0, nvar-1 do begin
    get_data, vars[i], t0, val, en, limits = lim
    den = en[1:*]-en[0:n_elements(en)-2] & den = [0,den]
    for j = 0, n_elements(t0)-1 do val[j,*]*= den
    store_data, envars[i], t0, val, en, limits = lim
endfor
options, 'fh_en', 'ztitle', 'Log Eflux!C(keV/cm!U2!N-s-sr-keV)'
options, 'fo_en', 'ztitle', 'Log Eflux!C(keV/cm!U2!N-s-sr-keV)'
options, 'fhe1_en', 'ztitle', 'Log Eflux!C(keV/cm!U2!N-s-sr-keV)'
options, 'fhe2_en', 'ztitle', 'Log Eflux!C(keV/cm!U2!N-s-sr-keV)'
zlim, envars, enz0, enz1, 1

vars = ['ne','nh','no']
options, vars, 'ytitle', 'number density!C(cm!U-3!N)'

tvars = ['fh','fo','nh2ospec','mh','mo','mh2ospec']
nvar = n_elements(tvars)
pos = (sgcalcpos(nvar, position = posl))
tplot, tvars, var_label = polabs, trange = tr, position = pos, $
    title = titl, /noerase
    
tvars = ['ne','nh','no','nh2o','rhoh2o']
nvar = n_elements(tvars)
pos = (sgcalcpos(nvar, position = posr))
tplot, tvars, var_label = polabs, trange = tr, position = pos, $
    title = '', /noerase

title = 'Polar cusp ion composition, event 2004-11-07, plotted on '+time_string(ut0)
xyouts, 0.5, 0.91, /normal, title, alignment = 0.5, charsize = 1.25
sgclose

end