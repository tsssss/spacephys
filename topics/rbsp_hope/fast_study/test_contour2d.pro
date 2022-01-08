;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

; date and time.
date = '2013-04-14'
date = '2014-02-19'
date = '2013-05-01'
ut0 = time_double(date+'/07:41')
utr = time_double(date)+[0,86400d]
probe = ['b']
ion = 'proton'

case ion of
    'proton':zrng = [3.5,6]
    'oxygen':zrng = [3,5]
    'helium':zrng = [2,5]
endcase

; suppress math exception.
!except = 0

; constants.
npxl = 5

rad = !dpi/180d
deg = 180d/!dpi

timespan, date, 1, /day

case ion of
    'proton': type0 = 'FPDU'
    'oxygen': type0 = 'FODU'
    'helium': type0 = 'FHEDU'
endcase

; load hope l3 data.
hopel3 = sread_rbsp_hope_l3(utr, probes = probe)

case ion of
    'proton': mass0 = 1
    'oxygen': mass0 = 16
    'helium': mass0 = 4
endcase
mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.


tmp = min(hopel3.epoch_ion-stoepoch(ut0,'unix'),rec,/absolute)
ut0 = sfmepoch(hopel3.epoch_ion[rec],'unix')
i0 = where(tag_names(hopel3) eq type0)
ion_dt = 12

; level 3.
datl3 = reform((hopel3.(i0))[rec,*,*])            ; in [nen,npa].
enl3s = reform(hopel3.hope_energy_ion[rec,*])
pal3s = reform(hopel3.pitch_angle)
npal3 = n_elements(pal3s)
idx = where(datl3 eq -1e31, cnt)
if cnt ne 0 then datl3[idx] = !values.d_nan

; remove duplicated energy bins.
enl3s = enl3s[0:*:2]
nenl3 = n_elements(enl3s)
datl3 = datl3[0:*:2,*]


; the data for polar contour, [nen,2*npa].
tdat = [[reverse(datl3,2)],[datl3]]
tang = [-reverse(pal3s),pal3s]
tdis = sqrt(2*enl3s/mass0)*1e-3

tang = tang ## ((bytarr(nenl3)+1)+smkarthm(0,0.001,nenl3,'n'))
tdis = tdis # (bytarr(2*npal3)+1)

tang = tang*rad

idx = where(finite(tdat,/nan))
tdat[idx] = 0


idx = where(tdat ne 0)
min0 = min(tdat[idx],/nan)
max0 = max(tdat[idx],/nan)
nztick = 10
if n_elements(zrng) eq 0 then zrng = [floor(alog10(min0)),ceil(alog10(max0))-2]
tpos = [0.15,0.15,0.85,0.85]
titl = 'RBSP-'+strupcase(probe)+' HOPE L3 '+ion+' Eflux!C'+$
    time_string(ut0)+' - '+time_string(ut0+ion_dt,tformat='hh:mm:ss')


wid = shomedir()+'/hope_l3_pitch2d_'+ion+'.eps'
;wid = 0
sgpsopen, wid, xsize = 5, ysize = 5, /inch
sgindexcolor, 43
sgdistr2d, tdat, tang, tdis, position = tpos, zrange = zrng, title = titl, ncolor = 10
sgpsclose, /pdf


end
