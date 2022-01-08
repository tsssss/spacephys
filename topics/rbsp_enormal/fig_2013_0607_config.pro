
_2013_0607_load_data

; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1

top = 254
maxcnt = 800
maxcnt = 500
mincnt = 50

sctrng = [3.5,-2.5]
sczrng = [1.5,4.5]

eqz0 = 0        ; equatorial plane.
eqxrng = -[5.5,3.5]
eqyrng = [3.5,-2.5]


utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
tut = time_double('2013-06-07/04:55')
tet = stoepoch(tut,'unix')


symsz = 0.5
tpos = [0.2,0.2,0.9,0.9]
xr = [0,10]
yr = [0,5]
ds = 0.1
dir = 1
tmodel = 't01'
;sgeopack_par, utr, tmodel, /delete


mlats = smkarthm(70d,80,1,'dx')
;mlats = [50,60,70,80,90,75]
nline = n_elements(mlats)
mlts = 22+dblarr(nline)

geopack_epoch, tet, year, mo, dy, hr, mi, sc, /breakdown_epoch
geopack_recalc, year, mo, dy, hr, mi, sc, /date, tilt = tilt

get_data, tmodel+'_par', uts, par
par = sinterpol(par, uts, tut)


ofn = 0
ofn = shomedir()+'/fig_config.pdf'
sgopen, ofn, xsize = 6, ysize = 3, /inch
device, decomposed = 0
loadct2, 43


xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


plot, xr, yr, /nodata, /iso, position = tpos, $
    xstyle = 1, xticks = 10, xminor = 5, xtitle = 'R at 22 MLT (Re)', $
    ystyle = 1, yticks = 5, yminor = 5, ytitle = 'Z GSM (Re)'


; add earth.
tmp = findgen(21)/20*!dpi/2
plots, cos(tmp), sin(tmp)

xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1.2, /normal, alignment = 0, $
    'Model: '+strupcase(tmodel)



for i = 0, nline-1 do begin
    tmp1 = (24-mlts[i])*15*rad  ; angle between mlt and midnight.
    tmp2 = mlats[i]*rad
    v0 = [-cos(tmp2)*cos(tmp1),cos(tmp2)*sin(tmp1),sin(tmp2)]
;    v0 = v0*r0      ; approx coord in gsm.
    
;    while 1 do begin
;        geopack_igrf_gsm, v0[0], v0[1], v0[2], bx, by, bz
;        geopack_t01, par, v0[0], v0[1], v0[2], bx1, by1, bz1, tilt = tilt
;        tmp = -sunitvec([bx+bx1,by+by1,bz+bz1])
;        v1 = v0+tmp*ds
;        tv0 = [-v0[0]*cos(tmp1)+v0[1]*sin(tmp1),v0[2]]
;        tv1 = [-v1[0]*cos(tmp1)+v1[1]*sin(tmp1),v1[2]]
;        
;        if tv1[0] le xr[0] then break
;        if tv1[0] ge xr[1] then break
;        if tv1[1] le yr[0] then break
;        if tv1[1] ge yr[1] then break
;        plots, [tv0[0],tv1[0]], [tv0[1],tv1[1]]
;;        plots, [sqrt(v0[0]^2+v0[1]^2),sqrt(v1[0]^2+v1[1]^2)], [v0[2],v1[2]]
;        v0 = v1
;    endwhile
    geopack_trace, v0[0], v0[1], v0[2], dir, par, xf, yf, zf, fline = fline, /t01
    srotate, fline, tmp1, 2
    oplot, -fline[*,0], fline[*,2]
endfor


; add RBSP.
get_data, 'rbspa_pos_gsm', uts, dat
pos = sinterpol(dat, uts, tut)
tx = -pos[0]*cos(tmp1)+pos[1]*sin(tmp1)
ty = pos[2]
plots, tx, ty, psym = 6, color = 6, symsize = symsz
tmp = convert_coord([tx,ty],/data,/to_normal)
xyouts, tmp[0], tmp[1]+ychsz*0.6, /normal, alignment = 0, 'RBSP-A', color = 6

get_data, 'rbspb_pos_gsm', uts, dat
pos = sinterpol(dat, uts, tut)
tx = -pos[0]*cos(tmp1)+pos[1]*sin(tmp1)
ty = pos[2]
plots, tx, ty, psym = 6, color = 2, symsize = symsz
tmp = convert_coord([tx,ty],/data,/to_normal)
xyouts, tmp[0], tmp[1]-ychsz*1.2, /normal, alignment = 0, 'RBSP-B', color = 2
sgclose
end
