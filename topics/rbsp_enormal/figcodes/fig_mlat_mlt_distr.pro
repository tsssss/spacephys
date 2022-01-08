
ifn = shomedir()+'/Google Drive/works/data/rbsp_de/32hz_event_info.tplot'
tplot_restore, filename = ifn
tvar = 'psbl_de_32hz_info'
get_data, tvar, tmp, infos



mlts = infos.fptmlt
mlats = abs(infos.fptmlat)
mlatvs = [80,70,60]
minlat = 60
mltvs = [0,6,12,18]


rad = !dpi/180
n0 = 50     ; # of azimuthal points for the circles.
xr = [-1,1]
yr = [-1,1]
sdeg = '!9'+string(176b)+'!X'


; convert to polar coord, and normalize.
rs = (90-mlats)/(90d - 60)
rvs = (90-mlatvs)/(90d - 60)
ts = mlts*15
ts = ts*rad
tvs = mltvs*15
tvs = tvs*rad


ofn = shomedir()+'/fig_mlat_mlt_distr.pdf'
;ofn = 0
sgopen, ofn, xsize = 4, ysize = 4, /inch
sgtruecolor

tpos = [0.15,0.15,0.85,0.85]

plot, xr, yr, /polar, /nodata, $
    xrange = xr, yrange = yr, /isotropic, xstyle = 5, ystyle = 5, $
    position = tpos, symsize = 0.8

; add axis.
oplot, xr, [0,0]
oplot, [0,0], yr

foreach tmlat, rvs do $
    oplot, smkarthm(tmlat,tmlat,n0,'n'), smkarthm(0,2*!dpi,n0,'n'), /polar

td = 5e-2   ; nudge value.
tr = 1.1
xyouts, tr*cos(tvs), tr*sin(tvs)-td, string(mltvs,format='(I2)'), /data, alignment = 0.5
tt = -120*rad
xyouts, rvs*cos(tt)+1.5*td, rvs*sin(tt), string(mlatvs,format='(I2)')+sdeg, /data, alignment = 0.5

tmp = findgen(11)*2*!dpi/10
txs = cos(tmp)
tys = sin(tmp)
usersym, txs, tys, /fill
plots, rs*cos(ts), rs*sin(ts), psym = 8, /data, symsize = 0.4

sgclose

end
