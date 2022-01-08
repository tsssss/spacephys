
idinfos = psbl_de_id('all_e')
nid = n_elements(idinfos)

utr = ['2012-09-25','2016-01-01']
;utr = ['2013-01-01','2016-01-01']
utr = time_double(utr)



rad = !dpi/180
deg = 180/!dpi
re = 6378d & re1 = 1d/re

device, decomposed = 0
loadct2, 43


; **** plot: event versus dst.
symhfn = shomedir()+'/Google Drive/works/data/rbsp_de/symh.tplot'
if file_test(symhfn) eq 1 then tplot_restore, filename = symhfn

mltfn = shomedir()+'/Google Drive/works/data/rbsp_de/apo_mlt.tplot'
if file_test(mltfn) eq 1 then tplot_restore, filename = mltfn


aefn = shomedir()+'/Google Drive/works/data/rbsp_de/ae.tplot'
if file_test(aefn) eq 1 then tplot_restore, filename = aefn


tvar = 'ae'
if tnames(tvar) eq '' then begin
    dat = sread_omni(utr, vars=['Epoch','AE_INDEX'])
    uts = sfmepoch(dat.epoch,'unix')
    store_data, tvar, uts, dat.ae_index, limits = $
        {constant:[0,500,1000], ytitle:'AE!C(nT)'}
    tplot_save, tvar, filename = aefn
endif



get_data, 'ae', uts, ae
get_data, 'symh', uts, dst

maxdt0 = 10*60    ; min.

maxde = idinfos.emax
dsts = dblarr(nid)
aes = dblarr(nid)
autr = []
butr = []
adat = []
bdat = []

for i = 0, nid-1 do begin
    tprobe = strmid(idinfos[i].id,0,1,/reverse)
    tutr = idinfos[i].utr
    if tprobe eq 'a' then begin
        adat = [adat,!values.d_nan, 0,2, !values.d_nan]
        autr = [autr, tutr[0]+[-1,0],tutr[1]+[0,1]]
    endif else begin
        bdat = [bdat,!values.d_nan, 0,2, !values.d_nan]
        butr = [butr, tutr[0]+[-1,0],tutr[1]+[0,1]]
    endelse

    idx = where(uts ge tutr[0]-maxdt0 and uts le tutr[1]+maxdt0, cnt)
    if cnt eq 0 then dsts[i] = !values.d_nan else dsts[i] = min(dst[idx])
    if cnt eq 0 then aes[i] = !values.d_nan else aes[i] = max(ae[idx])
endfor

store_data, 'rbspa_all_e', autr, adat, limits = {ytitle:'', ytickformat:'(A1)', yticks:1, yminor:0, colors:6, panel_size:0.3, yrange:[0,2], labels:'RBSP-A'}
store_data, 'rbspb_all_e', butr, bdat, limits = {ytitle:'', ytickformat:'(A1)', yticks:1, yminor:0, colors:2, panel_size:0.3, yrange:[0,2], labels:'RBSP-B'}


poss = sgcalcpos(2)
binsize = 90
binsize = 1
xr = [0,1000]
y1s = histogram(ae, binsize = binsize, locations = xs, min = 0, max = 2000)
plot, xs, y1s, psym = 10, ylog = 0, ystyle = 1, xrange = xr, position = poss[*,0]

y2s = histogram(aes, binsize = binsize, locations = xs, min = 0, max = 2000)
y2s = smooth(y2s, 5)
oplot, xs, y2s*1e4, psym = 10, color = 6

plot, xs, double(y2s)/y1s, position = poss[*,1], xrange = xr, /noerase, psym = -1


poss = sgcalcpos(2)
binsize = 5
binsize = 1
xr = [0,100]
y1s = histogram(-dst, binsize = binsize, locations = xs, min = 0, max = 200)
plot, xs, y1s, psym = 10, ylog = 0, ystyle = 1, xrange = xr, position = poss[*,0]

y2s = histogram(-dsts, binsize = binsize, locations = xs, min = 0, max = 200)
y2s = smooth(y2s, 5)
oplot, xs, y2s*1e3, psym = 10, color = 6

plot, xs, double(y2s)/y1s, position = poss[*,1], xrange = xr, /noerase, psym = -1



stop

ofn = 0
ofn = shomedir()+'/fig_all_e_dst_ae_overview.pdf'
sgopen, ofn, xsize = 50, ysize = 5, /inch

device, decomposed = 0
loadct2, 43

tplot_options, 'ystyle', 1
tplot_options, 'labflag', -1
tplot_options, 'version', 3


tvar = 'rb_apo_mlt'
options, tvar, 'colors', [6,2]

tvar = 'ae'
options, tvar, 'yrange', [0,1500]
options, tvar, 'yticks', 3
options, tvar, 'yminor', 5

tvar = 'symh'
options, tvar, 'yrange', [50,-100]
options, tvar, 'yticks', 3
options, tvar, 'yminor', 5

vars = ['symh', 'ae', 'rbspa_all_e', 'rbspb_all_e', 'rb_apo_mlt']
options, vars, 'yticklen', 0.001
options, vars, 'xticklen', 1
options, vars, 'xgridstyle', 1

tplot, vars, trange = minmax(uts)

sgclose
stop


; **** plot: event in MLT-L plane.

n0 = 50     ; # of azimuthal points for the circles.
xr =-7.1*[-1,1]
yr = 7.1*[-1,1]
xtickv = [-6,-4,-2,0,2,4,6]
ytickv = [-6,-4,-2,0,2,4,6]
tpos = [.15,.15,.85,.85]



ofn = 0
;ofn = shomedir()+'/fig_all_e_mlt_lshell_distr.pdf'
sgopen, ofn, xsize = 5, ysize = 5, /inch


plot, xr, yr, /polar, /nodata, $
    xrange = xr, yrange = yr, /isotropic, xstyle = 1, ystyle = 1, $
    xtitle = 'Re', ytitle = 'Re', xtickv = xtickv, ytickv = ytickv, $
    xtickname = string(abs(xtickv),format='(I0)'), $
    ytickname = string(abs(ytickv),format='(I0)'), $
    position = tpos, symsize = 0.8, $
    title = 'RBSP 32 Hz |E|>25 mV/m, >30 sec, '+sgnum2str(nid)+' events!CMLT-L plane'
oplot, xr, [0,0], linestyle = 2
oplot, [0,0], yr, linestyle = 2

; add earth.
tmp = smkarthm(0,2*!dpi,n0,'n')
txs = cos(tmp)
tys = sin(tmp)
polyfill, txs, tys, /data, color = sgcolor('white')
idx = where(txs le 0)
txs = [0,txs[idx],0]
tys = [1,tys[idx],-1]
polyfill, txs, tys, /data, color = sgcolor('black')

oplot, smkarthm(1,1,n0,'n'), smkarthm(0,2*!dpi,n0,'n'), /polar
oplot, smkarthm(3,3,n0,'n'), smkarthm(0,2*!dpi,n0,'n'), /polar, linestyle = 2
oplot, smkarthm(7,7,n0,'n'), smkarthm(0,2*!dpi,n0,'n'), /polar, linestyle = 2


;oplot, infos.lshell, (180-infos.mlt*15)*rad, psym = 1, symsize = 0.8, /polar


for i = 0, nid-1 do begin
    tpos = idinfos[i].posgsm
    tutr = idinfos[i].utr
    tprobe = idinfos[i].probe
    
    lshell = idinfos[i].lshell
    mlt = idinfos[i].mlt
    
    tx =-lshell*cos(mlt*15*rad)
    ty = lshell*sin(mlt*15*rad)
    
    plots, tx, ty, psym = 1, symsize = 0.5
endfor

sgclose



; **** plot: event in ILat-MLT plane.


mlatvs = [80,70,60]
minlat = 60
mltvs = [0,6,12,18]



r0 = 1+100*re1
n0 = 50     ; # of azimuthal points for the circles.
xr = [-1,1]
yr = [-1,1]
sdeg = '!9'+string(176b)+'!X'


; convert to polar coord, and normalize.
rvs = (90-mlatvs)/(90d - 60)
tvs = mltvs*15
tvs = tvs*rad


ofn = 0
ofn = shomedir()+'/fig_all_e_ilat_mlt_distr.pdf'
sgopen, ofn, xsize = 5, ysize = 5, /inch


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


for i = 0, nid-1 do begin
    tpos = idinfos[i].posgsm
    tutr = idinfos[i].utr
    tprobe = idinfos[i].probe
    
;    posgsm = idinfos[i].posgsm
;    tut = mean(tutr)
;    tet = stoepoch(tut,'unix')
;    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
;    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
;    
;    xp = posgsm[0] & yp = posgsm[1] & zp = posgsm[2]
;    if n_elements(dir) eq 0 then $
;        dir = (zp gt 0)? -1: 1      ; zp gt 0 > northern hemisphere.
;        
;    ; model parameters, t89.
;    par = 2
;    
;    geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
;        /refine, /ionosphere, /t89
;    
;    ; convert from gsm to mag.
;    geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
;    
;    fptmlat = asin(tzf/r0)*deg
;    fptmlon = atan(tyf,txf)*deg
;    fptmlt = slon2lt(fptmlon/15, tet, /mag)

    lshell = idinfos[i].lshell
    ilat = acos(sqrt(1d/lshell))*deg
    mlt = idinfos[i].mlt
    
    rs = (90-abs(ilat))/(90d - 60)
    tx = rs*cos(mlt*15*rad)
    ty = rs*sin(mlt*15*rad)
    
    plots, tx, ty, psym = 1, symsize = 0.5
endfor

sgclose



end