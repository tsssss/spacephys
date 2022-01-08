
_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'xgridstyle', 1
tplot_options, 'zcharsize', 0.8



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1







vars = ['b_gsm','keogram_t01','pf_fac_mat_map','vsc','e_en','e_en_mageis']
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

ofn = 0
ofn = shomedir()+'/fig_overview2.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43



poss = sgcalcpos(nvar*2+2)
pre0 = 'rbspa_'
tposs = poss[*,0:nvar-1]
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-A'
pre0 = 'rbspb_'
tposs = poss[*,nvar+1:nvar+1+nvar-1]
for i = 0, nvar-1 do tposs[[1,3],i] -= 0.03
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-B'

sgclose

stop



tutr = time_double(['2013-06-07/04:50','2013-06-07/05:15'])
tplot, pre0+['demat','dbmat'], trange = tutr

; plot all bands.
ofn = shomedir()+'/'+pre0+'pflux_eb_bands.pdf'
;    ofn = 0
sgopen, ofn, xsize = 11, ysize = 8.5, /inch

device, decomposed = 0
loadct2, 43

pos1 = [0.1,0.1,0.27,0.9]
pos2 = [0.4,0.1,0.57,0.9]
pos3 = [0.7,0.1,0.87,0.9]

suf = [pfids,'']
nvar = n_elements(suf)

tplot, pre0+'defac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos1)
tplot, pre0+'dbfac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos2)
tplot, pre0+'pffac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos3)

xyouts, 0.5, 0.92, /normal, alignment = 0.5, 'RBSP-'+strupcase(tprobe)+' bands of E, B and S in FAC'
sgclose









ofn = shomedir()+'/'+pre0+'pflux_eb_bands_scaled.pdf'
;    ofn = 0
sgopen, ofn, xsize = 11, ysize = 8.5, /inch

device, decomposed = 0
loadct2, 43

pos1 = [0.1,0.1,0.27,0.9]
pos2 = [0.4,0.1,0.57,0.9]
pos3 = [0.7,0.1,0.87,0.9]

nvar = n_elements(suf)

vars = ['de','db','pf']+'fac_mat'
foreach tvar, vars do stplot_minmax, pre0+tvar+suf

tplot, pre0+'defac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos1)
tplot, pre0+'dbfac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos2)
tplot, pre0+'pffac_mat'+suf, /noerase, position = sgcalcpos(nvar, position = pos3)

xyouts, 0.5, 0.92, /normal, alignment = 0.5, 'RBSP-'+strupcase(tprobe)+' bands of E, B and S in FAC'
sgclose



ofn = shomedir()+'/eb_mat_spec.pdf'
vars = ['demat','dbmat']
nvar = n_elements(vars)

sgopen, ofn, xsize = 8.5, ysize = 11, /inch

device, decomposed = 0
loadct2, 43

pos1 = [0.2,0.55,0.85,0.9]
pos2 = [0.2,0.1,0.85,0.45]
tplot, 'rbspa_'+vars, /noerase, position = sgcalcpos(nvar, position = pos1), title = 'RBSP-A E and B period spectrogram'
tplot, 'rbspb_'+vars, /noerase, position = sgcalcpos(nvar, position = pos2), title = 'RBSP-B E and B period spectrogram'

sgclose

stop


; photon count at footprint, 1x1 box.
foreach tprobe, probes do begin

    suf0 = '_'+model
    ; map to 100 km, northern hemisphere.
    scalc_map_coef, pre0+'pos_gse', pre0+'b0_gse', model = model, coord = 'gse', /igrf, prefix = pre0, suffix = suf0, dir = -1
    
    get_data, pre0+'fpt_mlon'+suf0, tuts, mlon
    mlt = slon2lt(mlon, stoepoch(tuts,'unix'), /mag, /deg)/15    ; in hour.
    mlt = (mlt+24) mod 24
    store_data, pre0+'fpt_mlt'+suf0, tuts, mlt, $
        limits = {ytitle:'MLT/fpt (hr)'}
        
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlon'+suf0, tuts, mlon
    get_data, pre0+'fpt_mlat'+suf0, tuts, mlat
    tvar = pre0+'fpt_lonlat'+suf0
    store_data, tvar, tuts, [[mlon],[mlat]]
    thm_asi_countrate, stoepoch(utr0,'unix'), site, tvar, del = 1, sc = 'rbsp'+tprobe
endforeach

stop





; **** load mltimg.
xsz = 800
ysz = xsz*0.5
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

utr = time_double(['2013-06-07/04:50','2013-06-07/05:05'])


if load then begin
    tvar = 'mltimg'
    
    asi = sread_thg_mlt(utr, sites, /half, type = 'asf')
    store_data, tvar, sfmepoch(asi.epoch, 'unix'), asi.mltimg
    
    asi = sread_thg_mlt(utr>time_double('2013-06-07/05:00'), sites, /half, type = 'asf')
    get_data, tvar, uts, imgs
    store_data, tvar, [uts,sfmepoch(asi.epoch,'unix')], [imgs,asi.mltimg]
    
    
    get_data, tvar, uts, imgs
    nrec = n_elements(uts)
    for i = 0, nrec-1 do begin
    
        ofn = shomedir()+'/thg_asi/thg_rb_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
        ;    ofn = 0
        sgopen, ofn, xsize = xsz, ysize = ysz
        
        device, decomposed = 0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        timg = reform(imgs[i,*,*])
        timg = congrid(timg, xsz, ysz)
        
        sgtv, timg, position = tpos
        sgset_map, position = tpos, color = white, xrange = xr
        
        xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
        
        for j = 0, n_elements(probes)-1 do begin
            pre0 = 'rbsp'+probes[j]+'_'
            get_data, pre0+'fpt_mlat', tuts, mlat
            get_data, pre0+'fpt_mlt', tuts, mlt
            mlat = interpol(mlat, tuts, uts[i], /quadratic)
            mlt = interpol(mlt, tuts, uts[i], /quadratic)
            plots, mlt*15, mlat, psym = syms[j], color = colors[j]
            xyouts, tpos[0], tpos[1]+ychsz*(j+1), /normal, 'RBSP-'+strupcase(probes[j]), color = colors[j], charsize = chsz
            
        endfor
        
        sgclose
    endfor
endif


stop



; **** check polarization of the waves.

tutr = time_double(['2013-06-07/04:56','2013-06-07/04:59'])

;tutr = time_double(['2013-06-07/04:57','2013-06-07/05:00'])
get_data, 'rbspa_de_fac', uts, defac
get_data, 'rbspa_db_fac', uts, dbfac

idx = where(uts ge tutr[0] and uts le tutr[1], cnt)
defac = defac[idx,*]
dbfac = dbfac[idx,*]

edge = 50
for i = 0, 2 do begin
    defac[*,i] = smooth(defac[*,i], edge)
    dbfac[*,i] = smooth(dbfac[*,i], edge)
endfor
defac = defac[edge:cnt-1-edge,*]
dbfac = dbfac[edge:cnt-1-edge,*]
cnt = cnt-edge*2

device, decomposed = 0
loadct2, 43

xr = [-40,40]
yr = xr
plot, xr, yr, /nodata, /isotropic, xtitle = 'p', ytitle = 'v'
for i = 1, cnt-1 do begin
    plots, defac[i-1:i,1], defac[i-1:i,2]
    plots, dbfac[i-1:i,1], dbfac[i-1:i,2], color = 6
endfor

stop

plot, xr, yr, /nodata, /isotropic, xtitle = 'b', ytitle = 'v'
for i = 1, cnt-1 do begin
    plots, defac[i-1:i,0], defac[i-1:i,2]
    plots, dbfac[i-1:i,0], dbfac[i-1:i,2], color = 6
endfor

stop


plot, xr, yr, /nodata, /isotropic, xtitle = 'b', ytitle = 'p'
for i = 1, cnt-1 do begin
    plots, defac[i-1:i,0], defac[i-1:i,1]
    plots, dbfac[i-1:i,0], dbfac[i-1:i,1], color = 6
endfor

stop


; **** cross correlation between -A and -B.
get_data, 'rbspa_vsc', uts, dat
idx = where(uts ge utr0[0] and uts le utr0[1])
auts = uts[idx]
avsc = dat[idx]

get_data, 'rbspb_vsc', uts, dat
idx = where(uts ge utr0[0] and uts le utr0[1])
buts = uts[idx]
bvsc = dat[idx]

; dt about 150 sec.
dt0 = -150
drec0 = dt0/dr0 ; roughly 150 sec.
tnrec = abs(drec0/4)
corr = dblarr(tnrec)
nrec = n_elements(avsc)
edge = abs(drec0)+tnrec
for i = 0, tnrec-1 do begin
    tmp1 = avsc[edge:nrec-edge-1]
    tmp2 = (shift(bvsc,drec0+i-(tnrec-1)/2))[edge:nrec-edge-1]
    corr[i] = total((tmp1-mean(tmp1))*(tmp2-mean(tmp2))/stddev(tmp1)/stddev(tmp2))
endfor
corr = corr/(nrec-edge*2)

tmp = max(corr,idx)
drec = drec0+idx-(tnrec-1)/2
dtvsc = drec*dr0

ofn = shomedir()+'/fig_vsc_correlation.pdf'
if test then ofn = 0
sgopen, ofn, xsize = 8, ysize = 4, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

xr = [200,2000]
yr = [-50,-0]
tpos = [0.1,0.2, 0.5,0.9]
dvsc = -20  ; shift -B by -20 V.
tuts = auts[edge:nrec-edge-1]
txs = tuts-tuts[0]
tys = avsc[edge:nrec-edge-1]
plot, txs, tys, xstyle = 1, ystyle = 9, yrange = yr, $
    ytitle = 'RBSP-A!CVsc (V)', position = tpos, $
    xrange = xr, xticks = 3, xminor = 6, xticklen = 1, xgridstyle = 1, $
    xtitle = 'Second from '+time_string(auts[edge])
axis, max(xr), /yaxis, yrange = yr-dvsc, ytitle = 'RBSP-B!CVsc (V)!Cshifted -20 V'
tys = (shift(bvsc, drec))[edge:nrec-edge-1]+dvsc
oplot, txs, tys, color = 6

xyouts, tpos[0]+xchsz*2, tpos[3]-ychsz*0.5-ychsz*1.2*1, /normal, 'RBSP-A'
xyouts, tpos[0]+xchsz*2, tpos[3]-ychsz*0.5-ychsz*1.2*2, /normal, 'RBSP-B, dT = '+sgnum2str(drec*dr0,nsgn=4)+' sec', color = 6


tpos = [0.7,0.2,0.95,0.9]
tx = (drec0+findgen(tnrec)-0.5*(tnrec-1))*dr0
ty = corr
plot, tx, ty, position = tpos, /noerase, /ynozero, $
    ytitle = 'Correlation', xtitle = 'Time shift (sec)'
    
xr = !y.crange
yr = !y.crange
plots, drec*dr0+[0,0], yr, linestyle = 1
tx = drec*dr0
ty = max(corr)

xyouts, tpos[0]+xchsz*10, tpos[1]+ychsz*3, /normal, 'dT = '+sgnum2str(drec*dr0,nsgn=4)+' sec'
sgclose


get_data, 'rbspa_pos_gse', uts, apos
get_data, 'rbspb_pos_gse', uts, bpos
ddis = snorm(apos-bpos)
tut = time_double('2013-06-07/04:54')
tdis = interpol(ddis, uts, tut, /quadratic)

print, sgnum2str(tdis/abs(dtvsc)*re)+' km/s'

stop



; **** load pitch 2d.
types = ['electron','proton']
log = 1
unit = 'energy'

uts = smkarthm(utr[0],utr[1],12,'dx')

foreach tprobe, probes do begin
    hopel3 = sread_rbsp_hope_l3(utr, probes = tprobe)
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
        plot_hope_l3_pitch2d, uts[i], types[j], unit = unit, $
        log = log, hopel3 = hopel3, probe = tprobe
endforeach

stop

; **** load pina asf.
tvar = site+'_asf'
if tnames(tvar) eq '' then begin
    asi = sread_thg_asi(utr, site, type = 'asf')
    store_data, tvar, asi.utsec, asi.img
endif


xsz = 800
ysz = xsz
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

get_data, tvar, uts, imgs
nrec = n_elements(uts)
for i = 0, nrec-1 do begin

    ofn = shomedir()+'/thg_asi/'+site+'/thg_asf_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
    ;    ofn = 0
    sgopen, ofn, xsize = xsz, ysize = ysz
    
    device, decomposed = 0
    loadct, 39
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    timg = reform(imgs[i,*,*])
    timg *= 64d/(median(timg) > 1)
    timg = bytscl(timg)
    timg = congrid(timg, xsz, ysz)
    
    sgtv, timg, position = tpos
    
    xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
    
    sgclose
endfor

stop


; **** load mltimg.
xsz = 800
ysz = xsz*0.5
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

if load then begin
    tvar = 'mltimg'
    
    asi = sread_thg_mlt(utr, sites, /half, type = 'asf')
    store_data, tvar, sfmepoch(asi.epoch, 'unix'), asi.mltimg
    
    asi = sread_thg_mlt(utr>time_double('2013-06-07/05:00'), sites, /half, type = 'asf')
    get_data, tvar, uts, imgs
    store_data, tvar, [uts,sfmepoch(asi.epoch,'unix')], [imgs,asi.mltimg]
    
    
    get_data, tvar, uts, imgs
    nrec = n_elements(uts)
    for i = 0, nrec-1 do begin
    
        ofn = shomedir()+'/thg_asi/thg_rb_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
        ;    ofn = 0
        sgopen, ofn, xsize = xsz, ysize = ysz
        
        device, decomposed = 0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        timg = reform(imgs[i,*,*])
        timg = congrid(timg, xsz, ysz)
        
        sgtv, timg, position = tpos
        sgset_map, position = tpos, color = white, xrange = xr
        
        xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
        
        for j = 0, n_elements(probes)-1 do begin
            pre0 = 'rbsp'+probes[j]+'_'
            get_data, pre0+'fpt_mlat', tuts, mlat
            get_data, pre0+'fpt_mlt', tuts, mlt
            mlat = interpol(mlat, tuts, uts[i], /quadratic)
            mlt = interpol(mlt, tuts, uts[i], /quadratic)
            plots, mlt*15, mlat, psym = syms[j], color = colors[j]
            xyouts, tpos[0], tpos[1]+ychsz*(j+1), /normal, 'RBSP-'+strupcase(probes[j]), color = colors[j], charsize = chsz
            
        endfor
        
        sgclose
    endfor
endif

end
