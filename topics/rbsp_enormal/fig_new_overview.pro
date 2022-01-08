; Vsc, HOPE density, S, keogram, 
; electron and ion energy spectrogram and MAGEIS electron energy spectrogram.


_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'yticklen', -0.01
tplot_options, 'xticklen', -0.02
tplot_options, 'zcharsize', 0.8
tplot_options, 'ycharsize', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1




models = ['t89','t04s','t01','t96']
modelidx = 2
probes = ['a','b']
nprobe = n_elements(probes)

diputa = time_double('2013-06-07/04:54:30')
diputb = time_double('2013-06-07/04:54:59')

; make gsm.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'b_gse', uts, bgse
    bgsm = sgse2gsm(bgse, stoepoch(uts,'unix'))
    bcomb = bgsm
    ndim = 3
    b0s = dblarr(ndim)
    b0s = [200,-100,0,180]
    b0s = round(bcomb[0,*]/5)*5
    labs = ['B!D'+['X','Y','Z']]
    for i = 0, ndim-1 do begin
        bcomb[*,i] = bcomb[*,i]-b0s[i]
        if b0s[i] eq 0 then continue
        if b0s[i] gt 0 then begin
            labs[i] = labs[i]+'!N-'+sgnum2str(b0s[i])+' nT'
        endif else begin
            labs[i] = labs[i]+'!N+'+sgnum2str(-b0s[i])+' nT'
        endelse
    endfor
    
    store_data, pre0+'b_combine', uts, bcomb, limits = $
        {ytitle:'(nT)', colors:[6,4,2,0], labels:labs, $
        yrange:[-20,60], ystyle:1, yticks:2, yminor:5}
endforeach


yr = [60,68]
mlat0s = smkarthm(yr[0],yr[1],0.1,'dx')
nmlat0 = n_elements(mlat0s)

get_data, 'asf_info', tmp, info
mlats = info.mlats
mlts = info.mlts
imgsz = info.imgsz
get_data, 'asf_mos', uts, mos, pxidx
nrec = n_elements(uts)

; make keogram.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlat', tuts, scmlats
    get_data, pre0+'fpt_mlt', tuts, scmlts
    
    scmlats = reform(scmlats[*,modelidx])
    scmlts = reform(scmlts[*,modelidx])
    
    scmlats = interpol(scmlats,tuts,uts)
    scmlts = interpol(scmlts,tuts,uts)*15-360
    
    keos = fltarr(nrec,nmlat0)
    cnts = fltarr(nrec)
    for i = 0, nrec-1 do begin
        idx = where(abs(mlts-scmlts[i]) le 0.1, cnt)
        if cnt lt 3 then continue
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        
        timg = timg[idx]
        tmlats = mlats[idx]
        idx = sort(tmlats)
        timg = timg[idx]
        tmlats = tmlats[idx]
        
        keos[i,*] = interpol(timg,tmlats,mlat0s)
        
        tmp = min(abs(mlat0s-scmlats[i]), idx)
        cnts[i] = keos[i,idx]
        keos[i,idx] = !values.f_nan
    endfor
    
    tvar = (pre0+'keogram_'+models[modelidx])[0]
    store_data, tvar, uts, keos, mlat0s, $
        limits = {spec:1, no_interp:1, $
        ytitle:'MLat!C(deg)', $
        yrange:yr, ztitle:'Photon count', yticks:2, zrange:[50,300], zticks:1, ystyle:1}
    tvar = (pre0+'count_'+models[modelidx])[0]
    store_data, tvar, uts, cnts, $
        limits = {ytitle:'Photon count',yrange:[100,500], yminor:5, labels:'Photon!C  @SC footpnt'}
endforeach


; calc mapped pflux.
pflabs = 'S!D'+['||','!9^!X,West','!9^!X,North']

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pf_fac', uts, dat, limits = lim
    get_data, pre0+'map_coef', tuts, tmp
    tmp = interpol(tmp[*,modelidx],tuts, uts)
    for i = 0, 2 do dat[*,i] *= tmp
    store_data, pre0+'pf_fac_mor_map', uts, dat, limits = lim
    options, pre0+'pf_fac_mor_map', 'labels', pflabs
    
    store_data, pre0+'pf_fac_mor_map', uts, dat[*,0], limits = $
        {ytitle:'(mW/m!U2)', labels:pflabs[0]+'!N @100km', colors:0}
endforeach




vars = tnames('rbsp?_e_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'zrange', [1e6,1e9]
options, vars, 'ytickname', '1e'+['2','3','4']
options, vars, 'ztickname', '1e'+['6','7','8','9']

vars = tnames('rbsp?_b_combine')
options, vars, 'constant', !values.d_nan

vars = tnames('rbsp?_n_combine')
options, vars, 'ytitle', '(cm!U-3!N)'

vars = tnames('rbsp?_pf_fac_mor_map')
options, vars, 'ytitle', '(mW/m!U2!N)'
options, vars, 'labels', 'Map to!C100 km'
options, vars, 'constant', !values.d_nan

vars = tnames('rbsp?_e_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'yrange', [200,2e4]
options, vars, 'ytickv', [1000,10000]
options, vars, 'ytickname', '1e'+['3','4']
options, vars, 'zrange', [1e6,1e9]
options, vars, 'zticks', 3
options, vars, 'ztickname', '1e'+['6','7','8','9']

vars = tnames('rbsp?_o_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'yrange', [40,4e4]
options, vars, 'ytickv', [100,1000,10000]
options, vars, 'ystyle', 1
options, vars, 'yminor', 5
options, vars, 'yticks', 2
options, vars, 'ytickname', '1e'+['2','3','4']
options, vars, 'zrange', [1e4,1e7]
options, vars, 'zticks', 3
options, vars, 'ztickname', '1e'+['4','5','6','7']

vars = 'rbsp?_keogram_t01'
options, vars, 'ytitle', 'MLat (deg)'
options, vars, 'yminor', 4


utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
;utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])


vars = ['e_en','n_combine','b_combine','keogram_t01','pf_fac_mor_map','o_en']
logidx = 4
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

ofn = 0
;ofn = shomedir()+'/fig_new_overview.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

lmarg = 15
cols = [6,4,2]
pos1 = [0.1,0.5,1,1]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1
pos2 = [0.1,0,1,0.5]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1

pre0 = 'rbspa_'
tposs = sgcalcpos(nvar, region=pos1)
tplot, pre0+vars[0:logidx-1], trange = utr, /nouttick, $
    /noerase, position = tposs[*,0:logidx-1], vlab_margin = lmarg, title = 'RBSP-A'
y1 = tposs[1,logidx] & y2 = (tposs[1,logidx]+tposs[3,logidx])*0.5
device, decomposed=1
polyfill, tposs[[0,2,2,0,0],logidx], [y1,y1,y2,y2,y1], /normal, color=sgcolor('silver')
device, decomposed=0
stplot_linlog, pre0+vars[logidx], trange = utr, $
    /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 3, $
    linyr = [-5,5], logyr = [5,200], logytickv = [10,100], linytickv = [-5,0,5]
xyouts, tposs[2,logidx]-xchsz*0.5, y2-ychsz*0.8, /normal, charsize=0.8, alignment=1, 'Linear'
xyouts, tposs[2,logidx]-xchsz*0.5, y2+ychsz*0.3, /normal, charsize=0.8, alignment=1, 'Log'
tplot, pre0+vars[logidx+1:*], trange = utr, /noerase, $
    position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs

lab2s = ['Electron','Density','B GSM', 'Keogram!C       @SC MLT', $
    'S parallel!C       1/4-1200s', 'Oxygen']

lab1s = ['a','b','c','d','e','f']+'-1.'
for i = 0, nvar-1 do begin
    xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
        lab1s[i]+' '+lab2s[i]
endfor

tidx = where(vars eq 'pf_fac_mor_map')
tpos = tposs[*,tidx]
xyouts, tpos[0]+xchsz*2, tpos[3]-ychsz*1.5, /normal, alignment = 0, 'Earthward', color = 0
plots, tpos[[0,2]], tpos[1]+(tpos[3]-tpos[1])*0.25, /normal, linestyle = 5

ttut = diputa
tx = (ttut-utr[0])/(utr[1]-utr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
plots, tx+[0,0], [ty1,ty2], /normal, color = 6
tidx = where(vars eq 'b_combine')
ty = tposs[3,tidx]-ychsz*1.5
xyouts, tx-xchsz*0.5, ty, /normal, alignment = 1, 'T_Dip', color = 0




pre0 = 'rbspb_'
tposs = sgcalcpos(nvar, region=pos2)
tplot, pre0+vars[0:logidx-1], trange = utr, /nouttick, $
    /noerase, position = tposs[*,0:logidx-1], vlab_margin = lmarg, title = 'RBSP-B'
y1 = tposs[1,logidx] & y2 = (tposs[1,logidx]+tposs[3,logidx])*0.5
device, decomposed=1
polyfill, tposs[[0,2,2,0,0],logidx], [y1,y1,y2,y2,y1], /normal, color=sgcolor('silver')
device, decomposed=0
stplot_linlog, pre0+vars[logidx], trange = utr, $
    /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 3, $
    linyr = [-5,5], logyr = [5,200], logytickv = [10,100], linytickv = [-5,0,5]
xyouts, tposs[2,logidx]-xchsz*0.5, y2-ychsz*0.8, /normal, charsize=0.8, alignment=1, 'Linear'
xyouts, tposs[2,logidx]-xchsz*0.5, y2+ychsz*0.3, /normal, charsize=0.8, alignment=1, 'Log'
tplot, pre0+vars[logidx+1:*], trange = utr, /noerase, $
    position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs

lab1s = ['a','b','c','d','e','f']+'-2.'
for i = 0, nvar-1 do begin
    xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
        lab1s[i]+' '+lab2s[i]
endfor

tidx = where(vars eq 'pf_fac_mor_map')
tpos = tposs[*,tidx]
xyouts, tpos[0]+xchsz*2, tpos[3]-ychsz*1.5, /normal, alignment = 0, 'Earthward', color = 0
plots, tpos[[0,2]], tpos[1]+(tpos[3]-tpos[1])*0.25, /normal, linestyle = 5

ttut = diputb
tx = (ttut-utr[0])/(utr[1]-utr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
;plots, tx+[0,0], [ty1,ty2], /normal, color = 6
tidx = where(vars eq 'b_combine')
ty = tposs[3,tidx]-ychsz*1.5
xyouts, tx-xchsz*0.5, ty+ychsz*0.5, /normal, alignment = 1, 'T_Dip?', color = 0
if size(ofn,/type) eq 7 then tmp = 100 else tmp = 5
arrow, tx, ty+ychsz, tx, ty-ychsz*0.2, /normalized, color=6, /solid, hsize=tmp

sgclose


end
