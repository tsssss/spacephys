; Vsc, HOPE density, S, keogram, 
; electron and ion energy spectrogram and MAGEIS electron energy spectrogram.


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


models = ['t89','t04s','t01','t96']
modelidx = 2
probes = ['a','b']
nprobe = n_elements(probes)

yr = [60,68]
mlat0s = smkarthm(yr[0],yr[1],0.1,'dx')
nmlat0 = n_elements(mlat0s)

get_data, 'asf_info', tmp, info
mlats = info.mlats
mlts = info.mlts
imgsz = info.imgsz
get_data, 'asf_mos', uts, mos, pxidx
nrec = n_elements(uts)


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


foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pf_fac_mat', uts, dat, limits = lim
    get_data, pre0+'map_coef', tuts, tmp
    tmp = interpol(tmp[*,modelidx],tuts, uts)
    for i = 0, 2 do dat[*,i] *= tmp
    store_data, pre0+'pf_fac_mat_map', uts, dat, limits = lim
    store_data, pre0+'pf_fac_mat_map2', uts, dat, limits = lim
endforeach





utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])


vars = ['keogram_t01','count_t01','pf_fac_mat_map2','pf_fac_mat_map']
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

ofn = 0
ofn = shomedir()+'/fig_pflux_vs_keogram_and_count.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


tvar = ['rbspa_','rbspb_']+'pf_fac_mat_map'
options, tvar, 'yrange', [-5,5]
options, tvar, 'yticks', 2
options, tvar, 'ytitle', ''
options, tvar, 'labels', ['','','']

tvar = ['rbspa_','rbspb_']+'pf_fac_mat_map2'
options, tvar, 'yrange', [5,100]
options, tvar, 'ylog', 1
options, tvar, 'ystyle', 1
options, tvar, 'ytitle', ''
options, tvar, 'labels', ['','','']


pflabs = 'S!D'+['||','!9^!X,West','!9^!X,North']
cols = [6,4,2]

poss = sgcalcpos(2)
pre0 = 'rbspa_'
tposs = sgcalcpos(nvar-1, region=[0,0.5,1,1])
tpos = tposs[*,nvar-2]
xyouts, tpos[0]-xchsz*4, 0.5*(tpos[1]+tpos[3]), /normal, alignment = 0.5, orientation = 90, '(mW/m!U2-!N)'
tys = tpos[3]-(tpos[3]-tpos[1])*0.25*[1,2,3]
for i = 0, 2 do $
    xyouts, tpos[2]+xchsz, tys[i], /normal, alignment = 0, pflabs[i], color = cols[i]
tposs = [[tposs[*,0:nvar-3]],[sgcalcpos(2,position=tposs[*,nvar-2],ypad=0)]]
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-A'
pre0 = 'rbspb_'
tposs = sgcalcpos(nvar-1, region=[0,0,1,0.5])
tpos = tposs[*,nvar-2]
xyouts, tpos[0]-xchsz*4, 0.5*(tpos[1]+tpos[3]), /normal, alignment = 0.5, orientation = 90, '(mW/m!U2-!N)'
tys = tpos[3]-(tpos[3]-tpos[1])*0.25*[1,2,3]
for i = 0, 2 do $
    xyouts, tpos[2]+xchsz, tys[i], /normal, alignment = 0, pflabs[i], color = cols[i]
tposs = [[tposs[*,0:nvar-3]],[sgcalcpos(2,position=tposs[*,nvar-2],ypad=0)]]
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-B'

sgclose


end
