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

yr = [60,70]
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
    endfor
    
    tvar = (pre0+'keogram_'+models[modelidx])[0]
    store_data, tvar, uts, keos, mlat0s, $
        limits = {spec:1, no_interp:1, $
        ytitle:'RBSP-'+strupcase(tprobe)+'!CFptMLat!C(deg)', $
        yrange:yr, ztitle:'Photon count', yticks:2, zrange:[50,300], zticks:1}
endforeach





utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])


vars = ['n_combine','e_en','o_en','keogram_t01','pf_fac_mat_map','e_en_mageis']
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

ofn = 0
ofn = shomedir()+'/fig_overview6.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43

options, 'rbspa_pf_fac_mat_map', 'yrange', [-20,90]
options, 'rbspa_pf_fac_mat_map', 'yticks', 2
options, 'rbspb_pf_fac_mat_map', 'yrange', [-20,180]
options, 'rbspb_pf_fac_mat_map', 'yticks', 2


poss = sgcalcpos(nvar*2+2)
pre0 = 'rbspa_'
tposs = poss[*,0:nvar-1]
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-A'
pre0 = 'rbspb_'
tposs = poss[*,nvar+1:nvar+1+nvar-1]
for i = 0, nvar-1 do tposs[[1,3],i] -= 0.05
tplot, pre0+vars, var_label = pre0+labs, trange = utr, $
    /noerase, position = tposs, vlab_margin = 12, title = 'RBSP-B'

sgclose


end
