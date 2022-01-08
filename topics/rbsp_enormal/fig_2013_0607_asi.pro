
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


get_data, 'asf_info', tmp, info
mlats = info.mlats
mlts = info.mlts
imgsz = info.imgsz
get_data, 'asf_mos', uts, mos, pxidx




; fpt photon count.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlt', tuts, scmlts
    get_data, pre0+'fpt_mlat', tuts, scmlats
    scmlts = interpol(scmlts[*,0], tuts, uts)*15-360
    scmlats = interpol(scmlats[*,0], tuts, uts)
    
    phcnts = fltarr(nrec)
    for i = 0, nrec-1 do begin
        tmlat = scmlats[i]
        tmlt = scmlts[i]
        idx = where(abs(mlats-tmlat) le 0.1 and abs(mlts-tmlt) le 0.1, cnt)
        if cnt eq 0 then continue
        
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        phcnts[i] = max(timg[idx])
    endfor
    store_data, pre0+'asf_cnt', uts, phcnts, limits = $
        {ytitle:'Photon cnt!C@Footprint', yrange:[0,300], yticks:2, ystyle:1, yminor:5}
    get_data, pre0+'pf_fac_mat_map', tuts, dat
    store_data, pre0+'pfb_fac_mat_map', tuts, dat[*,0], limits = $
        {ytitle:'S@100 km!C0.1-812 sec!C(mV/m!U2!N)', yrange:[-15,15], ystyle:1, yticks:2, yminor:5, labels:'S!D||!N'}
endforeach


vars = ['asf_cnt','pfb_fac_mat_map']
vars = ['rbspa_'+vars, 'rbspb_'+vars]
nvar = n_elements(vars)

ofn = shomedir()+'/fig_asf_pf.pdf'
sgopen, ofn, xsize = 8.5, ysize = 5, /inch
device, decomposed = 0

poss = sgcalcpos(nvar)
loadct2, 43
tplot, vars, vlab_margin = 12, position = poss
sgclose

stop



; ewogram.
tmlat = 65.1
yr = [21,23]*15
mlt0s = smkarthm(yr[0],yr[1],0.1,'dx')-360
nmlt0 = n_elements(mlt0s)


ewos = fltarr(nrec,nmlt0)
for i = 0, nrec-1 do begin
    idx = where(abs(mlats-tmlat) le 0.1, cnt)
    if cnt lt 3 then continue
    timg = fltarr(imgsz)
    timg[pxidx] = mos[i,*]
    
    timg = timg[idx]
    tmlts = mlts[idx]
    idx = sort(tmlts)
    timg = timg[idx]
    tmlts = tmlts[idx]
    
    ewos[i,*] = interpol(timg,tmlts,mlt0s)
endfor

tvar = 'ewogram_'+sgnum2str(tmlat)
store_data, tvar, uts, ewos, mlt0s/15+24, $
    limits = {spec:1, no_interpo:1, $
    ytitle:'Ewogram!C@'+sgnum2str(tmlat)+' deg!CMLT (hr)', yrange:yr/15, $
    yticks:2, yminor:5, $
    ztitle:'Norm.count', zrange:[0,256], zticks:2}





ofn = 0
ofn = shomedir()+'/fig_ewogram.pdf'
sgopen, ofn, xsize = 8.5, ysize = 5, /inch
device, decomposed = 0
loadct2, 43

tutr = time_double(['2013-06-07/04:54','2013-06-07/04:58'])

poss = sgcalcpos(2, ypad = 5)

tvar = 'ewogram_'+sgnum2str(tmlat)
tplot, tvar, trange = utr, position = poss[*,0], /noerase, vlab_margin = 12
timebar, tutr, color = 6

tplot, tvar, trange = tutr, position = poss[*,1], /noerase, vlab_margin = 12

txs = time_double(['2013-06-07/04:54:45','2013-06-07/04:55:53'])
tys = [22.59,21.65]
plot, tutr, yr/15, xstyle = 5, ystyle = 5, /nodata, /noerase, position = poss[*,1]
plots, txs, tys, color = 255, thick = 4

sgclose

end
