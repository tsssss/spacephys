
_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
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
nrec = n_elements(uts)

probes = ['a','b']

utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])

; fpt photon count.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlt', tuts, scmlts
    get_data, pre0+'fpt_mlat', tuts, scmlats
    scmlts = interpol(scmlts[*,2], tuts, uts)*15-360
    scmlats = interpol(scmlats[*,2], tuts, uts)
    
    phcnts = fltarr(nrec)
    for i = 0, nrec-1 do begin
        tmlat = scmlats[i]
        tmlt = scmlts[i]
        
        ; rbspa, 22.4802, 62.7122.
        ; rbspb, 22.1528, 61.3796.
        if tprobe eq 'a' then begin
            tmlat = 62.7122
            tmlt = 22.4802
        endif else begin
            tmlat = 61.3796
            tmlt = 22.1528
        endelse
        tmlt = (tmlt-24)*15
        
        idx = where(abs(mlats-tmlat) le 0.1 and abs(mlts-tmlt) le 0.1, cnt)
        if cnt eq 0 then continue
        
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        phcnts[i] = max(timg[idx])
    endfor
    store_data, pre0+'asf_cnt', uts, phcnts, limits = $
        {ytitle:'Photon cnt!C@Footprint', yrange:[0,650], yticks:2, ystyle:1, $
        yminor:5, constant:150, ytickv:[150,400,650]}
    
    tvar = pre0+'pfb_fac_mat_map'
    get_data, tvar, tuts, pflux
    get_data, 'asf_mos', uts
    idx = where(uts ge utr[0] and uts le utr[1], cnt)
    pflux = interpol(pflux,tuts,uts[idx])
    store_data, pre0+'pflux_tmp', uts[idx], pflux, limits = $
        {ytitle:'S@100 km!C3-812 sec!C(mV/m!U2!N)', yrange:[-10,40], $
        ystyle:1, yticks:2, yminor:5, labels:'S!D||!N', ytickv:[0,20,40]}
endforeach


vars = ['asf_cnt','pflux_tmp']
vars = ['rbspa_'+vars, 'rbspb_'+vars]
nvar = n_elements(vars)

ofn = shomedir()+'/fig_pflux_vs_asf.pdf'
; ofn = shomedir()+'/fig_pflux_vs_asf_bestcorr.pdf'
; ofn = 0
sgopen, ofn, xsize = 5, ysize = 6, /inch
device, decomposed = 0

poss = sgcalcpos(nvar)
loadct2, 43
tplot, vars, vlab_margin = 12, position = poss, trange = utr
sgclose

stop

tvar = 'rbspb_pfb_fac_mat_map'
get_data, tvar, tuts, pflux
get_data, 'asf_mos', uts, mos, pxidx
idx = where(uts ge utr[0] and uts le utr[1], cnt)
pflux = interpol(pflux,tuts,uts[idx])

irng = [200,300]
jrng = [90,150]
di = 10
dj = 10
is = smkarthm(irng[0],irng[1],di,'dx') & ni = n_elements(is)
js = smkarthm(jrng[0],jrng[1],dj,'dx') & nj = n_elements(js)
corrs = dblarr(ni-1,nj-1)


for i = 0, ni-2 do begin
    for j = 0, nj-2 do begin
        tdat = dblarr(cnt)
        for k = 0, cnt-1 do begin
            timg = fltarr(imgsz)
            timg[pxidx] = mos[idx[k],*]
            tdat[k] = max(timg[is[i]:is[i]+di,js[j]:js[j]+dj])
            if tdat[k] eq 0 then break
        endfor
        if min(tdat) eq 0 then continue
        corrs[i,j] = c_correlate(pflux,tdat,0)
    endfor
endfor

tut = time_double('2013-06-07/04:55:06')
timg = fltarr(imgsz)
timg[pxidx] = mos[where(uts eq tut),*]

window, 0, xsize = (is[-1]-is[0])*2, ysize = (js[-1]-js[0])*2
tvscl, timg[is[0]:is[-1],js[0]:js[-1]], 0
tvscl, congrid(corrs,is[-1]-is[0],js[-1]-js[0]), 1
tvscl, congrid(corrs,is[-1]-is[0],js[-1]-js[0]), 2


print, max(corrs, tpos)
tpos = array_indices(corrs, tpos)
irng = [is[tpos[0]],is[tpos[0]+1]]
jrng = [js[tpos[1]],js[tpos[1]+1]]
di = 2
dj = 2
is = smkarthm(irng[0],irng[1],di,'dx') & ni = n_elements(is)
js = smkarthm(jrng[0],jrng[1],dj,'dx') & nj = n_elements(js)
corrs = dblarr(ni-1,nj-1)


for i = 0, ni-2 do begin
    for j = 0, nj-2 do begin
        tdat = dblarr(cnt)
        for k = 0, cnt-1 do begin
            timg = fltarr(imgsz)
            timg[pxidx] = mos[idx[k],*]
            tdat[k] = max(timg[is[i]:is[i]+di,js[j]:js[j]+dj])
            if tdat[k] eq 0 then break
        endfor
        if min(tdat) eq 0 then continue
        corrs[i,j] = c_correlate(pflux,tdat,0)
    endfor
endfor

print, max(corrs, tpos)
tpos = array_indices(corrs, tpos)
ti = is[tpos[0]]
tj = js[tpos[1]]
print, mlts[ti,tj]/15+24, mlats[ti,tj]

; rbspa, 22.4802, 62.7122.
; rbspb, 22.1528, 61.3796.

stop

end
