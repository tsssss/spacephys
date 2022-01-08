
; plot vsc and bmag.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
probes = ['a','b']
spinrate = 11
perp = '!9'+string(94b)+'!X'
labfac = ['||',perp+',West',perp+',North']
dr0 = 1d/16

filters = [0,6,18,132,812]    ; determined by mat spectrogram.
filters = [0,6,132,812]    ; determined by mat spectrogram.
filters = [0,4,20,132,812]    ; determined by mat spectrogram.
nfilter = n_elements(filters)-1
sclinfo = [0.6,1200,50]
filterstrs = ['1/'+string(1d/dr0,format='(I0)'),string(filters[1:*],format='(I0)')]
matidstrs = string(findgen(nfilter)+1,format='(I0)')


device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



; **** break pos_gsm into pos_gsm_[xyz].
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pos_gsm', uts, posgsm
    store_data, pre0+'pos_gsm_z', uts, posgsm[*,2], limits = $
        {ytitle:'Z GSM'}
endforeach




; **** control tplot var parameters.
pre0 = 'rbspa_'
; pre0 = 'rbspb_'

tvar = ['rbspa_','rbspb_']+'de_fac'
options, tvar, 'labels', 'dE!D'+labfac
tvar = ['rbspa_','rbspb_']+'db_fac'
options, tvar, 'labels', 'dB!D'+labfac

tvar = pre0+'de_fac'
options, tvar, 'yrange', [-60,60]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5

tvar = pre0+'db_fac'
options, tvar, 'yrange', [-20,20]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5

tvar = pre0+['de','db']+'_mat_spec'
options, tvar, 'zticks', 2
options, tvar, 'zcharsize', 0.9
options, tvar, 'color', 0
options, tvar, 'yrange', sclinfo[0:1]
options, tvar, 'ystyle', 1
options, tvar, 'yticklen', -0.015
options, tvar, 'ygridstyle', 0
options, tvar, 'yticks', 3
options, tvar, 'ytickv', [1,10,100,1000]
options, tvar, 'ytickname', ['0','1','2','3']
options, tvar, 'yminor', 5
options, tvar, 'ytitle', 'Log!D10!NP!C(sec)'

tvar = pre0+'db_mat_spec'
options, tvar, 'zrange', [0,0.6]



ofn = 0
ofn = shomedir()+'/fig_rbspa_pflux.pdf'
sgopen, ofn, xsize = 11, ysize = 5, /inch
;sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43

!x.gridstyle = 1
!x.ticklen = 1

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

tvar = ['_fac','_mat_spec']

figlabs = ['a.','b.']
vars = [pre0+'de'+tvar]
nvar = n_elements(vars)
poss = sgcalcpos(nvar, region = [0,0.45,0.5,1]+0.05*[0,1,0,1])
tplot, vars, trange = utr, position = poss, vlab_marg = 12, /noerase
plot, utr, sclinfo[0:1], /nodata, /noerase, position = poss[*,1], $
    xstyle = 1, xtickformat='(A1)', xticklen = 0, $
    ystyle = 1, ytickformat='(A1)', yticklen = 0, ylog = 1
for i = 0, n_elements(filters)-1 do plots, utr, filters[i]+[0,0], color = 255
for i = 0, nvar-1 do begin
    tx = poss[0,i]-xchsz*10
    ty = poss[3,i]-ychsz*0.5
    xyouts, tx, ty, /normal, alignment = 0, figlabs[i]
endfor

figlabs = ['c.','d.']
vars = [pre0+'db'+tvar]
nvar = n_elements(vars)
;poss = sgcalcpos(nvar, 2, region = [0,0.5,1,1])
poss = sgcalcpos(nvar, region = [0.5,0.45,1,1]+0.05*[0,1,0,1])
tplot, vars, trange = utr, position = poss, vlab_marg = 12, /noerase
plot, utr, sclinfo[0:1], /nodata, /noerase, position = poss[*,1], $
    xstyle = 1, xtickformat='(A1)', xticklen = 0, $
    ystyle = 1, ytickformat='(A1)', yticklen = 0, ylog = 1
for i = 0, n_elements(filters)-1 do plots, utr, filters[i]+[0,0], color = 255
for i = 0, nvar-1 do begin
    tx = poss[0,i]-xchsz*10
    ty = poss[3,i]-ychsz*0.5
    xyouts, tx, ty, /normal, alignment = 0, figlabs[i]
endfor

tplot_options, 'ygridstyle', 0
tplot_options, 'yticklen', 0



; **** update pflux.
vars = pre0+['de','db','pf']+'_fac_mat?'
; store_data, vars, /delete
; stplot_calc_pflux_mat, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', filter = filters, scaleinfo = sclinfo


; map using any of the models.
get_data, pre0+'map_coef', tuts, coef
get_data, pre0+'pf_fac_mat'+matidstrs[0], uts
coef = interpol(coef[*,0], tuts, uts)

for i = 0, n_elements(matidstrs)-1 do begin
    tvar = pre0+'pf_fac_mat'+matidstrs[i]
    get_data, tvar, uts, dat, limits = lim
    tvar = pre0+'pf_map_fac_mat'+matidstrs[i]
    for j = 0, 2 do dat[*,j]*= coef
    store_data, tvar, uts, dat, limits = {colors:[6,4,2]}
endfor

; integrate the mapped pflux.
vars = pre0+'pf_map_fac_mat'+matidstrs
intpfluxparas = dblarr(nfilter)
for i = 0, nfilter-1 do begin
    tvar = vars[i]
    get_data, tvar, uts, dat
    val = total(dat[*,0])*sdatarate(uts)    ; in mJ/m^2.
    store_data, tvar, uts, dat, val
    intpfluxparas[i] = val
endfor


tvar = pre0+'de_fac_mat'+matidstrs
options, tvar, 'yrange', [-40,40]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 4

tvar = pre0+'db_fac_mat'+matidstrs
options, tvar, 'yrange', [-9,9]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 3

tvar = pre0+'pf_map_fac_mat'+matidstrs
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [-10,5,20]
options, tvar, 'yminor', 5
options, tvar, 'yrange', [-12,28]


poss = sgcalcpos(1,nfilter,position = [0,0,1,0.6], xpad = 0, ypad = 0)
figlabs = ['e','f','g']
for i = 0, nfilter-1 do begin
    tpos = reform(poss[*,i])
    vars = pre0+['de','db','pf_map']+'_fac_mat'+matidstrs[i]
    options, vars, 'ytitle', ''
    options, vars, 'labels', ['','','']
    nvar = n_elements(vars)
    tposs = sgcalcpos(nvar, region = tpos+[1,0,1,0]*xchsz*3.5*(nfilter-1-i), lmargin = 8, rmargin = 5, bmargin = 3)
    tplot, vars, trange = utr, position = tposs, /noerase, /novtitle, /nouttick
    xyouts, tposs[0,0], tposs[3,0]+ychsz*(0.2+1.0*0), alignment = 0, /normal, $
        'Bandpass: '+filterstrs[i]+'-'+filterstrs[i+1]+' sec'
    xyouts, tposs[0,0], tposs[3,0]+ychsz*(0.2+1.0*1), alignment = 0, /normal, $
        'Integrated parallel S: '+sgnum2str(intpfluxparas[i],ndec=0)+' mJ/m!U2!N'
    
    for j = 0, nvar-1 do $
        xyouts, tposs[0,0]-xchsz*5, tposs[3,j]-ychsz*0.5, alignment = 0, /normal, $
            figlabs[j]+string(i+1,format='(I0)')+'.'
    
    if i eq 0 then begin
        labs = ['dE (mV/m)','dB (nT)', 'S (mW/m!U2!N)']
        for j = 0, nvar-1 do begin
            xyouts, xchsz*3, (tposs[1,j]+tposs[3,j])*0.5-ychsz*0.5, /normal, alignment = 0, labs[j]
        endfor
    endif
    

    tmp = stplot_ebratio(vars[0], vars[1], deidx = 2, dbidx = 1, trange = utr, method = 'minmax')
    ebratio = tmp[0]
    
    xyouts, tposs[0,0], tposs[3,0]+ychsz*(0.2+1.0*2), alignment = 0, /normal, $
        'R(E/B): '+sgnum2str(ebratio,nsgn=2)+' km/s'
endfor
  

sgclose

end
