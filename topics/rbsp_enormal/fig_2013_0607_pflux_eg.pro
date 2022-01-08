
; plot vsc and bmag.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
probes = ['a','b']
spinrate = 11
perp = '!9'+string(94b)+'!X'
labfac = ['||',perp+',West',perp+',North']
dr0 = 1d/16

filters = [0,16,132,812]    ; determined by mat spectrogram.
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
tplot_options, 'xticklen', 0
tplot_options, 'yticklen', 0
tplot_options, 'ygridstyle', 0
tplot_options, 'xgridstyle', 0
!x.ticklen = 0
!x.gridstyle = 0



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
pre0 = 'rbspb_'

tvar = ['de_fac','db_fac']
tvar = ['rbspa_'+tvar,'rbspb_'+tvar]
options, tvar, 'labels', 'dE!D'+labfac

tvar = pre0+'de_fac'
options, tvar, 'yrange', [-150,150]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5

tvar = pre0+'db_fac'
options, tvar, 'yrange', [-20,40]
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

tvar = pre0+'de_mat_spec'
options, tvar, 'zrange', [0,1.6]

tvar = pre0+'db_mat_spec'
options, tvar, 'zrange', [0,1.2]



ofn = 0
ofn = shomedir()+'/fig_pflux_eg.pdf'
sgopen, ofn, xsize = 11, ysize = 5, /inch
;sgopen, ofn, xsize = 11, ysize = 8.5, /inch

rgb = sgcolor(['red','green','blue'])

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


figlabs = ['a.','b.']
vars = pre0+['de_fac','db_fac']
options, vars, 'colors', rgb
nvar = n_elements(vars)
poss = sgcalcpos(3, region = [0,0,0.5,1])
tplot, vars, trange = utr, position = poss, vlab_marg = 12, /noerase, /novtitle, title = 'Original dE dot0 and dB in FAC'
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
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5

tvar = pre0+'db_fac_mat'+matidstrs
options, tvar, 'yrange', [-12,12]
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 4

tvar = pre0+'pf_map_fac_mat'+matidstrs
options, tvar, 'ystyle', 1
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5


options, pre0+'de_fac_mat1', 'yrange', [-120,120]
options, pre0+'de_fac_mat2', 'yrange', [-40,40]
options, pre0+'de_fac_mat3', 'yrange', [-40,40]
options, pre0+'de_fac_mat4', 'yrange', [-15,15]


options, pre0+'pf_map_fac_mat1', 'yrange', [-50,100]
options, pre0+'pf_map_fac_mat2', 'yrange', [-50,100]
options, pre0+'pf_map_fac_mat3', 'yrange', [-10,70]
options, pre0+'pf_map_fac_mat4', 'yrange', [-8,16]


poss = sgcalcpos(3,region = [0.5,0,1,1])
figlabs = ['e.','f.','g.']
vars = pre0+['de','db','pf_map']+'_fac_mat3'
nvar = n_elements(vars)
options, vars, 'colors', rgb


options, vars[0], 'ytitle', 'dE FAC!C(mV/m)'
options, vars[1], 'ytitle', 'dB FAC!C(nT)'
options, vars[2], 'ytitle', 'S FAC!C(mW/m!U2!N)'



tplot, vars, trange = utr, position = poss, /noerase, /novtitle, title = 'Fields filtered in 20-132 sec'
for i = 0, nvar-1 do begin
    tx = poss[0,i]-xchsz*10
    ty = poss[3,i]-ychsz*0.5
    xyouts, tx, ty, /normal, alignment = 0, figlabs[i]
endfor


sgclose

end
