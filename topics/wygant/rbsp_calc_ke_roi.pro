
; test time.
reload = 0

utr = ['2013-05-01/07:00','2013-05-01/08:00']
prb = 'b'

;utr = ['2012-11-14/04:00','2012-11-14/05:00']
;prb = 'b'

utr = time_double(utr)
pre0 = 'rbsp'+prb+'_'



; settings.
device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1


; vars to load.
allvars = []

; load from data file.
var1s = pre0+['h_spec','o_spec','e_spec','o_180','h_180','hopel3']
allvars = [allvars,var1s]

; get from drawing lines.
var2s = pre0+['o_en_lim','h_en_lim']
allvars = [allvars,var2s]

rootdir = srootdir()
ofn = rootdir+'/hope_l3_ion.dat.tplot'
if file_test(ofn) eq 1 then tplot_restore, filename = ofn

load = 0
foreach tvar, var1s do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    vars = []
    ; read hope l3.
    hopel3 = sread_rbsp_hope_l3(utr, probe = prb)
    store_data, pre0+'hopel3', 0, hopel3
    
    pas = hopel3.pitch_angle
    
    uts = sfmepoch(hopel3.epoch_ele,'unix')
    val = hopel3.hope_energy_ele
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e10], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!Ce (eV)'}
    dat = total(hopel3.fedu, 3, /nan)
    tvar = pre0+'e_spec'
    store_data, tvar, uts, dat, val, limits = lim
    
    uts = sfmepoch(hopel3.epoch_ion,'unix')
    val = hopel3.hope_energy_ion
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e7], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!CH (eV)'}
    dat = total(hopel3.fpdu, 3, /nan)
    tvar = pre0+'h_spec'
    store_data, tvar, uts, dat, val, limits = lim

    idx = where(pas ge 140)
    dat = total(hopel3.fpdu[*,*,idx], 3, /nan)
    tvar = pre0+'h_180'
    store_data, tvar, uts, dat, val, limits = lim
    options, tvar, 'ytitle', lim.ytitle+'!Canti-para'
    
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e7], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!CO (eV)'}
    dat = total(hopel3.fodu, 3, /nan)
    tvar = pre0+'o_spec'
    store_data, tvar, uts, dat, val, limits = lim
    
    idx = where(pas ge 140)
    dat = total(hopel3.fodu[*,*,idx], 3, /nan)
    tvar = pre0+'o_180'
    store_data, tvar, uts, dat, val, limits = lim
    options, tvar, 'ytitle', lim.ytitle+'!Canti-para'

    tplot_save, allvars, filename = ofn
endif



load = 0
foreach tvar, var1s do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    tvar = pre0+['o_180','o_spec','h_spec','e_spec']
    tplot, tvar, trange = utr
    stplot_add_line, vxs, vys
    vxs = [utr[0],vxs,utr[1]]
    vys = [vys[0],vys,vys[-1]]
    get_data, pre0+'o_spec', uts
    enlims = sinterpol(vys,vxs,uts)
    store_data, pre0+'o_en_lim', uts, enlims, limits = $
        {ytitle:'O Energy Limit!C(eV)', ylog:1, yrange:[20,5e4]}
    store_data, pre0+'h_en_lim', uts, enlims, limits = $
        {ytitle:'H Energy Limit!C(eV)', ylog:1, yrange:[20,5e4]}
    
    tplot_save, allvars, filename = ofn
endif


get_data, pre0+'hopel3', tmp, hopel3

get_data, pre0+'o_en_lim', uts, enlims
plot_hope_l3_keflux, utr, probe = prb, 'oxygen', hopel3 = hopel3, max_energy = enlims

get_data, pre0+'h_en_lim', uts, enlims
plot_hope_l3_keflux, utr, probe = prb, 'proton', hopel3 = hopel3, max_energy = enlims




ofn = 0
ofn = shomedir()+'/rbspb_keflux_o_h_2013_0501_insitu.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


vars = pre0+['o_spec','o_180','keflux_oxygen']
nvar = n_elements(vars)
pos1 = [0,0.5,1,1]
poss = sgcalcpos(nvar, region = pos1, lmargin=15)
tplot, vars, position = poss, /noerase, title = 'Oxygen'

; add energy limit.
idx = [0,1]
get_data, pre0+'o_en_lim', uts, enlims, limit = lims
for i = 0, n_elements(idx)-1 do begin
    plot, utr, lims.yrange, /noerase, /nodata, position = poss[*,idx[i]], $
        ylog = 1, ystyle = 5, xstyle = 5
    oplot, uts, enlims, color = 6, thick = 2
endfor
xyouts, poss[0,2]+xchsz*3, poss[1,2]+ychsz*4, /normal, alignment = 0, 'upward', charsize = 1.2



vars = pre0+['h_spec','h_180','keflux_proton']
nvar = n_elements(vars)
pos2 = [0,0,1,0.5]
poss = sgcalcpos(nvar, region = pos2, lmargin=15)
tplot, vars, position = poss, /noerase, title = 'Proton'

; add energy limit.
idx = [0,1]
get_data, pre0+'h_en_lim', uts, enlims, limit = lims
for i = 0, n_elements(idx)-1 do begin
    plot, utr, lims.yrange, /noerase, /nodata, position = poss[*,idx[i]], $
        ylog = 1, ystyle = 5, xstyle = 5
    oplot, uts, enlims, color = 6, thick = 2
endfor
xyouts, poss[0,2]+xchsz*3, poss[3,2]-ychsz*5, /normal, alignment = 0, 'earthward', charsize = 1.2

sgclose

end
