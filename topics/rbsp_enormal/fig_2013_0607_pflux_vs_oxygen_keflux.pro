
; plot vsc and bmag.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
;utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
probes = ['a','b']
spinrate = 11
perp = '!9'+string(94b)+'!X'
labfac = ['||',perp+',West',perp+',North']
dr0 = 1d/16



device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 0
tplot_options, 'ygridstyle', 0
tplot_options, 'xticklen', 0
tplot_options, 'xgridstyle', 0


; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1


foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    plot_hope_l3_keflux, utr, 'oxygen', probe = tprobe
    get_data, pre0+'pf_fac_mat', uts, dat, limits = lim
    store_data, pre0+'pfb_fac_mat', uts, dat[*,0], limits = {ytitle:lim.ytitle}
endforeach



tvar = 'rbsp?_keflux_oxygen'
options, tvar, 'labels', 'O!U+!N KE!C  in situ'
options, tvar, 'yminor', 5
tvar = 'rbspa_keflux_oxygen'
options, tvar, 'yrange', [-0.04,0.01]
tvar = 'rbspb_keflux_oxygen'
options, tvar, 'yrange', [-0.04,0.01]

tvar = 'rbsp?_pfb_fac_mat'
options, tvar, 'labels', 'S!D||!N!C  in situ'

tvar = 'rbspa_pfb_fac_mat'
options, tvar, 'yrange', [-0.05,0.2]

tvar = 'rbspb_pfb_fac_mat'
options, tvar, 'yrange', [-0.2,0.8]


ofn = 0
ofn = shomedir()+'/fig_pflux_vs_o_ke.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
;sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43


vars = ['o_en','o_pa','keflux_oxygen','pfb_fac_mat']
nvar = n_elements(vars)

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

poss = sgcalcpos(nvar, region = [0,0.5,1,1])
tplot, 'rbspa_'+vars, trange = utr, /novtitle, position = poss, /noerase, title = 'RBSP-A'

xyouts, poss[2,-2]-xchsz*2, poss[1,-2]+ychsz*0.5, /normal, alignment = 1, 'Upward'
xyouts, poss[2,-1]-xchsz*2, poss[3,-1]-ychsz*1.2, /normal, alignment = 1, 'Downward'


poss = sgcalcpos(nvar, region = [0,0,1,0.5])
tplot, 'rbspb_'+vars, trange = utr, /novtitle, position = poss, /noerase, title = 'RBSP-B'

xyouts, poss[2,-2]-xchsz*2, poss[1,-2]+ychsz*0.5, /normal, alignment = 1, 'Upward'
xyouts, poss[2,-1]-xchsz*2, poss[3,-1]-ychsz*1.2, /normal, alignment = 1, 'Downward'

sgclose

end
