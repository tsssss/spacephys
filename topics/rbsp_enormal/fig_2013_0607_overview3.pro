; plot particle related parameters, Vsc, HOPE density, 
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





utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])


vars = ['vsc','n','e_en','h_en','e_en_mageis']
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

ofn = 0
ofn = shomedir()+'/fig_overview3.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43

tvars = ['rbspa','rbspb']+'_n'
foreach tvar, tvars do begin
    options, tvar, 'yrange', [0.02,2]
    options, tvar, 'ystyle', 1
endforeach


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


end
