
_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
tplot_options, 'yticklen', 0
tplot_options, 'ygridstyle', 0
tplot_options, 'xticklen', 0
tplot_options, 'xgridstyle', 0
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



foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'e_en_mageis', uts, dat, val
    idx = smkarthm(0,2,7,'x0')
    nch = n_elements(idx)
    store_data, pre0+'mageis', uts, dat[*,idx], val[*,idx], limits = $
        {ytitle:'MAGEIS e!C(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)', $
        colors:findgen(nch), labels:reform(string(round(val[0,idx]),format='(I0)')+'keV'), $
        spec:0, yrange:[1,1e7], ylog:1, ystyle:1, yticks:3, ytickv:[1e1,1e3,1e5,1e7]}
endforeach

options, 'rbspa_de_fac', 'yrange', [-100,100]

vars = ['de_fac','pf_fac_mat','mageis']

ofn = 0
ofn = shomedir()+'/fig_pflux_vs_mageis.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
;sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43

poss = sgcalcpos(nvar, region = [0,0.5,1,1], bmargin=-4)
tplot, 'rbspa_'+vars, trange = utr, /novtitle, position = poss, /noerase, title = 'RBSP-A'

poss = sgcalcpos(nvar, region = [0,0,1,0.5], bmargin=-4)
tplot, 'rbspb_'+vars, trange = utr, /novtitle, position = poss, /noerase, title = 'RBSP-B'

sgclose

end
