;+
; MLT image keo, ewo, etc.
;-


time_range = time_double(['2015-02-18/02:05','2015-02-18/02:15'])
sites = ['gbay']
mlt_range = [-2.5,-0.5]
mlat_range = [61.5,64.5]
zrange = [60,180]
best_model = 't89'

probe = 'a'
rbsp_read_orbit, time_range, probe=probe, coord='gsm'
prefix = 'rbsp'+probe+'_'
r_var = prefix+'r_gsm'
read_geopack_info, r_var, model=best_model, direction=-1, refine=1


; EWOgram.
themis_read_mltimg_ewo, time_range, mlat_range=mlat_range, sites=sites
ewo_var = 'thg_asf_ewo'

yrange = mlt_range
ystep = 1
ytickv = make_bins(yrange, ystep, /inner)
yticks = n_elements(ytickv)-1
yminor = 5
options, ewo_var, 'yrange', yrange
options, ewo_var, 'yticks', yticks
options, ewo_var, 'ytickv', ytickv
options, ewo_var, 'yminor', yminor
options, ewo_var, 'zrange', zrange


; KEOgram.
themis_read_mltimg_keo, time_range, mlt_range=mlt_range, sites=sites
keo_var = 'thg_asf_keo'

yrange = mlat_range
ystep = 1
ytickv = make_bins(yrange, ystep, /inner)
yticks = n_elements(ytickv)-1
yminor = 4
options, keo_var, 'yrange', yrange
options, keo_var, 'yticks', yticks
options, keo_var, 'ytickv', ytickv
options, keo_var, 'yminor', yminor
options, keo_var, 'zrange', zrange


;    pflux_grant_read_pflux, time_range, probe=probe
;    pflux_grant_read_efield, time_range, probe=probe
;    pflux_grant_read_bfield, time_range, probe=probe
;    rbsp_read_en_spec, time_range, probe=probe
;    rbsp_read_pa_spec, time_range, probe=probe


plot_file = 0
sgopen, plot_file, xsize=4, ysize=6

vars = [prefix+['e_mgse','b1_gsm','e_en_spec'], $
    [ewo_var,keo_var]]
nvar = n_elements(vars)

margins = [8,4,8,1]
poss = sgcalcpos(nvar, margins=margins)
tplot, vars, trange=time_range, position=poss

pos_vars = prefix+['fmlt','fmlat']+'_'+best_model
foreach var_id, nvar+[-2,-1] do begin
    var = vars[var_id]
    tpos = poss[*,var_id]
    xrange = time_range
    yrange = get_setting(var, 'yrange')
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, $
        position=tpos, nodata=1, noerase=1
    
    get_data, pos_vars[var_id-nvar+2], xxs, yys
    oplot, xxs, yys, color=sgcolor('white')
endforeach


end