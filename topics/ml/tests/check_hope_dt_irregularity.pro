;+
; Xiangning said HOPE time tags have irregularity. Let's check.
;-

time_range = time_double(['2016','2016-01-01/01:00'])
probe = 'b'
test = 0

prefix = 'rbsp'+probe+'_'

resolution = 5*60d
suffix = '_'+string(resolution/60,format='(I0)')+'min'
foreach species, ['e','p'] do begin
    flux_var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, $
        species=species, resolution=resolution, errmsg=errmsg)
    rename_var, flux_var, to=flux_var+suffix
endforeach

rbsp_read_en_spec, time_range, probe=probe
vars = prefix+['e','p']+'_en_spec'
options, vars, 'zlog', 1
options, vars, 'ylog', 1
foreach var, vars do begin
    get_data, var, times
    dtimes = times[1:-1]-times[0:-2]
    store_data, var+'_dt', times, [0,dtimes], limits={$
        ytitle:'(sec)', ylog:1, yrange:[10,100]}
endforeach


vars = prefix+'e_en_spec'+['','_5min']
zlim, vars, 1e4, 1e10, 1
ylim, vars, 4e0, 4e4, 1
ztickv = [1e4,1e7,1e10]
zticks = n_elements(ztickv)
options, vars, 'ztickv', ztickv
options, vars, 'zticks', zticks
options, vars, 'ztitle', 'e- (#/cm!U-2!N-s-sr-keV)'
options, vars, 'ytitle', 'Energy (eV)'

vars = prefix+'p_en_spec'+['','_5min']
zlim, vars, 1e0, 1e7, 1
ylim, vars, 4, 4e4, 1
ztickv = [1e0,1e3,1e6]
zticks = n_elements(ztickv)
options, vars, 'ztickv', ztickv
options, vars, 'zticks', zticks
options, vars, 'ztitle', 'H+ (#/cm!U-2!N-s-sr-keV)'
options, vars, 'ytitle', 'Energy (eV)'



vars = prefix+['e_en_spec'+['','_dt','_5min'],'p_en_spec'+['','_dt','_5min']]
nvar = n_elements(vars)
labels = letters(nvar)+') '+['e- '+['orig','dt','5min'],'H+ '+['orig','dt','5min']]
ypans = [1,0.5,1,1,0.5,1]


plot_file = join_path([googledir(),'works','ml_hope','plots','check_hope_dt_irregularity',$
prefix+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_to_')+'_check_hope_dt_irregularity_v01.pdf'])
if keyword_set(test) then plot_file = test
sgopen, plot_file, xsize=7, ysize=7, xchsz=xchsz, ychsz=ychsz
margins = [12,4,10,1]
poss = sgcalcpos(nvar, ypans=ypans, margins=margins)
tplot, vars, position=poss, trange=time_range
timebar, make_bins(time_range, 5*60), linestyle=1
for ii=0,nvar-1 do begin
    tpos = poss[*,ii]
    tx = xchsz*2
    ty = tpos[3]-ychsz*0.8
    xyouts, tx,ty,normal=1, labels[ii]
endfor

if keyword_set(test) then stop
sgclose

end