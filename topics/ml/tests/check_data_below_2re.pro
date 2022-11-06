;+
; Data below 2 Re.
;-

time_range = time_double(['2017','2017-01-02'])
probe = 'a'
test = 0


prefix = 'rbsp'+probe+'_'

resolution = 5*60d
suffix = '_'+string(resolution/60,format='(I0)')+'min'
foreach species, ['e','p','o'] do begin
    flux_var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, $
        species=species, resolution=resolution, errmsg=errmsg)
    rename_var, flux_var, to=flux_var+suffix
endforeach

dis_var = ml_rbsp_read_dis(time_range, probe=probe)
                

vars = prefix+'e_en_spec_5min'
zlim, vars, 1e4, 1e10, 1
ylim, vars, 4e0, 4e4, 1
ztickv = [1e4,1e7,1e10]
zticks = n_elements(ztickv)
options, vars, 'ztickv', ztickv
options, vars, 'zticks', zticks
options, vars, 'ztitle', 'e- (#/cm!U-2!N-s-sr-keV)'
options, vars, 'ytitle', 'Energy (eV)'

vars = prefix+'p_en_spec_5min'
zlim, vars, 1e0, 1e7, 1
ylim, vars, 4, 4e4, 1
ztickv = [1e0,1e3,1e6]
zticks = n_elements(ztickv)
options, vars, 'ztickv', ztickv
options, vars, 'zticks', zticks
options, vars, 'ztitle', 'H+ (#/cm!U-2!N-s-sr-keV)'
options, vars, 'ytitle', 'Energy (eV)'

vars = prefix+'o_en_spec_5min'
zlim, vars, 1e0, 1e7, 1
ylim, vars, 4, 4e4, 1
ztickv = [1e0,1e3,1e6]
zticks = n_elements(ztickv)
options, vars, 'ztickv', ztickv
options, vars, 'zticks', zticks
options, vars, 'ztitle', 'O+ (#/cm!U-2!N-s-sr-keV)'
options, vars, 'ytitle', 'Energy (eV)'


options, dis_var, 'ytickv', [1,3,5]
options, dis_var, 'yticks', 2
options, dis_var, 'yminor', 2
options, dis_var, 'yrange', [1,6]
options, dis_var, 'constant', 2

plot_file = join_path([googledir(),'works','ml_hope','plots','check_data_below_2re',$
    prefix+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_to_')+'_check_data_below_2re_v01.pdf'])
if keyword_set(test) then plot_file = test
sgopen, plot_file, xsize=7, ysize=5, xchsz=xchsz, ychsz=ychsz
vars = [prefix+['e','p','o']+'_en_spec_5min',dis_var]
nvar = n_elements(vars)
ypans = [1,1,1,0.5]
labels = letters(nvar)+') '+['e-','H+','O+','|R|']
margins = [12,4,10,1]
poss = sgcalcpos(nvar, ypans=ypans, margins=margins)
tplot, vars, trange=time_range, position=poss

get_data, dis_var, times, dis
index = where(dis le 2)
perigee_times = times[time_to_range(index,time_step=1)]
timebar, perigee_times, linestyle=0

for ii=0,nvar-1 do begin
    tpos = poss[*,ii]
    tx = xchsz*2
    ty = tpos[3]-ychsz*0.8
    xyouts, tx,ty,normal=1, labels[ii]
endfor

if keyword_set(test) then stop
sgclose

end