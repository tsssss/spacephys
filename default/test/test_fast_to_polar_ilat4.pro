; test on real data.
; use tplot.

; load fast spec.
fn = sdiskdir('Works')+'/data/cusp/fa_sdt_esa_1998_1001_08342.tplot'
tplot_restore, filename = fn

fn = sdiskdir('Works')+'/confs/seminar/seminar_2014_1021/data/seminar_dat.tplot'
tplot_restore, filename = fn
ofn = shomedir()+'/agu_2014_fig_eflux.eps'
eventid = '1998_1001_02'
logfile = sdiskdir('Works')+'/works/cusp/cusp_list_of_conjun.log'
info = cusp_read_conjun_list(logfile, event = eventid)
potr = info.polar.plot_time
fatr = info.fast.plot_time
potrcusp = info.polar.cusp_time
fatrcusp = info.fast.cusp_time
tplot_options, 'num_lab_min', 6
tplot_options, 'version', 2


var = 'pf_fac_mat_para_map'
get_data, 'po_'+var, pot0, pvar
get_data, 'fa_'+var, fat0, fvar

get_data, 'po_ilat', tmp, poilat & poilat = interpol(poilat,tmp,pot0)
get_data, 'fa_ilat', tmp, failat & failat = interpol(failat,tmp,fat0)


sgtruecolor
red = sgcolor('red')
green = sgcolor('green')
white = sgcolor('white')
black = sgcolor('black')

; fast var.
fy = fvar
fx = failat      ; x.
ft = fat0        ; t.

; ensure that fast ilat is monotonic.
idx = where(ft ge fatr[0] and ft le fatr[1])
ft = ft[idx]
fx = fx[idx]
fy = fy[idx]
fidx = idx

; polar var.
py = pvar        ; y.
px = poilat      ; x.
pt = pot0        ; t.

yr = [0,9]
xr = [55,80]

; **** map ft to pt using x as bridge.
ftp = interpol(pt,px,fx, /spline)
ftpp = interpol(ft,fx,px, /spline)
xr = sgcalcrange(minmax([ftp,pt]))
yr = sgcalcrange(minmax([fy,py]))
fatr = potr;[ftp[0],ftp[n_elements(ftp)-1]]

vars = ['fa_ele_keflux_map','fa_ion_keflux_map','fa_pf_fac_mat_para_map','fa_ilat']
fvars = ['feeflux','fieflux','fpflux','fa_ilat']
for i = 0, n_elements(vars)-1 do begin
    get_data, vars[i], t0, val, limits = lim
    t0 = interpol(ftp,ft,t0)
    store_data, fvars[i], t0, val, limits = lim
endfor
get_data, 'ion_en_spec', t0, val, tmp, limits = lim
t0 = interpol(ftp,ft,t0)
store_data, 'fien', t0, val, tmp, limits = lim
fvars = ['fien','feeflux','fieflux','fpflux']
flabs = ['fa_ilat']

x1 = pt
y1 = py
x2 = ftp
y2 = fy
tt = [[pt],[px],[ftpp]]

store_data, 'ftpp', pt, ftpp

; polar variables.
stplot_renew, 'po_ele_keflux_map', newname = 'peeflux'
stplot_renew, 'po_ion_keflux_map', newname = 'pieflux'
stplot_renew, 'po_pf_fac_mat_para_map', newname = 'ppflux'

pregion = [0d,0,1,0.5]
pvars = ['peeflux','pieflux','ppflux']
nvar = n_elements(pvars)
ppos = transpose(sgcalcpos(nvar, region = pregion, margins = [30,6,15,0]))
plabs = ['po_ilat']

; fast variables.

fregion = [0d,0.5,1,1]
fvars = ['fien','feeflux','fieflux','fpflux']
nvar = n_elements(fvars)
fpos = transpose(sgcalcpos(nvar, region = fregion, margins = [30,6,15,5]))
flabs = ['fa_ilat']

sgpsopen, shomedir()+'/test_fast2polar_ilat4.eps', xsize = 600, ysize = 900
sgtruecolor & !p.color = black & !p.background = white
sgindexcolor, 43
tplot, pvars, var_label = plabs, position = ppos, trange = potr, /noerase
tplot, fvars, var_label = flabs, position = fpos, trange = fatr, /noerase, uttick = 'ftpp'
sgpsclose,/pdf

end
