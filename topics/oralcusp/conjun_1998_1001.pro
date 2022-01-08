rt = sdiskdir('Works')+'/data/cusp/'

fn_fa_fld = rt+'fa_sdt_fld_1998_1001_08343.tplot'
fn_fa_esa = rt+'fa_sdt_esa_1998_1001_08343.tplot'
tr_fa = time_double(['1998-10-01/04:18','1998-10-01/04:32'])

fn_po_fld = rt+'po_sdt_fld_1998_1001_03.sdt'
tr_po = time_double(['1998-10-01/02:00','1998-10-01/04:15'])

colors = [0,6]
load = 1
ps = 0

sgwindow
device, decomposed = 0
loadct2, 43
!p.font = 1
tplot_options, 'ygap', 0.25
tplot_options, 'ynozero', 1
tplot_options, 'version', 2
tplot_options, 'num_lab_min', 8
tplot_options, 'labflag', 1
tplot_options, 'zcharsize', 0.9
time_stamp, /off

; plot 1: fast esa.
if keyword_set(load) then tplot_restore, filename = fn_fa_esa
coef = 1.6e-9   ; convert # cm^-2 s^-1 to uA m^-2.
get_data, 'ion_nflux', t0, ion
get_data, 'ele_nflux', tmp, ele
if n_elements(ele) ne n_elements(ion) then ele = sinterpol(ele,tmp,t0)
store_data, 'fa_nflux', t0, [[ion],[ele]]*coef, $
    limits = {colors:colors, labels:['Ji','Je'], ytitle:'(uA/m!U-2!N)'}

get_data, 'ion_eflux', t0, ion
get_data, 'ele_eflux', tmp, ele
if n_elements(ele) ne n_elements(ion) then ele = sinterpol(ele,tmp,t0)
store_data, 'fa_eflux', t0, [[ion],[ele]], $
    limits = {colors:colors, labels:['KEi','KEe'], ytitle:'(mW/m!U-2!N)'}

get_data, 'ion_n', t0, ion
get_data, 'ele_n', tmp, ele
if n_elements(ele) ne n_elements(ion) then ele = sinterpol(ele,tmp,t0)
store_data, 'fa_density', t0, [[ion],[ele]], $
    limits = {colors:colors, labels:['ni','ne'], ytitle:'(cm!U-3!N)', ylog:1}
vars = ['ion_para_spec','ion_perp_spec','ion_anti_spec', $
    'fa_nflux','fa_eflux','fa_density','ele_en_spec']
labs = ['ilat','mlt','dis']
titl = 'FAST ESA, '+time_string(t0[0],tformat='YYYY-MM-DD/hh')
tplot, vars, trange = tr_fa, var_label = labs, title = titl
fn = shomedir()+'/'+time_string(t0[0],tformat='YYYY_MMDD_hh')+'_esa.eps'
if keyword_set(ps) then begin
    sgpsopen, fn
    tplot, vars, trange = tr_fa, var_label = labs, title = titl
    sgpsclose
endif

; plot 2: fast field.
if keyword_set(load) then begin
    fast_sdt_prep_poynting_flux, fn_fa_fld
    stplot_calc_poynting_flux, 'fa_de_fac', 'fa_db_fac', 'fa_pf_fac', $
        method = 'mat', filter = [3,10,30,100], scaleinfo = [0.5,250,60]
endif

; plot 3: polar field.
if keyword_set(load) then begin
    polar_sdt_prep_poynting_flux, fn_po_fld
    stplot_calc_poynting_flux, 'po_de_fac', 'po_db_fac', 'po_pf_fac', $
        method = 'mat', filter = [10,40,180,1500], scaleinfo = [0.5,2000,60]
endif

end