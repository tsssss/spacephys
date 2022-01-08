
fac = ['b','e','n']
rgb = [6,4,2]
tplot_options, 'labflag', -1
tplot_options, 'xticklen', -0.02
tplot_options, 'yticklen', -0.01


tprobe = 'd'
pre0 = 'th'+tprobe+'_'
savelist = []
utr0 = time_double(['2014-08-28/09:30','2014-08-28/11:30'])

; read density, and pos, vbulk, E and B in GSM.
;read_themis_pos, utr0, probe=tprobe, addto=savelist
;read_themis_density, utr0, probe=tprobe, addto=savelist
;read_themis_vbulk, utr0, probe=tprobe, addto=savelist
;read_themis_efield, utr0, probe=tprobe, addto=savelist
;read_themis_bfield, utr0, probe=tprobe, addto=savelist

; preprocess B field: separate background B and dB, calc |B|.
; will return pre_[db_gsm,b0_gsm,bmod_gsm,bmag]
posvar = pre0+'pos_gsm'
stplot_prep_bfield, pre0+'b_gsm', posvar=posvar, addto=savelist
; preprocess E field: remove background.
stplot_prep_efield, pre0+'edot0_gsm', newname=pre0+'de_gsm', addto=savelisttlimit

; convert to fac.
vars = pre0+['de_gsm','db_gsm','vbulk_gsm']
newnames = pre0+['de_fac','db_fac','vbulk_fac']
labels = ['FAC dE','FAC dB','FAC V']
foreach tvar, vars, i do stplot_tofac, tvar, posvar=posvar, bvar=pre0+'bmod_gsm', addto=savelist, newname=newnames[i], label=labels[i]
; calc pflux.
stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scinfo
options, pre0+'pf_fac', 'labels', 'In-situ S'+fac
; calc mapcoef.
scalc_map_coef, pre0+'pos_gsm', pre0+'bmod_gsm', model='t89', coord='gsm', /igrf, prefix=pre0, dir=-1
; map pflux
get_data, pre0+'pf_fac', uts, dat
get_data, pre0+'map_coef', uts, mapc
for i=0, 2 do dat[*,i] *= mapc
store_data, pre0+'pf_para_map', uts, dat[*,0], limits={ytitle:'(mW/m!U2!N)', labels:'S!D||!N@100km', constant:0}
store_data, pre0+'pf_fac_map', uts, dat, limits={ytitle:'(mW/m!U2!N)', labels:'Mapped S'+fac, colors:rgb}

strmu = '!9'+string(109b)+'!X'
tvar = pre0+'pf_fac_mor_spec_1'
get_data, tvar, uts, dat, val, limits=lim
tvar = pre0+'pf_spec'
store_data, tvar, uts, dat*1e3, val, limits=lim
zrng = max(abs(dat))*0.5*1e3
options, tvar, 'zrange', [-1,1]*zrng
options, tvar, 'ztitle', '('+strmu+'W/m!U2!N)'

sgopen, srootdir()+'/plot_th'+tprobe+'_overview.pdf'
device, decomposed=0
loadct2, 43

vars = pre0+['density','bmag','vbulk_fac','de_fac','db_fac','pf_fac','pf_spec','pf_para_map']
tplot, vars, trange=utr0, /novtitle
sgclose

end