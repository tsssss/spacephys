
; plot vsc correlation.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
probes = ['a','b']
nprobe = n_elements(probes)
models = ['t89','t96','t01','t04s','ts07d']
nmodel = n_elements(models)


device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 2
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1

ofn = shomedir()+'/fig_bmod_vs_b.pdf'
;ofn = 0
sgopen, ofn, xsize = 11, ysize = 6, /inch

device, decomposed = 0
loadct2, 43


; data from http://rbspgway.jhuapl.edu/SGMagEphem.
tuts = smkarthm(utr[0],utr[1],5*60,'dx')
bxs = [159.299780d, 160.840912d, 163.935090d, 167.240134d, 169.697084d, 171.900571d, 175.095153d]
bys = [-81.089260d, -78.983080d, -77.698481d, -76.301046d, -74.826114d, -73.264838d, -71.746279d]
bzs = [22.380453d, 22.812183d, 22.674540d, 22.816816d, 23.126006d, 24.059139d, 25.160107d]
store_data, 'rbspb_bmod_gsm_ts07d', tuts, [[bxs],[bys],[bzs]]
store_data, 'rbspa_bmod_gsm_ts07d', tuts, dblarr(n_elements(tuts),3)+!values.d_nan


nvar = nmodel+1
poss = sgcalcpos(nprobe,nvar, ypad = 1, xpad = 2)

for i = 0, nprobe-1 do begin
    tprobe = probes[i]
    pre0 = 'rbsp'+tprobe+'_'

    vars = pre0+['b_gse','bmod_gse_'+models]
    nvar = n_elements(vars)
    
    for j = 0, nvar-1 do begin
        tvar = vars[j]
        if strpos(tvar,'ts07d') ne -1 then continue
        get_data, tvar, uts, dat
        tpos = strpos(tvar,'gse')
        tvar = strmid(tvar,0,tpos)+'gsm'+strmid(tvar,tpos+3)
        print, tvar
        store_data, tvar, uts, sgse2gsm(dat,stoepoch(uts,'unix'))
    endfor

    vars = pre0+['b_gsm','bmod_gsm_*']
    options, vars, 'yrange', [-100,300]
    options, vars, 'ystyle', 1
    options, pre0+'b_gsm', 'ytitle', 'RBSP-'+strupcase(tprobe)+'!C(nT)'
    options, pre0+'bmod_gsm_*', 'ytitle', ''
    options, vars, 'colors', [6,4,2]
        
endfor


pre0 = 'rbsp'+['a','b']+'_'
vars = ['b_gsm','bmod_gsm_'+models]
tits = ['B GSM',strupcase(models)]
for i = 0, nvar-1 do begin
    tplot, pre0+vars[i], vlab_margin = 10, trange = utr, $
        position = reform(poss[*,i,*]), /noerase, title = tits[i], /novtitle

endfor

sgclose


end
