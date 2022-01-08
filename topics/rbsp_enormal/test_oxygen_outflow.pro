
utr1 = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
utr2 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
tut = time_double('2013-06-07/04:55')

ionminen = 30
eleminen = 200

; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1




models = ['t89','t04s','t01','t96']
modelidx = 2
probes = ['a','b']
nprobe = n_elements(probes)




_2013_0607_load_data
device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'yticklen', 0.01
tplot_options, 'xticklen', 0.02
tplot_options, 'zcharsize', 0.8




vars = ['o_pa','o_en','h_pa','h_en','e_pa','e_en']
specvars = ['rbspa_'+vars,'rbspb_'+vars]


device, decomposed = 0
loadct2, 43

tplot, specvars, trange = utr2

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    
    get_data, pre0+'o_pa', uts, dat, val
    nrec = n_elements(uts)
    
    alphas = dblarr(nrec)
    for i = 0, nrec-1 do begin
        tmp = max(reform(dat[i,*]), idx)
        alphas[i] = val[idx]
    endfor
    store_data, pre0+'o_alpha', uts, alphas
    
    
    get_data, pre0+'o_en', uts, dat, val
    nrec = n_elements(uts)
    
    ens = dblarr(nrec)
    for i = 0, nrec-1 do begin
        tidx = where(val[i,*] ge 100)
        tmp = max(reform(dat[i,tidx]), idx)
        ens[i] = val[i,tidx[idx]]
    endfor
    v0s = sqrt(2*ens*1.6e-19/(16*1.67e-27))*1e-3    ; in km/s.
    store_data, pre0+'o_energy', uts, ens
    
    
    
    ; trace use the middle time and pos.
    tet = stoepoch(tut,'unix')
    ; get sc position.
    get_data, pre0+'pos_gsm', tuts, dat
    tpos = sinterpol(dat, tuts, tut)
    get_data, 't01_par', tuts, pars
    par = sinterpol(pars, tuts, tut)
    
    r0 = 2d
    geopack_epoch, tet, year, mo, dy, hr, mi, sc, /breakdown_epoch
    geopack_recalc, year, mo, dy, hr, mi, sc, /date, tilt = tilt
    geopack_trace, tpos[0], tpos[1], tpos[2], -1, par, xf, yf, zf, $
        fline = fline, /igrf, r0 = r0, /t01
    geopack_igrf_gsm, fline[*,0], fline[*,1], fline[*,2], bx0, by0, bz0
    geopack_t01, par, fline[*,0], fline[*,1], fline[*,2], dbx, dby, dbz
    
    b0 = 240d 
    bgsms = [[bx0+dbx],[by0+dby],[bz0+dbz]]
    bmags = snorm(bgsms)
    bcoef = 1+(b0-min(bmags))/bmags
    for i = 0, 2 do bgsms[*,i] *= bcoef
    bmags = snorm(bgsms)
    
    bs = b0/sin(alphas*rad)^2
    
    r0s = dblarr(nrec,3)
    t0s = dblarr(nrec)
    for i = 0, nrec-1 do begin
        for j = 0, 2 do r0s[i,j] = sinterpol(fline[*,j], bmags, bs[i])
        idx = where(fline[*,0] lt r0s[i,0])
        tfline = [fline[idx,*],r0s[i,*]]
        tbmags = bmags[idx]
        tnrec = n_elements(tfline[*,0])
        t0s[i] = 0
        for j = 1, tnrec-2 do begin
            tds = snorm(tfline[j,*]-tfline[j-1,*])*re
            tdt = tds/(v0s[i]*sqrt(1-tbmags[j-1]^2/bs[i]^2))
            t0s[i] = t0s[i]+tdt
        endfor
    endfor
    
    t1s = smooth(t0s,3)
    r1s = snorm(r0s)
    get_data, pre0+'o_en', uts, dat
    for i = 0, nrec-1 do begin
        if max(dat[i,*]) lt 0 then begin
            t1s[i] = !values.d_nan
            r1s[i] = !values.d_nan
        endif
    endfor
    
    store_data, pre0+'dt', uts, t1s, limits = {psym:-1, ytitle:'Time lag!C(sec)'}
    store_data, pre0+'ds', uts, r1s, limits = {psym:-1, ytitle:'Energized R!C(Re)',yrange:[1,6]}

endforeach

;if tnames('rbspa_keflux_oxygen') eq '' then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        plot_hope_l3_keflux, utr2, 'oxygen', probe = tprobe, min_energy = ionminen
    endforeach
;endif


; calc mapped pflux.
pflabs = 'S!D'+['||','!9^!X,West','!9^!X,North']

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pf_fac_mat', uts, dat, limits = lim
    get_data, pre0+'map_coef', tuts, tmp
    tmp = interpol(tmp[*,modelidx],tuts, uts)
    for i = 0, 2 do dat[*,i] *= tmp
    store_data, pre0+'pf_fac_mat_map', uts, dat, limits = lim
    options, pre0+'pf_fac_mat_map', 'labels', pflabs
    
    store_data, pre0+'pf_fac_mat_map', uts, dat[*,0], limits = $
        {ytitle:'(mW/m!U2!N)', labels:pflabs[0]+'!N!C@100km', colors:0}
    
    get_data, pre0+'keflux_oxygen', uts, dat, limits = lim
    get_data, pre0+'map_coef', tuts, tmp
    tmp = interpol(tmp[*,modelidx],tuts, uts)
    dat = dat*tmp
    store_data, pre0+'keflux_oxygen_map', uts, dat, limits = $
        {ytitle:'(mW/m!U2!N)', labels:'KEflux O!U+!N!C  @100km', colors:0}
    
endforeach



pre0 = 'rbspb_'
diput = time_double('2013-06-07/04:54:40')


options, pre0+'o_en', 'yrange', [40,4e4]

tvar = pre0+'keflux_oxygen_map'
options, tvar, 'yrange', [-9.5,0.5]
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [-9,-4.5,0]
options, tvar, 'yminor', 5
options, tvar, 'ystyle', 1

tvar = pre0+'o_pa'
options, tvar, 'yticklen', -0.01
options, tvar, 'ytitle', 'Pitch!C(deg)'


tvar = pre0+'o_en'
options, tvar, 'yticklen', -0.01
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [100,1000,10000]
options, tvar, 'yminor', 5
options, tvar, 'ytickname', '10!U'+['2','3','4']
options, tvar, 'ytitle', 'Energy!C(eV)'

tvar = pre0+'dt'
get_data, tvar, uts, dat
store_data, tvar, uts, dat/60d
options, tvar, 'ytitle', '(min)'
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [0,4,8]
options, tvar, 'yminor', 4


tvar = pre0+'ds'
options, tvar, 'ytitle', '(Re)'
options, tvar, 'R!Dmax', '(Re)'
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [1,4,7]
options, tvar, 'yminor', 3



ofn = shomedir()+'/'+pre0+'o_outflow_vs_pflux.pdf'
ofn = 0
sgopen, ofn, xsize = 8, ysize = 6, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

lmarg = 12
pos1 = [0d,0,1,1]

vars = ['pf_fac_mat_map','keflux_oxygen_map','o_pa','o_en']
logidx = 0
nvar = n_elements(vars)
labs = ['mlt','lshell','mlat']

tposs = sgcalcpos(nvar, region=pos1)


tutr = utr1
stplot_linlog, pre0+vars[logidx], trange = tutr, $
    /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 2, $
    linyr = [-5,5], logyr = [5,100], logytickv = [10,200], linytickv = [-5,0,5]
tplot, pre0+vars[logidx+1:*], trange = tutr, /noerase, $
    position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs

ttut = diput
tx = (ttut-tutr[0])/(tutr[1]-tutr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
plots, tx+[0,0], [ty1,ty2], /normal, color = 6


get_data, pre0+'keflux_oxygen_map', tuts, dat
okeut = time_double('2013-06-07/04:57:17')
idx = where(tuts ge tutr[0] and tuts le tutr[1])
tuts = tuts[idx]
dat = dat[idx]
tmp = min(dat, idx)
okeut = tuts[idx]
ttut = okeut
tx = (ttut-tutr[0])/(tutr[1]-tutr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
plots, tx+[0,0], [ty1,ty2], /normal, color = 6, linestyle = 2

tidx = 2
ttx = tx+xchsz*1
tty = tposs[1,tidx]+ychsz*0.5
xyouts, ttx, tty+ychsz*0, /normal, alignment = 0, color = 255, $
    'dT = 3 min, Source below 2.7 Re'
xyouts, ttx, tty+ychsz*1.1, /normal, alignment = 0, color = 255, $
    'T_Max|KE| - T_Dip = '+sgnum2str((okeut-diput)/60,ndec=1)+' min'

figlabs = ['a. S!D||', 'b. KEflux O!U+', 'c. O!U+!N PA', 'd. O!U+!N EN']
for i = 0, nvar-1 do begin
    xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
        figlabs[i]
endfor

sgclose


stop



pre0 = 'rbspa_'
diput = time_double('2013-06-07/04:54:30')


options, pre0+'o_en', 'yrange', [40,4e4]

tvar = pre0+'keflux_oxygen_map'
options, tvar, 'yrange', [-3.5,0.5]
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [-3,-1.5,0]
options, tvar, 'yminor', 5
options, tvar, 'ystyle', 1

tvar = pre0+'o_pa'
options, tvar, 'yticklen', -0.01
options, tvar, 'ytitle', 'Pitch!C(deg)'


tvar = pre0+'o_en'
options, tvar, 'yticklen', -0.01
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [100,1000,10000]
options, tvar, 'yminor', 5
options, tvar, 'ytickname', '10!U'+['2','3','4']
options, tvar, 'ytitle', 'Energy!C(eV)'

tvar = pre0+'dt'
get_data, tvar, uts, dat
store_data, tvar, uts, dat/60d
options, tvar, 'ytitle', '(min)'
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [0,4,8]
options, tvar, 'yminor', 4


tvar = pre0+'ds'
options, tvar, 'ytitle', '(Re)'
options, tvar, 'R!Dmax', '(Re)'
options, tvar, 'yticks', 2
options, tvar, 'ytickv', [1,4,7]
options, tvar, 'yminor', 3



ofn = shomedir()+'/'+pre0+'o_outflow_vs_pflux.pdf'
ofn = 0
sgopen, ofn, xsize = 8, ysize = 6, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

lmarg = 12
pos1 = [0d,0,1,1]

vars = ['pf_fac_mat_map','keflux_oxygen_map','o_pa','o_en']
logidx = 0
nvar = n_elements(vars)
labs = ['mlt','lshell','mlat']

tposs = sgcalcpos(nvar, region=pos1)


tutr = utr1
stplot_linlog, pre0+vars[logidx], trange = tutr, $
    /noerase, position = tposs[*,logidx], lastpanel = 0, ytitlepos = 2, $
    linyr = [-5,5], logyr = [5,100], logytickv = [10,100], linytickv = [-5,0,5]
tplot, pre0+vars[logidx+1:*], trange = tutr, /noerase, $
    position = tposs[*,logidx+1:*], vlab_margin = lmarg, var_label = pre0+labs
    
ttut = diput
tx = (ttut-tutr[0])/(tutr[1]-tutr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
plots, tx+[0,0], [ty1,ty2], /normal, color = 6


get_data, pre0+'keflux_oxygen_map', tuts, dat
okeut = time_double('2013-06-07/04:57:17')
idx = where(tuts ge tutr[0] and tuts le tutr[1])
tuts = tuts[idx]
dat = dat[idx]
tmp = min(dat, idx)
okeut = tuts[idx]
ttut = okeut
tx = (ttut-tutr[0])/(tutr[1]-tutr[0])
tx = tposs[0,0]+(tposs[2,0]-tposs[0,0])*tx
ty1 = tposs[3,0]
ty2 = tposs[1,nvar-1]
plots, tx+[0,0], [ty1,ty2], /normal, color = 6, linestyle = 2

tidx = 2
ttx = tx+xchsz*1
tty = tposs[1,tidx]+ychsz*0.5
xyouts, ttx, tty+ychsz*0, /normal, alignment = 0, color = 255, $
    'dT = 3 min, Source below 2.7 Re'
xyouts, ttx, tty+ychsz*1.1, /normal, alignment = 0, color = 255, $
    'T_Max|KE| - T_Dip = '+sgnum2str((okeut-diput)/60,ndec=1)+' min'
    
figlabs = ['a. S!D||', 'b. KEflux O!U+', 'c. O!U+!N PA', 'd. O!U+!N EN']
for i = 0, nvar-1 do begin
    xyouts, tposs[0,i]-xchsz*lmarg, tposs[3,i]-ychsz*1, /normal, alignment = 0, $
        figlabs[i]
endfor

sgclose



end