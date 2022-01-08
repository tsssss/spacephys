


; settings.
utr = time_double(['2015-04-10/00:00','2015-04-10/01:00'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
asiutr = time_double(['2015-04-10/00:05','2015-04-10/00:35'])
load = 0
reload = 0



; one settings.
vars = ['bgse','egse','rgse','pos_gse','mlt','lshell','dis','vsc','q_uvw2gse']
datvars = ['rbspa_'+vars, 'rbspb_'+vars]    ; basic variables.

vars = ['defac','dbfac','dbgse','b0gse','fpt_mlt','fpt_mlat','fpt_mlon','bmodgse','map_coef']
dervars = ['rbspa_'+vars, 'rbspb_'+vars]    ; derived variables.

model = 't89'
model = 't04s'


dr0 = 1d/16
site = 'nrsq'                   ; the site under conjunction.
sites = ['nrsq']                ; sites for asf.
probes = ['a','b']
syms = [1,1]                    ; symbols for rbsp-a and -b.
colors = [6,255]                ; colors for rbsp-a and -b.

device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



datfn = shomedir()+'/psbl_de_32hz/'+time_string(utr[0],tformat='YYYY_MMDD')+'_data.tplot'
if tnames(datvars[0]) eq '' then if file_test(datfn) ne 0 then begin
    store_data, '*', /delete
    tplot_restore, filename = datfn
endif
foreach tvar, datvars do if tnames(tvar) eq '' then reload = 1
if load eq 1 then reload = 1

if reload eq 1 then begin
    if file_test(datfn) ne 0 then begin
        store_data, '*', /delete
        tplot_restore, filename = datfn
    endif
    rbsp_load_spice_kernels, trange = utr
endif


; load high res E/B field.
foreach tprobe, probes do begin
    if reload eq 0 then continue
    
    pre0 = 'rbsp'+tprobe+'_'
    
    emfisis = sread_rbsp_emfisis_l3(utr0, type = 'hires', probe = tprobe, coord = 'gse')
    store_data, pre0+'bgse', sfmepoch(emfisis.epoch,'unix'), emfisis.mag, $
        limit = {ytitle:'B GSE!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
        
    rbsp_preprocess_efield, utr, probes = tprobe, /no_spice_load
endforeach


; update data file.
if reload eq 1 then tplot_save, datvars, filename = datfn





; check if need to re-derive variables.
foreach tvar, dervars do if tnames(tvar) eq '' then reload = 1


; decompose b into b0 and db, convert to FAC, map coeff.

; uniform time for E/B.
uts = smkarthm(utr0[0], utr0[1], dr0, 'dx')
foreach tprobe, probes do begin
    if reload eq 0 then continue
    
    pre0 = 'rbsp'+tprobe+'_'
    
    get_data, pre0+'pos_gse', tuts, dat
    store_data, pre0+'rgse', tuts, dat      ; make a copy to be interpolated.
    
    ; get uniform times.
    vars = pre0+['egse','bgse','rgse']
    foreach tvar, vars do begin
        get_data, tvar, tuts, dat
        store_data, tvar, uts, sinterpol(dat, tuts, uts)
    endforeach
    
    
    ; decompose B into B0 and dB.
    get_data, pre0+'bgse', uts, bgse
    b0gse = bgse
    for i = 0, 2 do b0gse[*,i] = scalcbg(bgse[*,i])
    dbgse = bgse-b0gse
    store_data, pre0+'b0gse', uts, b0gse, limits = {ytitle:'B0 GSE!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
    store_data, pre0+'dbgse', uts, dbgse, limits = {ytitle:'dB GSE!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
    
    
    ; convert to FAC.
    get_data, pre0+'rgse', uts, posgse
    rhat = sunitvec(posgse)
    bhat = sunitvec(b0gse)
    phat = sunitvec(scross(rhat,bhat))
    vhat = scross(bhat,phat)
    
    dbfac = [[sdot(dbgse,bhat)],[sdot(dbgse,phat)],[sdot(dbgse,vhat)]]
    store_data, pre0+'dbfac', uts, dbfac, limits = {ytitle:'dB FAC!C(nT)', labels:['b','p','v'], colors:[6,4,2]}
    
    get_data, pre0+'egse', uts, egse
    defac = [[sdot(egse,bhat)],[sdot(egse,phat)],[sdot(egse,vhat)]]
    store_data, pre0+'defac', uts, defac, limits = {ytitle:'dE FAC!C(mV/m)', labels:['b','p','v'], colors:[6,4,2]}
    
    
    ; map to 100 km, northern hemisphere.
    scalc_map_coef, pre0+'pos_gse', pre0+'b0gse', model = model, coord = 'gse', /igrf, prefix = pre0, dir = -1
    
    get_data, pre0+'fpt_mlon', tuts, mlon
    mlt = slon2lt(mlon, stoepoch(tuts,'unix'), /mag, /deg)/15    ; in hour.
    mlt = (mlt+24) mod 24
    store_data, pre0+'fpt_mlt', tuts, mlt, $
        limits = {ytitle:'MLT/fpt (hr)'}
        
    get_data, pre0+'bmod_gsm', tuts, dat
    store_data, pre0+'bmodgse', tuts, sgsm2gse(dat,stoepoch(tuts,'unix')), $
        limit = {colors:[6,4,2], labels:['x','y','z'], ytitle:'B model GSE!C(nT)'}
endforeach

; update data file.
if reload eq 1 then tplot_save, [datvars,dervars], filename = datfn





stop

; **** load mltimg.
xsz = 800
ysz = xsz*0.5
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

reload = 1

cprobes = ['c1','c2','c3','c4']
ccolors = [6,5,4,3]
csyms = [1,1,1,1]


if reload then begin

    defsysv,'!rbsp_spice', exists=flag
    if flag eq 0 then rbsp_load_spice_kernels, trange = utr0

    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        tvar = pre0+'pos_gse'
        if tnames(tvar) ne '' then continue

        uts = smkarthm(utr0[0],utr0[1],30,'dx')
        rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = uts, /no_spice_load
        get_data, pre0+'state_pos_gse', tmp, posgse
        posgse = posgse*re1
        store_data, tvar, tmp, posgse

        if tnames(pre0+'_fpt_mlon') ne '' then continue
        scalc_map_coef, tvar, model = model, coord = 'gse', /igrf, prefix = pre0, dir = -1
        get_data, pre0+'fpt_mlon', tuts, mlon
        mlt = slon2lt(mlon, stoepoch(tuts,'unix'), /mag, /deg)/15    ; in hour.
        mlt = (mlt+24) mod 24
        store_data, pre0+'fpt_mlt', tuts, mlt, $
            limits = {ytitle:'MLT/fpt (hr)'}
    endforeach
    

    foreach tprobe, cprobes do begin
        pre0 = tprobe+'_'
        tvar = pre0+'pos_gse'
        if tnames(tvar) ne '' then continue

        ifn = shomedir()+'/Downloads/'+pre0+'pos_gse.cdf'
        cdf = scdfread(ifn)
        uts = sfmepoch(*cdf[0].value,'unix')
        posgse = (*cdf[1].value)*re1
        store_data, tvar, uts, posgse
        if tnames(pre0+'_fpt_mlon') ne '' then continue
        scalc_map_coef, tvar, model = model, coord = 'gse', /igrf, prefix = pre0, dir = -1
        get_data, pre0+'fpt_mlon', tuts, mlon
        mlt = slon2lt(mlon, stoepoch(tuts,'unix'), /mag, /deg)/15    ; in hour.
        mlt = (mlt+24) mod 24
        store_data, pre0+'fpt_mlt', tuts, mlt, $
            limits = {ytitle:'MLT/fpt (hr)'}
    endforeach


    tvar = 'mltimg'
    
    if tnames(tvar) eq '' then begin
        asi = sread_thg_mlt(asiutr, sites, /half, type = 'asf')
        store_data, tvar, sfmepoch(asi.epoch, 'unix'), asi.mltimg
    endif
    get_data, tvar, uts, imgs
    nrec = n_elements(uts)
    for i = 0, nrec-1 do begin
    
        ofn = shomedir()+'/thg_asi/'+time_string(utr[0],tformat='YYYY_MMDD')+ $
            '/thg_rb_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
; ofn = 0
        sgopen, ofn, xsize = xsz, ysize = ysz
        
        device, decomposed = 0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        timg = reform(imgs[i,*,*])
        timg = congrid(timg, xsz, ysz)
        
        sgtv, timg, position = tpos
        sgset_map, position = tpos, color = white, xrange = xr
        
        xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
        
        for j = 0, n_elements(probes)-1 do begin
            pre0 = 'rbsp'+probes[j]+'_'
            get_data, pre0+'fpt_mlat', tuts, mlat
            get_data, pre0+'fpt_mlt', tuts, mlt
            mlat = interpol(mlat, tuts, uts[i], /quadratic)
            mlt = interpol(mlt, tuts, uts[i], /quadratic)
            plots, mlt*15, mlat, psym = syms[j], color = colors[j]
            xyouts, tpos[0], tpos[1]+ychsz*(j+1), /normal, 'RBSP-'+strupcase(probes[j]), color = colors[j], charsize = chsz
        endfor
        
        
        for j = 0, n_elements(cprobes)-1 do begin
            pre0 = cprobes[j]+'_'
            get_data, pre0+'fpt_mlat', tuts, mlat
            get_data, pre0+'fpt_mlt', tuts, mlt
            mlat = interpol(mlat, tuts, uts[i], /quadratic)
            mlt = interpol(mlt, tuts, uts[i], /quadratic)
            plots, mlt*15, mlat, psym = csyms[j], color = ccolors[j]
            xyouts, tpos[0], tpos[1]+ychsz*(j+3), /normal, strupcase(cprobes[j]), color = ccolors[j], charsize = chsz
        endfor
        
        sgclose
    endfor
endif







end
