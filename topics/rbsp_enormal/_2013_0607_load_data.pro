;+
; loads RBSP and ASI data.
; RBSP E/B at 16 Hz, in situ and mapped position, HOPE and MAGEIS spectrogram
; ASI keogram.
;-
pro _2013_0607_load_data

; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
models = ['t89','t96','t01','t04s']  ; the B0 model ~ measured. tried t04s, t96, t89, B0 are different.
nmodel = n_elements(models)
model0 = 't01'
probes = ['a','b']
allvars = []
dr0 = 1d/16
edot0 = 1
ionminen = 30
eleminen = 200

reload = 0



; one time settings.
vars = ['b_gse','e_gse','pos_gse','pos_gsm','mlt','lshell','mlat','dis','vsc','q_uvw2gse']
datvars = ['rbspa_'+vars, 'rbspb_'+vars]    ; basic variables.

vars = ['e_en_mageis','e_en','e_pa','h_en','h_pa','o_en','o_pa','ni','n','t','n_efw','n_combine']
ptlvars = ['rbspa_'+vars, 'rbspb_'+vars]    ; HOPE and MAGEIS data.

vars = []
vars = ['fpt_'+['mlat','mlon','mlt'],'map_coef','bmod_gse_'+models]
vars = [vars,'de_fac','db_fac','db_gse','b0_gse']
dervars = ['rbspa_'+vars, 'rbspb_'+vars]    ; derived variables.


asivars = 'asf_'+['info','mos']

mapvars = []
vars = ['grid_vertical','grid_equatorial','par']
foreach model, models do mapvars = [mapvars,model+'_'+vars]



; on pflux calc.
scaleinfo = {s0:0.25d, s1:1200d, dj:1d/8, ns:0d}
vars = [['de','db']+'_mor_spec','pf_fac','pf_para_mor_spec']
pfvars = ['rbspa_'+vars, 'rbspb_'+vars]


allvars = [datvars,ptlvars,dervars,pfvars,asivars,mapvars]




device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
dir = -1



datfn = shomedir()+'/psbl_de_32hz/2013_0607_data.tplot'
load = 0

foreach tvar, allvars do if tnames(tvar) eq '' then load = 1
if load then if file_test(datfn) ne 0 then begin
    store_data, '*', /delete
    tplot_restore, filename = datfn
endif


load = 0
if reload then load = 1
if load eq 1 then begin
    if file_test(datfn) ne 0 then begin
        store_data, '*', /delete
        tplot_restore, filename = datfn
    endif
    defsysv,'!rbsp_spice', exists=flag
    if flag eq 0 then rbsp_load_spice_kernels, trange = utr0
endif


; load high res E/B field.
load = 0
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        

        
        ; load position vars.
        tuts = smkarthm(utr0[0],utr0[1],60,'dx')
        rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = tuts, /no_spice_load
        get_data, pre0+'state_pos_gse', tmp, posgse
        posgse = posgse*re1
        tdat = posgse
        store_data, pre0+'state_*', /delete
        

        tvar = pre0+'mlt'
        tdat = atan(posgse[*,1],posgse[*,0])*deg
        tdat = (tdat+360) mod 360   ; convert to 0-360.
        tdat = (tdat/15 + 12) mod 24
        store_data, tvar, tuts, tdat, limits = {ytitle:'MLT (hr)'}
        
        
        tvar = pre0+'lshell'
        tets = stoepoch(tuts,'unix')
        possm = sgse2sm(posgse, tets)
        mlat = atan(possm[*,2],sqrt(possm[*,0]^2+possm[*,1]^2)) ; in rad.
        dis = sqrt(possm[*,0]^2+possm[*,1]^2+possm[*,2]^2)
        tdat = dis/(cos(mlat)^2)
        store_data, tvar, tuts, tdat, limits = {ytitle:'L-shell'}
        
        store_data, pre0+'mlat', tuts, mlat*deg, limits = {ytitle:'MLat (deg)'}
        store_data, pre0+'dis', tuts, dis, limits = {ytitle:'Dist (Re)'}
        store_data, pre0+'pos_gse', tuts, posgse, limits = $
            {ytitle:'R GSE (Re)', colors:[6,4,2], labels:'GSE '+['x','y','z']}
        store_data, pre0+'pos_gsm', tuts, sgse2gsm(posgse,tets), limits = $
            {ytitle:'R GSM (Re)', colors:[6,4,2], labels:'GSM '+['x','y','z']}
        
        
        
        ; uniform time for E/B.
        uts = smkarthm(utr0[0], utr0[1], dr0, 'dx')
        
        
        ; load B.
        emfisis = sread_rbsp_emfisis_l3(utr0, type = 'hires', probe = tprobe, coord = 'gse')
        tuts = sfmepoch(emfisis.epoch,'unix')
        store_data, pre0+'b_gse', uts, sinterpol(emfisis.mag,tuts,uts), $
            limit = {ytitle:'B GSE!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
        
        
        ; load vsvy and calc e.
        tvar = pre0+'vsc'
        dat = sread_rbsp_efw_l2(utr0, probes = tprobe, type = 'vsvy')
        tuts = sfmepoch(dat.epoch, 'unix')
        nrec = n_elements(uts)
        vsvy = sinterpol(dat.vsvy, tuts, uts)
        ; calibrate vsvy, tried rbsp_efw_get_cal_params, but really no effect
        ; since offset is 0, gain is not used.
        ; end up using V12, may consider more complex ways.
        vsc = mean(vsvy[*,0:1], dimension = 2)
        store_data, tvar, uts, vsc, limits = {ytitle:'Vsc!C(V)', labels:'Vsc'}
        store_data, pre0+'vsvy', uts, vsvy[*,0:3], limits = $
            {ytitle:'Vsvy!C(V)', colors:[1,2,3,4], labels:'V'+['1','2','3','4']}
        
        
        tvar = pre0+'e_uvw'
        ; the following is based on linear regression on v1,v2 and eu, v3,v4 and ev.
        eu = (vsvy[*,0]-vsvy[*,1])/100*1e3   ; V -> V/m -> mV/m.
        ev = (vsvy[*,2]-vsvy[*,3])/100*1e3
        ew = (vsvy[*,4]-vsvy[*,5])/12*1e3
        ew[*] = 0
        ; remove offset.
        spinrate = 11
        eu = eu-smooth(eu, spinrate/dr0, /edge_truncate, /nan)
        ev = ev-smooth(ev, spinrate/dr0, /edge_truncate, /nan)
        store_data, tvar, uts, [[eu],[ev],[ew]], limits = $
            {ytitle:'E UVW!C(mV/m)', labels:['Eu','Ev','Ew'], $
            colors:[6,4,2], yrange:[-200,200]}

        
        
        ; calc epoch times for spice (different than the epoch in cdfs).
        tmp = time_string(uts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+uts-uts[0]
        
        scid = strupcase(pre0+'science')
        cspice_pxform, scid, 'GSE', tets, muvw2gse
        quvw2gse = mtoq(transpose(muvw2gse))
        store_data, pre0+'q_uvw2gse', uts, quvw2gse

        
        ex = eu*muvw2gse[0,0,*] + ev*muvw2gse[1,0,*] + ew*muvw2gse[2,0,*]
        ey = eu*muvw2gse[0,1,*] + ev*muvw2gse[1,1,*] + ew*muvw2gse[2,1,*]
        ez = eu*muvw2gse[0,2,*] + ev*muvw2gse[1,2,*] + ew*muvw2gse[2,2,*]
        
        store_data, pre0+'e_gse', uts, [[ex],[ey],[ez]], limits = $
            {ytitle:'dE GSE!C(mV/m)', colors:[6,4,2], labels:'GSE '+['x','y','z']}

    endforeach
    tplot_save, allvars, filename = datfn
endif



load = 0
foreach tvar, ptlvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        mageis = sread_rbsp_mageis_l3(utr0, probes = tprobe)
        uts = sfmepoch(mageis.epoch,'unix')
        fedu = total(mageis.fedu,3,/nan)
        store_data, pre0+'e_en_mageis', uts, fedu, mageis.fedu_energy, limits = $
            {spec:1, no_interp:1, ylog:1, zlog:1, yrange:[4e1,4e3], ystyle:1, $
            ytitle:'MAGEIS!CElectron!CEenergy!C(keV)', $
            zticks:3, zrange:[1e-2,1e5], $
            ztitle:'(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'}
        
        vars = ['PITCH_ANGLE',$
            'Epoch_Ele','HOPE_ENERGY_Ele',$
            'Epoch_Ion','HOPE_ENERGY_Ion',$
            'FEDU','FPDU','FODU']
            
        hopel3 = sread_rbsp_hope_l3(utr0, probes = tprobe, vars = vars)
        pas = hopel3.pitch_angle

        ; ion spec.
        uts = sfmepoch(hopel3.epoch_ion,'unix')
        nrec = n_elements(uts)
        ens = hopel3.hope_energy_ion

        dat = hopel3.fodu
        type = 'o'
        zr = [1e4,1e7]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le ionminen, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor
        
        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7
        
        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}
            
        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7
        
        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[3,3e4], ystyle:1, ylog:1, $
            ytitle:'Energy (eV)'}
        
        
        dat = hopel3.fpdu
        type = 'h'
        zr = [1e4,1e7]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le ionminen, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor
        
        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7
        
        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}
            
        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7
        
        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[30,3e4], ystyle:1, ylog:1, $
            ytitle:'Energy (eV)'}
        
        
        ; electron.
        uts = sfmepoch(hopel3.epoch_ele,'unix')
        nrec = n_elements(uts)
        ens = hopel3.hope_energy_ele
        
        dat = hopel3.fedu
        type = 'e'
        zr = [1e7,1e10]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le ionminen, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor
        
        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7
        
        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}
            
        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7
        
        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[200,2e4], ylog:1, ystyle:1, $
            ytitle:'Energy (eV)'}
        
        vars = pre0+['h_en','h_pa','e_en','e_pa','o_en','o_pa']
        options, vars, 'spec', 1
        options, vars, 'no_interp', 1
        options, vars, 'zlog', 1
        options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'
        
        vars = pre0+['h_pa','e_pa','o_pa']
        options, vars, 'yrange', [0,180]
        options, vars, 'ystyle', 1
        
        tvar = pre0+'e_en'
        options, tvar, 'ytitle', 'e- energy!C(eV)'
        tvar = pre0+'h_en'
        options, tvar, 'ytitle', 'H+ energy!C(eV)'
        
        tvar = pre0+'e_pa'
        options, tvar, 'ytitle', 'e- pitch!C(deg)'
        tvar = pre0+'h_pa'
        options, tvar, 'ytitle', 'H+ pitch!C(deg)'
        tvar = pre0+'o_pa'
        options, tvar, 'ytitle', 'O+ pitch!C(deg)'
        
        
        hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom', $
            vars = ['Epoch_Ion','Dens_o_30','Dens_p_30'])
        if size(hopemom,/type) ne 8 then return
        uts = sfmepoch(hopemom.epoch_ion,'unix')
        store_data, pre0+'ni', uts, [[hopemom.dens_p_30],[hopemom.dens_o_30]], $
            limits = {ytitle:'n!Di!N!C(cm!U-3!N)', ylog:1, constant:1, labels:$
        ['n!DH!N>30eV','n!DO!N>30eV'], colors:[0,6]}
        
        
        hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
        if size(hopemom,/type) ne 8 then return
        uts = sfmepoch(hopemom.epoch_ele,'unix')
        store_data, pre0+'n', uts, hopemom.dens_e_200, $
            limits = {ytitle:'n!De!N!C(cm!U-3!N)', ylog:1, constant:1, labels:'n!De!N>200eV'}
        store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
            limits = {ytitle:'T!De!N!C(eV)', ylog:1, colors:[6,0], labels:['Tperp','Tpara'], constant:1000}
        tvar = pre0+'t'
        get_data, tvar, uts, dat
        idx = where(dat eq 1e20, cnt)
        if cnt ne 0 then dat[idx] = !values.d_nan
        store_data, tvar, uts, dat
        

        ; fit Vsc to density.
        get_data, pre0+'vsc', uts, vsc
        get_data, pre0+'n', tuts, hopen
        idx = where(finite(hopen))
        tuts = tuts[idx]
        hopen = hopen[idx]
        hopen = sinterpol(hopen, tuts, uts, /nan)
        
        
        ; hopen = a*exp(vsc*b)
        ; alog(hopen) = alog(a)+vsc*b
        
        tx0 = vsc
        ty0 = alog(hopen)
        
        a = regress(tx0, ty0, const = b)
        ty1 = a[0]*vsc+b[0]
        
        ybar = mean(ty0)
        sstot = total((ty0-ybar)^2)
        ssres = total((ty0-ty1)^2)
        
        r2 = 1-ssres/sstot
        print, a, b, r2
        
        
        efwn = exp(ty1)
        store_data, pre0+'n_efw', uts, efwn, limits = $
            {ytitle:'Density!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
            labels:'n!De,EFW Vsc'}
        store_data, pre0+'n_combine', uts, [[efwn],[hopen]], limits = $
            {ytitle:'RBSP-'+strupcase(tprobe)+'!CDensity!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
            labels:['n!De,EFW Vsc!N', 'n!De,HOPE > 200 eV'], colors:[0,6]}
        
    endforeach
    tplot_save, allvars, filename = datfn
endif




; check if need to re-derive variables.
load = 0
foreach tvar, dervars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        
        ; make a copy to be interpolated.  
        get_data, pre0+'e_gse', uts      
        get_data, pre0+'pos_gse', tuts, dat
        store_data, pre0+'r_gse', uts, sinterpol(dat, tuts, uts)
        
        
        ; decompose B into B0 and dB.
        get_data, pre0+'b_gse', uts, bgse
        b0gse = bgse
        for i = 0, 2 do b0gse[*,i] = scalcbg(bgse[*,i])
        dbgse = bgse-b0gse
        store_data, pre0+'b0_gse', uts, b0gse, limits = {ytitle:'B0 GSE!C(nT)', $
            labels:['x','y','z'], colors:[6,4,2]}
        store_data, pre0+'db_gse', uts, dbgse, limits = {ytitle:'dB GSE!C(nT)', $
            labels:['x','y','z'], colors:[6,4,2]}
        
        
        ; calc edot0.
        if edot0 eq 1 then begin
            ; convert B into uvw.
            get_data, pre0+'b0_gse', uts, bgse
            
            xs = bgse[*,0]
            ys = bgse[*,1]
            zs = bgse[*,2]
            
            ; the matrix is effectively transposed
            bu = xs*muvw2gse[0,0,*] + ys*muvw2gse[0,1,*] + zs*muvw2gse[0,2,*]
            bv = xs*muvw2gse[1,0,*] + ys*muvw2gse[1,1,*] + zs*muvw2gse[1,2,*]
            bw = xs*muvw2gse[2,0,*] + ys*muvw2gse[2,1,*] + zs*muvw2gse[2,2,*]
            
            get_data, pre0+'e_uvw', uts, euvw
            eu = euvw[*,0]
            ev = euvw[*,1]
            ew = -(bu*eu+bv*ev)/bw
            
            store_data, pre0+'edot0_uvw', uts, [[eu],[ev],[ew]], limits = $
                {ytitle:'dE dot0 UVW!C(mV/m)', colors:[6,4,2], labels:['u','v','w']}
            
            
            ex = eu*muvw2gse[0,0,*] + ev*muvw2gse[1,0,*] + ew*muvw2gse[2,0,*]
            ey = eu*muvw2gse[0,1,*] + ev*muvw2gse[1,1,*] + ew*muvw2gse[2,1,*]
            ez = eu*muvw2gse[0,2,*] + ev*muvw2gse[1,2,*] + ew*muvw2gse[2,2,*]
            
            store_data, pre0+'e_gse', uts, [[ex],[ey],[ez]], limits = $
                {ytitle:'dE dot0 GSE!C(mV/m)', colors:[6,4,2], labels:'GSE '+['x','y','z']}
        endif
        
        
        
        ; convert to FAC.
        get_data, pre0+'r_gse', uts, rgse
        rhat = sunitvec(rgse)
        bhat = sunitvec(b0gse)
        phat = sunitvec(scross(rhat,bhat))
        vhat = scross(bhat,phat)
        
        dbfac = [[sdot(dbgse,bhat)],[sdot(dbgse,phat)],[sdot(dbgse,vhat)]]
        store_data, pre0+'db_fac', uts, dbfac, limits = {ytitle:'dB FAC!C(nT)', $
            labels:['b','p','v'], colors:[6,4,2]}
        
        get_data, pre0+'e_gse', uts, egse
        defac = [[sdot(egse,bhat)],[sdot(egse,phat)],[sdot(egse,vhat)]]
        store_data, pre0+'de_fac', uts, defac, limits = {ytitle:'dE FAC!C(mV/m)', $
            labels:['b','p','v'], colors:[6,4,2]}
        
        
        
        ; map for each model.
        foreach model, models do begin
            suf0 = '_'+model
            ; map to 100 km, northern hemisphere.
            scalc_map_coef, pre0+'pos_gse', pre0+'b0_gse', model = model, $
                coord = 'gse', /igrf, prefix = pre0, suffix = suf0, dir = -1
            
            get_data, pre0+'bmod_gsm_'+model, tuts, dat
            store_data, pre0+'bmod_gse_'+model, tuts, sgsm2gse(dat,stoepoch(tuts,'unix')), $
                limit = {colors:[6,4,2], labels:['x','y','z'], ytitle:'B GSE '+model+'!C(nT)'}
        endforeach
        
        vars = pre0+['map_coef','fpt_'+['mlat','mlon','mlt']]
        foreach tvar, vars do begin
            get_data, tvar+'_'+models[0], limits = lim
            stplot_merge, tvar+'_'+models, newname = tvar, limits = $
                {ytitle:lim.ytitle, labels:models, colors:findgen(nmodel)}
        endforeach
        
    endforeach
    
    ; update data file.
    if load eq 1 then tplot_save, allvars, filename = datfn
endif




; **** calc Poynting flux.
load = 0
foreach tvar, pfvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        tvar = pre0+'de_mor_spec'
        get_data, pre0+'de_fac', uts, vec
        store_data, pre0+'demag', uts, snorm(vec)
        stplot_mor, pre0+'demag', scaleinfo=scaleinfo, newname = tvar
        options, tvar, 'zrange', [0,1]
        options, tvar, 'ztitle', '|dE|!U2!N (mV/m)'
        options, tvar, 'ytitle', 'Period!C(sec)'
        
        tvar = pre0+'db_mor_spec'
        get_data, pre0+'db_fac', uts, vec
        store_data, pre0+'dbmag', uts, snorm(vec)
        stplot_mor, pre0+'dbmag', scaleinfo = scaleinfo, newname = tvar
        options, tvar, 'zrange', [0,1]
        options, tvar, 'ztitle', '|dB|!U2!N (nT)'
        options, tvar, 'ytitle', 'Period!C(sec)'
        
        ; instataneous pflux.
        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', $
            scaleinfo=scaleinfo
        
        ; save the parallel pflux spectrogram.
        stplot_renew, pre0+'pf_fac_mor_spec_1', newname=pre0+'pf_para_mor_spec'
    endforeach
    
    ; update data file.
    tplot_save, allvars, filename = datfn
endif





if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        get_data, pre0+'map_coef', tuts, tdat
        tdat = reform(tdat[*,where(models eq model0)])
        get_data, pre0+'pf_fac', uts, pf
        tdat = sinterpol(tdat, tuts, uts)
        for i = 0, 2 do pf[*,i] *= tdat
        store_data, pre0+'pf_fac_map', uts, pf, limits = $
            {ytitle:'S@100 km!C(mW/m!U2!N)', colors:[6,4,2], labels:['b','p','v']}
    
        get_data, pre0+'b_gse', uts, bgse
        bgsm = sgse2gsm(bgse,stoepoch(uts,'unix'))
        store_data, pre0+'b_gsm', uts, bgsm, limits = $
            {ytitle:'B GSM!C(nT)', colors:[6,4,2], labels:'B GSM '+['x','y','z']}
            
        tvar = pre0+'e_pa'
        options, tvar, 'ytitle', 'HOPE!CElectron!CPitch Angle!C(deg)'
        
        tvar = pre0+'e_en'
        options, tvar, 'ytitle', 'HOPE!CElectron!CEnergy!C(eV)'
        options, tvar, 'yticks', 2
        options, tvar, 'ytickv', [1e2,1e3,1e4]
        options, tvar, 'yminor', 10
        options, tvar, 'ystyle', 1
        options, tvar, 'yrange', [200,2e4]
        options, tvar, 'zrange', [1e4,1e10]
        options, tvar, 'zticks', 3
        
        tvar = pre0+'h_en'
        options, tvar, 'ytitle', 'HOPE!CProton!CEnergy!C(eV)'
        options, tvar, 'yticks', 2
        options, tvar, 'ytickv', [1e2,1e3,1e4]
        options, tvar, 'yminor', 10
        options, tvar, 'ystyle', 1
        options, tvar, 'yrange', [40,4e4]
        options, tvar, 'zrange', [1e4,1e7]
        options, tvar, 'zticks', 3
        
        tvar = pre0+'o_en'
        options, tvar, 'ytitle', 'HOPE!COxygen!CEnergy!C(eV)'
        options, tvar, 'yticks', 2
        options, tvar, 'ytickv', [1e2,1e3,1e4]
        options, tvar, 'yminor', 10
        options, tvar, 'ystyle', 1
        options, tvar, 'yrange', [40,4e4]
        options, tvar, 'zrange', [1e4,1e7]
        options, tvar, 'zticks', 3
        
        tvar = pre0+'vsc'
        options, tvar, 'yticks', 2
        
        tvar = pre0+'b_gsm'
        options, tvar, 'ystyle', 1
        options, tvar, 'yrange', [-150,300]
        options, tvar, 'yticks', 3
        options, tvar, 'yminor', 5
        
        tvar = pre0+'pf_fac_map'
        perp = '!9'+string(94b)+'!X'
        options, tvar, 'ystyle', 1
        options, tvar, 'yrange', [-90,180]
        options, tvar, 'yticks', 3
        options, tvar, 'yminor', 5
        options, tvar, 'ytitle', 'S @100 km!C(mW/m!U2!N)'
        options, tvar, 'labels', 'S!D'+['||',perp+',West',perp+',North']
        
    endforeach
endif



; **** load ASI.
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445_3sites.cdf'

load = 0
foreach tvar, asivars do if tnames(tvar) eq '' then load = 1


if load then begin
    if file_test(asifn) eq 0 then begin
        tfn = sread_thg_mosaic(utr, sites, type = 'asf', minlat = 55, dark = 0)
        file_copy, tfn, asifn
        file_delete, tfn
    endif
    
    cdfs = scdfread(asifn)
    
    uts = (*cdfs[0].value)
    mos = (*cdfs[1].value)
    midn= (*cdfs[2].value)
    mlt = (*cdfs[3].value)
    mlat= (*cdfs[4].value)
    imgsz = (*cdfs[5].value)
    pxidx = (*cdfs[6].value)
    minlat = (*cdfs[7].value)[0]
    
    
    img = bytarr(imgsz)
    nrec = n_elements(uts)
    
    txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
    tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
    txs = txs/imgsz[0]*2
    tys = tys/imgsz[0]*2
    mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
    mlts = atan(tys, txs)*deg+90    ; in deg.
    
    store_data, 'asf_mos', uts, mos, pxidx
    store_data, 'asf_info', 0, {imgsz:imgsz, minlat:minlat, mlats:mlats, mlts:mlts}
    
    ; update data file.
    tplot_save, allvars, filename = datfn
endif



load = 0
foreach tvar, mapvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1


;load = 1
if load then begin
    
    ; prepare model parameters.
    tutr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    foreach tmodel, models do sgeopack_par, tutr, tmodel


    get_data, 'rbspa_pos_gsm', tuts, argsm
    get_data, 'rbspb_pos_gsm', tuts, brgsm
    
    idx = where(tuts ge tutr[0] and tuts le tutr[1])
    argsm = argsm[idx,*]
    brgsm = brgsm[idx,*]
    tuts = tuts[idx]
    
    rs = (snorm(argsm[*,0:1])+snorm(brgsm[*,0:1]))*0.5
    
    ; use the start, middle and end times.
    dr0 = 3     ; data rate for aurora images.
    uts = smkarthm(tutr[0], tutr[1], dr0, 'dx')
    nrec = n_elements(uts)
    utidx = fix([0,(nrec-1)/2,nrec-1])
    uts = uts[utidx]
    nrec = n_elements(uts)
    
    
    ; ** map the vertical grid to ionosphere.
    nsct = 60+1   ; # of grid in azimuthal direction.
    nscz = 10+1   ; # of grid in z direction.
    
    scr0 = mean(rs)     ; the L of vertical plane.
    sctrng = [2.8,-1.7] ; 4.5 MLT wide.
    sczrng = [1.5,4.5]
    scts = smkarthm(sctrng[0],sctrng[1], nsct, 'n')*15*rad
    sczs = smkarthm(sczrng[0],sczrng[1], nscz, 'n')
    scxs = -scr0*cos(scts)
    scys =  scr0*sin(scts)
    
    scfmlats = dblarr(nrec,nsct,nscz)
    scfmlons = dblarr(nrec,nsct,nscz)
    scfmlts = dblarr(nrec,nsct,nscz)

    
    foreach tmodel, models do begin
        get_data, tmodel+'_par', tmp, dat
        pars = sinterpol(dat, tmp, uts)
        ets = stoepoch(uts,'unix')
        
        t89 = (tmodel eq 't89')? 1: 0
        t96 = (tmodel eq 't96')? 1: 0
        t01 = (tmodel eq 't01')? 1: 0
        t04s = (tmodel eq 't04s')? 1: 0
        storm = (tmodel eq 't04s')? 1: 0

        ; map in situ grid to ionosphere.
        for i = 0, nrec-1 do begin
            tet = ets[i]
            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
            
            par = reform(pars[i,*])
            for j = 0, nsct-1 do begin
                for k = 0, nscz-1 do begin
                    xp = scxs[j] & yp = scys[j] & zp = sczs[k]
                    
                    geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
                        /refine, /ionosphere, $
                        t89 = t89, t96 = t96, t01 = t01, ts04 = ts04, storm = storm
                        
                    ; convert from gsm to mag.
                    geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
                    
                    scfmlats[i,j,k] = asin(tzf/r0)*deg
                    scfmlons[i,j,k] = atan(tyf,txf)*deg
                    
                    scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
                endfor
            endfor
        endfor
        scfmlts = (scfmlts+24) mod 24
        ; grid mlt and mlat in ionosphere; grid r0, azimuthal angle, z in tail.
        info = {fmlt:scfmlts, fmlat:scfmlats, r0:scr0, ts:scts, zs:sczs}
        
        tvar = tmodel+'_grid_vertical'
        store_data, tvar, uts, info
    endforeach
    
    
    
    ; ** map the equatorial grid to ionosphere.
    neqt = 60+1   ; # of grid in MLT.
    neqr = 10+1   ; # of grid in R.

    tilt0 = 13      ; deg, based both on model (12.99) and eyeballing (12).    
    eqlrng = [4.5,6.5]      ; L-shell range.
    eqtrng = 24-[2.8,-1.7]  ; MLT range.
    eqrs = smkarthm(eqlrng[0],eqlrng[1],neqr, 'n')
    eqts = smkarthm(eqtrng[0],eqtrng[1],neqt, 'n')*15*rad
    
    scfmlats = dblarr(nrec,neqt,neqr)
    scfmlons = dblarr(nrec,neqt,neqr)
    scfmlts = dblarr(nrec,neqt,neqr)
    
    
    foreach tmodel, models do begin
        get_data, tmodel+'_par', tmp, dat
        pars = sinterpol(dat, tmp, uts)
        ets = stoepoch(uts,'unix')
        
        t89 = (tmodel eq 't89')? 1: 0
        t96 = (tmodel eq 't96')? 1: 0
        t01 = (tmodel eq 't01')? 1: 0
        t04s = (tmodel eq 't04s')? 1: 0
        storm = (tmodel eq 't04s')? 1: 0
        
        ; map in situ grid to ionosphere.
        for i = 0, nrec-1 do begin
            tet = ets[i]
            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
            
            par = reform(pars[i,*])
            for k = 0, neqr-1 do begin    
                for j = 0, neqt-1 do begin
                    xp =-eqrs[k]*cos(tilt0*rad)*cos(eqts[j])
                    yp =-eqrs[k]*cos(tilt0*rad)*sin(eqts[j])
                    zp = eqrs[k]*sin(tilt0*rad)
                    
                    geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
                        /refine, /ionosphere, $
                        t89 = t89, t96 = t96, t01 = t01, ts04 = ts04, storm = storm
                        
                    ; convert from gsm to mag.
                    geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
                    
                    scfmlats[i,j,k] = asin(tzf/r0)*deg
                    scfmlons[i,j,k] = atan(tyf,txf)*deg
                    
                    scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
                endfor
            endfor
        endfor
        scfmlts = (scfmlts+24) mod 24
        ; grid mlt and mlat in ionosphere; grid x, y, z in tail.
        info = {fmlt:scfmlts, fmlat:scfmlats, z0:tilt0, rs:eqrs, ts:eqts}
        
        tvar = tmodel+'_grid_equatorial'
        store_data, tvar, uts, info
    endforeach
    
    ; update data file.
    tplot_save, allvars, filename = datfn
endif

end
