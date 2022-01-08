


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
omniutr = time_double(['2013-06-06/12:00','2013-06-08/00:00'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
probes = ['a','b']
load = 0
longutr = 1




device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    emfisis = sread_rbsp_emfisis_l3(utr0, type = 'hires', probe = tprobe, coord = 'gse')
    tvar = pre0+'bmag'
    store_data, tvar, sfmepoch(emfisis.epoch,'unix'), snorm(emfisis.mag), $
        limit = {ytitle:'|B|!C(nT)'}
    get_data, tvar, uts, dat
    idx = where(uts ge utr0[0] and uts le utr0[1])
    store_data, tvar, uts[idx], dat[idx,*]
endforeach

ofn = shomedir()+'/fig_bmag.pdf'
if longutr then begin
    tutr = utr0
    ofn = shomedir()+'/fig_bmag_long.pdf'
endif else tutr = utr

;ofn = 0
sgopen, ofn, xsize = 8.5, ysize = 6, /inch

vars = ['bmag','vsc']
vars = ['rbspa_'+vars,'rbspb_'+vars]
nvar = n_elements(vars)
poss = sgcalcpos(nvar)

tplot, vars, trange = tutr, position = poss, vlab_margin = 10
sgclose



stop


; **** load mageis data and plot them.
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        mageis = sread_rbsp_mageis_l3(utr0, probes = tprobe)
        uts = sfmepoch(mageis.epoch,'unix')
        fedu = total(mageis.fedu,3,/nan)
        store_data, pre0+'e_en_mageis', uts, fedu, mageis.fedu_energy, limits = $
            {spec:1, no_interp:1, ylog:1, zlog:1, yrange:[4e1,4e3], ystyle:1, $
            ytitle:'MAGEIS!CRBSP-'+strupcase(tprobe)+'!Ce- energy!C(keV)', $
            zrange:[1e-2,1e5], ztitle:'(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'}
    endforeach
    
    
    
    ofn = shomedir()+'/fig_mageis.pdf'
    if longutr then begin
        tutr = utr0
        ofn = shomedir()+'/fig_mageis_long.pdf'
    endif else tutr = utr
    
    ;ofn = 0
    sgopen, ofn, xsize = 8.5, ysize = 3, /inch
    
    vars = ['rbspa_','rbspb_']+'e_en_mageis'
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar)
    
    tplot, vars, trange = tutr, position = poss, vlab_margin = 10
    sgclose



; load spice kernel.
defsysv,'!rbsp_spice', exists=flag
if flag eq 0 then rbsp_load_spice_kernels, trange = utr0


; load AE.
omni = sread_omni(omniutr)
store_data, 'ae', sfmepoch(omni.epoch,'unix'), omni.ae_index, limits = $
    {ytitle:'AE!C(nT)'}
store_data, 'symh', sfmepoch(omni.epoch,'unix'), omni.sym_h, limits = $
    {ytitle:'Sym/H!C(nT)'}


; load hope data.
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        vars = ['PITCH_ANGLE',$
            'Epoch_Ele','HOPE_ENERGY_Ele',$
            'Epoch_Ion','HOPE_ENERGY_Ion',$
            'FEDU','FPDU']
        hopel2 = sread_rbsp_hope_l3(utr0, probes = tprobe, vars = vars)
        if size(hopel2,/type) ne 8 then return
        uts = sfmepoch(hopel2.epoch_ion,'unix')
        store_data, pre0+'h_en', uts, total(hopel2.fpdu,3,/nan), hopel2.hope_energy_ion
        store_data, pre0+'h_pa', uts, total(hopel2.fpdu,2,/nan), hopel2.pitch_angle
        uts = sfmepoch(hopel2.epoch_ele,'unix')
        store_data, pre0+'e_en', uts, total(hopel2.fedu,3,/nan), hopel2.hope_energy_ele
        store_data, pre0+'e_pa', uts, total(hopel2.fedu,2,/nan), hopel2.pitch_angle
    
        vars = pre0+['h_en','h_pa','e_en','e_pa']
        options, vars, 'spec', 1
        options, vars, 'no_interp', 1
        options, vars, 'zlog', 1
        options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'
    
        vars = pre0+['h_en','e_en']
        options, vars, 'ylog', 1
        options, vars, 'yrange', [10,4e4]
    
        vars = pre0+['e_en','e_pa']
        options, vars, 'zrange', [1e4,1e10]
        vars = pre0+['h_en','h_pa']
        options, vars, 'zrange', [1e4,1e7]
    
        vars = pre0+['h_pa','e_pa']
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
    
        hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
        if size(hopemom,/type) ne 8 then return
        uts = sfmepoch(hopemom.epoch_ele,'unix')
        store_data, pre0+'n', uts, hopemom.dens_e_200, $
            limits = {ytitle:'N!Ie!N!C(cm!E-3!N)', ylog:1, constant:1}
        store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
            limits = {ytitle:'T!I!N!C(eV)', ylog:1, colors:[6,0], labels:['Tperp','Tpara'], constant:1000}
        tvar = pre0+'t'
        get_data, tvar, uts, dat
        idx = where(dat eq 1e20, cnt)
        if cnt ne 0 then dat[idx] = !values.d_nan
        store_data, tvar, uts, dat
    
    
        rbsp_preprocess_efield, utr, probes = tprobe, /no_spice_load
    
    endforeach
endif




vars = 'ae'
options, vars, 'labels', 'AE'
options, vars, 'yrange', [0,2000]
options, vars, 'yticks', 2
options, vars, 'yminor', 10

vars = ['rbspa_','rbspb_']+'n'
options, vars, 'ytitle', 'n!De!N!C(cm!U-3!N)'
options, vars, 'labels', 'n!De!N>200eV'

vars = ['rbspa_','rbspb_']+'vsc'
options, vars, 'yrange', [-30,0]
options, vars, 'yticks', 3

vars = ['rbspa_','rbspb_']+'h_en'
options, vars, 'yticks', 2
options, vars, 'ytickv', [1e2,1e3,1e4]
options, vars, 'yminor', 10
options, vars, 'ystyle', 1
options, vars, 'yrange', [40,4e4]


vars = ['rbspa_','rbspb_']+'e_en'
options, vars, 'yticks', 2
options, vars, 'ytickv', [1e2,1e3,1e4]
options, vars, 'yminor', 10
options, vars, 'ystyle', 1
options, vars, 'yrange', [40,4e4]
options, vars, 'zrange', [1e4,1e10]
options, vars, 'zticks', 3


vars = ['h_en','e_en','n','vsc']
nvar = n_elements(vars)
labs = ['mlat','lshell','mlt']

ofn = shomedir()+'/fig_overview.pdf'
if longutr then begin
    tutr = utr0
    ofn = shomedir()+'/fig_overview_long.pdf'
endif else tutr = utr
;ofn = 0
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43

poss = sgcalcpos(nvar*2+3+1)
tplot, 'ae', trange = tutr, position = poss[*,11], /noerase, vlab_margin = 10

pre0 = 'rbspa_'
tplot, pre0+vars, var_label = pre0+labs, trange = tutr, $
    position = poss[*,0:-1+nvar], /noerase, title = 'RBSP-A', /nouttick, vlab_margin = 10

tpos = poss[*,5:4+nvar]
for i = 0, nvar-1 do tpos[[1,3],i]-= 0.04
pre0 = 'rbspb_'
tplot, pre0+vars, var_label = pre0+labs, trange = tutr, $
    position = tpos, /noerase, title = 'RBSP-B', /nouttick, vlab_margin = 10

sgclose







end