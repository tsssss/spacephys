;+
; loads RBSP burst data.
; RBSP-A VB1 in 4096 Hz, -B VB1 in 1024 Hz.
;-

pro _2013_0607_load_burst_data


; settings.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:01'])
autr = time_double(['2013-06-07/04:52','2013-06-07/04:58:30'])
butr = time_double(['2013-06-07/04:54:30','2013-06-07/05:01'])
probes = ['a','b']

; load burst data.
reload = 0
datfn = sdiskdir('GoogleDrive')+'/My Drive/works/works/psbl_de_32hz/2013_0607_b1_data.tplot'

; one time settings.
vars = ['vb1']
datvars = ['rbspa_'+vars,'rbspb_'+vars]
vars = ['u','v','w']+'_gsm'
posvars = ['rbspa_'+vars,'rbspb_'+vars]

allvars = [datvars,posvars]




device, decomposed = 0
loadct2, 43
!p.background = 255
!p.color = 0
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 3


load = 0
foreach tvar, datvars do if tnames(tvar) eq '' then load = 1

if load eq 0 then if file_test(datfn) ne 0 then begin
    store_data, '*', /delete
    tplot_restore, filename = datfn
endif

load = 0
foreach tvar, datvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load eq 1 then begin
    if file_test(datfn) ne 0 then begin
        store_data, '*', /delete
        tplot_restore, filename = datfn
    endif
endif


load = 0
if reload then load = 1
if load then begin
    defsysv,'!rbsp_spice', exists=flag
    if flag eq 0 then rbsp_load_spice_kernels, trange = utr0
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        tutr = (tprobe eq 'a')? autr: butr
        dr0 = (tprobe eq 'a')? 1d/4096: 1d/1024
        psbl_de_trim_b1_data, tutr, tprobe
        
        tvar = pre0+'vb1'
        get_data, tvar, tuts, dat
        uts = smkarthm(tutr[0], tutr[1], dr0, 'dx')
        dat = sinterpol(dat, tuts, uts)
        store_data, tvar, uts, dat
    endforeach
    
    options, 'rbspa_vb1', 'yrange', [-35,-5]
    options, 'rbspb_vb1', 'yrange', [-50,0]
    
    tplot_save, allvars, filename = datfn
endif



load = 0
foreach tvar, posvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load eq 1 then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        get_data, pre0+'vb1', tuts
        
        dr0 = 0.1   ; 110 points per spin.
        uts = smkarthm(min(tuts),max(tuts),dr0,'dx')
        
        ; calc epoch times for spice (different than the epoch in cdfs).
        tmp = time_string(uts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+uts-uts[0]
        
        scid = strupcase(pre0+'science')
        cspice_pxform, scid, 'GSM', tets, muvw2gsm
        
        ; muvw2gse[0,*,*] is u in GSM.
        ; muvw2gse[1,*,*] is v in GSM.
        ; muvw2gse[2,*,*] is w in GSM.
        store_data, pre0+'u_gsm', uts, transpose(reform(muvw2gsm[0,*,*])), $
            limits = {ytitle:'U GSM', labels:['x','y','z'], colors:[6,4,2]}
        store_data, pre0+'v_gsm', uts, transpose(reform(muvw2gsm[1,*,*])), $
            limits = {ytitle:'V GSM', labels:['x','y','z'], colors:[6,4,2]}
        store_data, pre0+'w_gsm', uts, transpose(reform(muvw2gsm[2,*,*])), $
            limits = {ytitle:'W GSM', labels:['x','y','z'], colors:[6,4,2]}
    endforeach
    
    tplot_save, allvars, filename = datfn
endif

end