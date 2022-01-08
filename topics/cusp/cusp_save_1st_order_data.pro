;+
; Read and save first order data, including
;   dE FAC, dB FAC, ion, ele KE flux, position and mapping info.
;-
pro cusp_save_1st_order_data, eventid, test=test, $
    no_plot=noplot, save_data=save_data, reload=reload
    
    
    if n_elements(eventid) eq 0 then message, 'no event id ...'
    
;---constant.
    re = 6378d & re1 = 1d/re
    va0 = 22d
    rad = !dpi/180d
    deg = 180d/!dpi
    
    
;---settings.
    ; rootdir to save the data and plots.
    rootdir = sdiskdir('GoogleDrive')+'/My Drive/works'
    if file_test(rootdir) eq 0 then rootdir = shomedir()

    ; dir to save plots and data.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'
    if ~file_test(figdir,/directory) then file_mkdir, figdir
    if ~file_test(datdir,/directory) then file_mkdir, datdir

    ; to prevent overwriting useful data without confirm.
    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif
    
    ; log file contains all the conjunc events.
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'


;---plot settings and generate plots to disk.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    falabs = 'fa_'+['ilat','mlt','dis']
    polabs = 'po_'+['ilat','mlt','dis']
    rgb = [6,4,2] & red = 6
    ct = 43
    charsz = 1
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'constant', 0
    time_stamp, /off
    
    

;---basic info for loading data.
; [po,fa]plottr. time range for plots.
; [po,fa]cusptr. cusp entry and exit times.

    ; read basic info from list of conjunc events.
    loginfo = cusp_read_conjun_list(logfile, event=eventid)
    if size(loginfo,/type) ne 8 then message, 'no id found ...'

    poinfo = loginfo.polar
    fainfo = loginfo.fast
    
    ; plot time range.
    poplotutr = poinfo.plot_time
    faplotutr = fainfo.plot_time
    if poplotutr[1] lt poplotutr[0] then poplotutr[1]+= 86400d
    if faplotutr[1] lt faplotutr[0] then faplotutr[1]+= 86400d

    ; cusp time range.
    pocusputr = poinfo.cusp_time
    facusputr = fainfo.cusp_time
    if pocusputr[1] lt pocusputr[0] then pocusputr[1]+= 86400d
    if facusputr[1] lt facusputr[0] then facusputr[1]+= 86400d



;---calculate data and save to tplot.
; [po,fa]_[ion,ele]_keflux
; [po,fa]_[de,db]_fac
; [po,fa]_map_coef

    ; KE flux, mapping coef.
    store_data, '*', /delete
    
    
    ;---fast particle and position data.
    ; fa_pos_gsm, fa_ele_n.
    faorb = string(loginfo.fast.orbit,format='(I05)')
    fn = rootdir+'/data/cusp/fa_sdt_esa_'+eventid+'_'+faorb+'.tplot'
    tplot_restore, filename = fn
    stplot_renew, 'ele_n', newname='fa_ele_n'
    get_data, 'fa_pos', uts, posgei
    ets = stoepoch(uts, 'unix')
    posgsm = sgse2gsm(sgei2gse(posgei, ets), ets)
    store_data, 'fa_pos_gsm', uts, posgsm*re1, $
        limits={colors:rgb,labels:['x','y','z']}


    ;---read kinetic energy flux.
    ; po_[ion,ele]_keflux.
    fn = rootdir+'/works/cusp/cusp list conjun/'+eventid+'/'+ $
        eventid+'_ke_special.svg'
    fn1 = rootdir+'/data/cusp/'+eventid+'_ke_special.svg'
    if file_test(fn1) eq 0 then file_copy, fn, fn1
    yrange = poinfo.keflux  ; [ele,ion].
    type = poinfo.ketype    ; etp,p2p,ptp.
    case type of
        'etp':ex = {epstopdf:1}
        'ptp':ex = {pstopdf:1}
        'p2p':ex = {ps2pdf:1}
    endcase
    polar_read_ke_flux, fn1, poplotutr, yrange[0]*[-1,1], $
        yrange[1]*[-1,1], _extra = ex

    ; fa_[ion,ele]_keflux.
    fn = rootdir+'/data/cusp/fa_sdt_esa_'+eventid+'_'+ $
        string(fainfo.orbit,format='(I05)')+'.tplot'
    tplot_restore, filename = fn
    stplot_renew, 'ele_eflux', newname='fa_ele_keflux', /delete
    stplot_renew, 'ion_eflux', newname='fa_ion_keflux', /delete
    vars = ['ion_*','ele_*']
    store_data, vars, /delete
    ; despike for ion ke. !!! need better despike algorithm.
    get_data, 'fa_ion_keflux', t0, dat
    idx = where(abs(dat) le 1, cnt)
    store_data, 'fa_ion_keflux', t0, interpol(dat[idx],t0[idx],t0)

    ;---read field data.
    ; [po,fa]_[de,db]_fac, original field.
    pre0 = 'po_'
    fn = rootdir+'/data/cusp/po_sdt_fld_'+eventid+'.sdt'
    ; preprocess field, remove spike, separate model B and dB, etc.
    polar_sdt_prep_poynting_flux, fn, e56=poinfo.e56, tr=poplotutr, $
        eventid=eventid, cusptr=pocusputr, noplot=noplot, $
        orootdir=figdir+'/'+eventid, titpre='Event ID: '+eventid+'    '
    vars = pre0+['b0_spc','spc2fac','db_spc','de_spc','mlat']
    store_data, vars, /delete

    pre0 = 'fa_'
    fn = rootdir+'/data/cusp/fa_sdt_fld_'+eventid+'_'+ $
        string(fainfo.orbit,format='(I05)')+'.tplot'
    if file_search(fn) eq '' then $
        fn = rootdir+'/data/cusp/fa_sdt_fld_'+eventid+'_'+$
        string(fainfo.orbit,format='(I05)')+'.tplot'
    fast_sdt_prep_poynting_flux, fn, $  ; plot time add padding time.
        trplot = faplotutr+(faplotutr[1]-faplotutr[0])*[-1,1]
    vars = [pre0+['pos','vel','alt','b0_gei','b0_gei'],'alt']
    store_data, vars, /delete
    
        
        
    ;---map coef.
    ; [po,fa]_fpt_[mlat,mlt]
    ; [po,fa]_pos_gsm
    ; [po,fa]_map_coef
    model = 't89'
    sats = ['po','fa']
    foreach tsat, sats do begin
        pre0 = tsat+'_'
        scalc_map_coef, pre0+'pos_gsm', pre0+'b', coord='gsm', $
            /igrf, pre=pre0, model=model
    endforeach
        


;---structure to hold useful info.
    scinfo = cusp_scinfo()
    
    scinfo.id = eventid
    scinfo.model = model


    ; polar.
    pre0 = 'po_'
    tinfo = scinfo.polar
    tloginfo = loginfo.polar
    cusputr = pocusputr
    ; time range.
    tinfo.plot_time = poplotutr
    tinfo.cusp_time = pocusputr
    ; values at the entry and exit.
    vars = pre0+['dis','mlt','ilat','fpt_mlt','fpt_mlat']
    nvar = n_elements(vars)
    tval = dblarr(nvar,2)
    for i=0, nvar-1 do begin
        get_data, vars[i], uts, dat
        tval[i,*] = interpol(dat, uts, cusputr)
    endfor
    tinfo.cusp_dis = reform(tval[0,*])
    tinfo.cusp_mlt = reform(tval[1,*])
    tinfo.cusp_ilat = reform(tval[2,*])
    tinfo.cusp_fpt_mlt = reform(tval[3,*])
    tinfo.cusp_fpt_mlat = reform(tval[4,*])
    ; characteristic values within cusp.
    vars = pre0+['b','dis','mlt','ilat']
    nvar = n_elements(vars)
    tval = dblarr(nvar)
    for i=0, nvar-1 do begin
        get_data, vars[i], uts, dat
        idx = where(uts ge cusputr[0] and uts le cusputr[1])
        tval[i] = median(dat[idx])
    endfor
    tinfo.b0 = tval[0]
    tinfo.dis0 = tval[1]
    tinfo.mlt0 = tval[2]
    tinfo.ilat0 = tval[3]
    ; data rate.
    get_data, pre0+'de_fac', uts
    tinfo.dr0 = sdatarate(uts)
    ; density.
    tinfo.n0 = max(tloginfo.n)
    ; derived quantities.
    tinfo.va = tinfo.b0/sqrt(tinfo.n0)*va0
    tinfo.vsc = (tinfo.cusp_ilat[1]-tinfo.cusp_ilat[0])*$
        rad*tinfo.dis0/(tinfo.cusp_time[1]-tinfo.cusp_time[0])*re
    tinfo.vdir = tinfo.cusp_ilat[1] gt tinfo.cusp_ilat[0]
    ; save values.
    scinfo.polar = tinfo

    
    

    ; fast.
    pre0 = 'fa_'
    tinfo = scinfo.fast
    tloginfo = loginfo.fast
    cusputr = facusputr
    ; time range.
    tinfo.plot_time = faplotutr
    tinfo.cusp_time = facusputr
    ; values at the entry and exit.
    vars = pre0+['dis','mlt','ilat','fpt_mlt','fpt_mlat']
    nvar = n_elements(vars)
    tval = dblarr(nvar,2)
    for i=0, nvar-1 do begin
        get_data, vars[i], uts, dat
        tval[i,*] = interpol(dat, uts, cusputr)
    endfor
    tinfo.cusp_dis = reform(tval[0,*])
    tinfo.cusp_mlt = reform(tval[1,*])
    tinfo.cusp_ilat = reform(tval[2,*])
    tinfo.cusp_fpt_mlt = reform(tval[3,*])
    tinfo.cusp_fpt_mlat = reform(tval[4,*])
    ; characteristic values within cusp.
    vars = pre0+['b','dis','mlt','ilat']
    nvar = n_elements(vars)
    tval = dblarr(nvar)
    for i=0, nvar-1 do begin
        get_data, vars[i], uts, dat
        idx = where(uts ge cusputr[0] and uts le cusputr[1])
        tval[i] = median(dat[idx])
    endfor
    tinfo.b0 = tval[0]
    tinfo.dis0 = tval[1]
    tinfo.mlt0 = tval[2]
    tinfo.ilat0 = tval[3]
    ; data rate.
    get_data, pre0+'de_fac', uts
    tinfo.dr0 = sdatarate(uts)
    ; density.
    get_data, pre0+'ele_n', uts, dat
    idx = where(uts ge cusputr[0] and uts le cusputr[1])
    tinfo.n0 = max(dat[idx])
    ; derived quantities.
    tinfo.va = tinfo.b0/sqrt(tinfo.n0)*va0
    tinfo.vsc = (tinfo.cusp_ilat[1]-tinfo.cusp_ilat[0])*$
        rad*tinfo.dis0/(tinfo.cusp_time[1]-tinfo.cusp_time[0])*re
    tinfo.vdir = tinfo.cusp_ilat[1] gt tinfo.cusp_ilat[0]
    ; save values.
    scinfo.fast = tinfo
    
    
    ; scinfo.
    tutr = minmax([pocusputr,facusputr])
    omni = sread_omni(tutr)
    uts = sfmepoch(omni.epoch, 'unix')
    idx = where(uts ge tutr[0] and uts le tutr[1])
    scinfo.ae = max((omni.ae_index)[idx])
    scinfo.dst = min((omni.sym_h)[idx])
    scinfo.fast_orbit = double(faorb)
    
    
    
    
    ;---calc ratio of O/H.

    ; polar.
    pre0 = 'po_'
    tinfo = scinfo.polar
    plotutr = tinfo.plot_time
    cusputr = tinfo.cusp_time
    
    dat = sread_polar_timas(plotutr, type='h0')
    if size(dat,/type) ne 8 then begin
        tinfo.r_oh_cusp = 0d
        tinfo.r_oh_pcap = 0d
    endif else begin
        uts = sfmepoch(dat.epoch_h, 'unix')
        utidx = where(uts ge plotutr[0] and uts le plotutr[1], nrec)
        if nrec eq 0 then begin
            tinfo.r_oh_cusp = 0d
            tinfo.r_oh_pcap = 0d
        endif else begin
            uts = uts[utidx]
            tsz = size(dat.flux_h,/dimensions)
            
            ens = dat.energy        ; energy bins, in eV.
            nen = n_elements(ens)
            if nen ne tsz[1] then begin ; some data has energy bin twice.
                ens = ens[0:tsz[1]-1]
                nen = n_elements(ens)
            endif
            dens = [ens[0],ens[1:nen-1]-ens[0:nen-2]]; dE at each energy bin, in eV.
            
            
            pas = dat.angle         ; pitch angle, in deg.
            npa = n_elements(pas)
            if npa ne tsz[2] then begin
                pas = pas[0:tsz[2]-1]
                npa = n_elements(pas)
            endif
            dpas = dblarr(npa)+180d/npa
            dpas = dpas*rad
            
            ospec = (dat.flux_o)[utidx,*,*]
            for i=0, nen-1 do ospec[*,i,*] *= (1d-3*dens[i]) ; convert to '#/(cm!U2!N-s-sr)'.
            hspec = (dat.flux_h)[utidx,*,*]
            for i=0, nen-1 do hspec[*,i,*] *= (1d-3*dens[i]) ; convert to '#/(cm!U2!N-s-sr)'.
            
            
            ;---oxygen.
            tspec = ospec
            pre1 = 'po_o_'
            tmass = omass
            
            idx = where(finite(tspec,/nan),cnt)
            if cnt ne 0 then tspec[idx] = 0
            
            tvar = pre1+'n'
            ns = dblarr(nrec)
            for i = 0, nrec-1 do begin
                nflux = reform(tspec[i,*,*])    ; in #/s-cm^2-sr.
                vs = sqrt(2*ens/tmass)*1e3      ; energy converted to velocity, in km/s.
                for j = 0, nen-1 do begin       ; integrate over the energy bins.
                    for k = 0, npa-1 do $       ; integrate each pitch angle.
                        ns[i]+= nflux[j,k]/vs[j]* $
                        2*!dpi*sin(pas[k])*dpas[k]
                endfor
            endfor
            store_data, pre1+'n', uts, ns
                
                
            ;---proton.
            tspec = hspec
            pre1 = 'po_h_'
            tmass = hmass
            
            idx = where(finite(tspec,/nan),cnt)
            if cnt ne 0 then tspec[idx] = 0
            
            tvar = pre1+'n'
            ns = dblarr(nrec)
            for i = 0, nrec-1 do begin
                nflux = reform(tspec[i,*,*])    ; in #/s-cm^2-sr.
                vs = sqrt(2*ens/tmass)*1e3      ; energy converted to velocity, in km/s.
                for j = 0, nen-1 do begin       ; integrate over the energy bins.
                    for k = 0, npa-1 do $       ; integrate each pitch angle.
                        ns[i]+= nflux[j,k]/vs[j]* $
                        2*!dpi*sin(pas[k])*dpas[k]
                endfor
            endfor
            store_data, pre1+'n', uts, ns
                
                
            ;---calc numbers.
            ; proton.
            get_data, pre0+'h_n', uts, dat
            idx = where(dat ne 0)
            dat = interpol(dat[idx],uts[idx],uts)
            
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            nh1 = mean(dat[idx])
            
            if tinfo.vdir eq 0 then begin
                idx = where(uts le cusputr[0])
                nh2 = mean(dat[idx])
            endif else begin
                idx = where(uts ge cusputr[1])
                nh2 = mean(dat[idx])
            endelse
            
            ; oxygen.
            get_data, pre0+'o_n', uts, dat
            idx = where(dat ne 0)
            dat = interpol(dat[idx],uts[idx],uts)
            
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            no1 = mean(dat[idx])
            
            if tinfo.vdir eq 0 then begin
                idx = where(uts le cusputr[0])
                no2 = mean(dat[idx])
            endif else begin
                idx = where(uts ge cusputr[1])
                no2 = mean(dat[idx])
            endelse
            
            
            tinfo.r_oh_cusp = double(no1/nh1)
            tinfo.r_oh_pcap = double(no2/nh2)
        endelse
    endelse
    scinfo.polar = tinfo
    
    
    ; fast.
    pre0 = 'fa_'
    tinfo = scinfo.fast
    plotutr = tinfo.plot_time
    cusputr = tinfo.cusp_time

    dat = sread_fast_tms(plotutr)
    if size(dat,/type) ne 8 then begin
        tinfo.r_oh_cusp = 0d
        tinfo.r_oh_pcap = 0d
    endif else begin    ; have data in that day.
        uts = sfmepoch(dat.epoch, 'unix')
        utidx = where(uts ge plotutr[0] and uts le plotutr[1], nrec)
        if nrec eq 0 then begin
            tinfo.r_oh_cusp = 0d
            tinfo.r_oh_pcap = 0d
        endif else begin
            uts = uts[utidx]
                
            ospec = (dat.o_)[utidx,*]
            o_ens = (dat.o__en)[utidx,*]
            hspec = (dat.h_)[utidx,*]
            h_ens = (dat.h__en)[utidx,*]
                
                
            ;---oxygen.
            tspec = ospec
            t_ens = o_ens
            pre1 = 'fa_o_'
            tmass = omass
            
            idx = where(finite(tspec,/nan),cnt)
            if cnt ne 0 then tspec[idx] = 0
                
            tvar = pre1+'n'
            ns = dblarr(nrec)
            for i = 0, nrec-1 do begin
                nflux = reform(tspec[i,*])      ; in eV/s-cm^2-sr-eV.
                tens = reform(t_ens[i,*])       ; in eV
                idx = sort(tens)
                nflux = nflux[idx]
                tens = tens[idx]
                nten = n_elements(tens)
                tden = [tens[0],tens[1:nten-1]-tens[0:nten-2]]
                vs = sqrt(2*tens/tmass)*1e3      ; energy converted to velocity, in km/s?
                for j = 0, nten-1 do begin      ; integrate over the energy bins.
                    ns[i]+= nflux[j]/tens[j]*tden[j]/vs[j]*4*!dpi
                endfor
            endfor
            store_data, pre1+'n', uts, ns
                    
                    
            ;---proton.
            tspec = hspec
            t_ens = h_ens
            pre1 = 'fa_h_'
            tmass = hmass
            
            idx = where(finite(tspec,/nan),cnt)
            if cnt ne 0 then tspec[idx] = 0
            
            tvar = pre1+'n'
            ns = dblarr(nrec)
            for i = 0, nrec-1 do begin
                nflux = reform(tspec[i,*])      ; in eV/s-cm^2-sr-eV.
                tens = reform(t_ens[i,*])       ; in eV
                idx = sort(tens)
                nflux = nflux[idx]
                tens = tens[idx]
                nten = n_elements(tens)
                tden = [tens[0],tens[1:nten-1]-tens[0:nten-2]]
                vs = sqrt(2*tens/tmass)*1e3     ; energy converted to velocity, in km/s?
                for j = 0, nten-1 do begin      ; integrate over the energy bins.
                    ns[i]+= nflux[j]/tens[j]*tden[j]/vs[j]*4*!dpi
                endfor
            endfor
            store_data, pre1+'n', uts, ns
                    
                    
            ;---calc numbers.
            ; proton.
            get_data, pre0+'h_n', uts, dat
            idx = where(dat ne 0)
            dat = interpol(dat[idx],uts[idx],uts)
            
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            nh1 = mean(dat[idx])
            
            if tinfo.vdir eq 0 then begin
                idx = where(uts le cusputr[0])
                nh2 = mean(dat[idx],/nan)
            endif else begin
                idx = where(uts ge cusputr[1])
                nh2 = mean(dat[idx],/nan)
            endelse
            
            ; oxygen.
            get_data, pre0+'o_n', uts, dat
            idx = where(dat ne 0)
            dat = interpol(dat[idx],uts[idx],uts)
            
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            no1 = mean(dat[idx])
            
            if tinfo.vdir eq 0 then begin
                idx = where(uts le cusputr[0])
                no2 = mean(dat[idx],/nan)
            endif else begin
                idx = where(uts ge cusputr[1])
                no2 = mean(dat[idx],/nan)
            endelse
                
            tinfo.r_oh_cusp = double(no1/nh1)
            tinfo.r_oh_pcap = double(no2/nh2)
        endelse
    endelse
    scinfo.fast = tinfo
    

    infofn = datdir+'/'+eventid+'_scinfo.tplot'
    store_data, 'scinfo', tutr, scinfo
    tplot_save, 'scinfo', filename = infofn




;---save data to disk.
    ; save labels, de_fac, db_fac, ele_keflux, ion_keflux, in situ.
    ; save event info to (1) include other info and (2) info for calc
    ; quantities are not saved (e.g., pflux).
    if keyword_set(save_data) then begin
        vars = ['ilat','mlt','dis','fpt_mlat','fpt_mlt','map_coef','pos_gsm', $
            'de_fac','db_fac','ele_keflux','ion_keflux']
        vars = ['po_'+vars,'fa_'+vars]
        ofn = datdir+'/'+eventid+'_1st_order_data.tplot'
        tplot_save, vars, filename = ofn
    endif
    while !d.window ne -1 do wdelete, !d.window



    if keyword_set(noplot) then return


    options, ['fa_','po_']+'dis', 'ytitle', 'R'
    ct = 43
    red = 6
    

;---plot 1: field preprocessing (several plots).
    ; done by polar_sdt_prep_poynting_flux.




;---plot 1, FAST ESA.
    mu = '!9'+string(109b)+'!X'
    gamma = '!9'+string(71b)+'!X'
    para = '||'
    perp = '!9'+string(94b)+'!X'
    orb = string(loginfo.fast.orbit,format='(I05)')
    fn = rootdir+'/data/cusp/fa_sdt_esa_'+eventid+'_'+orb+'.tplot'
    tplot_restore, filename = fn
    get_data, 'fa_dis', t0
    options, ['ion_n','ion_p','ele_n','ele_p'], 'ylog', 1
    vars = ['ion_p','ele_p']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, tmp
        store_data, vars[i], t0, tmp*1e-8
        options, vars[i], 'ytitle', '(ba)!Ccgs'
    endfor
    coef = 1.6e-9   ; convert # cm^-2 s^-1 to uA m^-2.
    pre = 'ion_'
    vars = pre+['para_spec','perp_spec','anti_spec',$
        'j','eflux','n','p','en_spec','pa_spec']
    fn = shomedir()+'/'+eventid+'_esa.pdf'
    get_data, 'ele_n', t0, elen
    get_data, 'ion_n', tmp, ionn
    ionn = interpol(ionn, tmp, t0)
    store_data, 'density', t0, [[ionn],[elen]], limits = $
        {labels:['Ni','Ne'], colors:[-1,6], ytitle:'(cm!U-3!N)',ylog:1}
    get_data, 'ele_j', t0, elej
    get_data, 'ion_j', tmp, ionj
    ionj = interpol(ionj, tmp, t0)
    store_data, 'nflux', t0, [[ionj],[elej]], limits = $
        {labels:['Ji','Je'], colors:[-1,6], ytitle:'(!9'+mu+'A/m!U2!N)'}
    vars = ['ion_en_spec','ion_pa_spec','nflux','fa_ion_keflux','fa_ele_keflux', $
        'density','ele_en_spec']
    ofn = figdir+'/'+eventid+'/'+eventid+'_esa.pdf'

    sgopen, ofn, xsize = 6, ysize = 8, /inch
    device, decomposed = 0
    loadct2, ct
    titl = 'Event ID: '+eventid+'    FAST ESA'
    tplot, vars, trange = faplotutr, var_label = falabs, title = titl
    timebar, facusputr, color = red, thick = 2

    sgclose
    

;---plot 2, footpoint.
    alltr = []
    satnames = ['polar (*,b)','fast (+,r)']
    satcolors = [sgcolor('blue'),sgcolor('red')]
    satsyms = [2,1]     ; [*,+].
    pos = [0.15,0.1,0.85,0.8]
    thick = 2
    ofn = figdir+'/'+eventid+'/'+eventid+'_footpoint.pdf'
    sgopen, ofn, xsize = 6, ysize = 3, /inch
    sgtruecolor
    ; determine hemisphere.
    get_data, 'po_ilat', t0, poilat
    if interpol(poilat, t0, pocusputr[0]) ge 0 then begin   ; north hem.
        sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
            ytickv = [50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
            xtickpos = 47
    endif else begin                                        ; south hem.
        sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
            ytickv = -[50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
            xtickpos = -47, yrange = [-50,-90]
    endelse
    ; use cusp time, expand 2.5 cusp time on both sides.
    ; polar ilat, mlt.
    j = 0
    dt = 600        ; 10 min.
    get_data, 'po_dis', t0, tmp
    if min(interpol(tmp, t0, pocusputr)) le 2.5 then dt = 120
    tmp = pocusputr
    ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
    ttr = ttr-(ttr mod dt) & ttr[1]+= dt
    alltr = [min([ttr,alltr]),max([ttr,alltr])]
    tuts = smkarthm(ttr[0], ttr[1], dt, 'dx')
    get_data, 'po_ilat', t0, tmp
    tilat = interpol(tmp, t0, tuts)
    get_data, 'po_mlt', t0, tmp
    tmlt = interpol(tmp, t0, tuts)*15
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j]
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j], $
        psym = satsyms[j], symsize = 0.4
    arrow, tmlt[-2], tilat[-2], tmlt[-1], tilat[-1], /data, $
        color = satcolors[j], /solid, thick = thick
    plots, interpol(tmlt, tuts, pocusputr), interpol(tilat, tuts, pocusputr), $
        color = satcolors[j], thick = thick*5
    get_data, 'po_dis', t0, tmp
    tmp = string(interpol(tmp, t0, pocusputr[0]), format='(F3.1)')
    xyouts, 0.6, 0.25, /normal, 'Polar: *, R = '+tmp+' Re', color = satcolors[j]
    ; fast ilat, mlt.
    j = 1
    dt = 60        ; 1 min.
    tmp = facusputr
    ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
    ttr = ttr-(ttr mod dt) & ttr[1]+= dt
    alltr = [min([ttr,alltr]),max([ttr,alltr])]
    tuts = smkarthm(ttr[0], ttr[1], dt, 'dx')
    get_data, 'fa_ilat', t0, tmp
    tilat = interpol(tmp, t0, tuts)
    get_data, 'fa_mlt', t0, tmp
    tmlt = interpol(tmp, t0, tuts)*15
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j]
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j], $
        psym = satsyms[j], symsize = 0.4
    arrow, tmlt[-2], tilat[-2], tmlt[-1], tilat[-1], /data, $
        color = satcolors[j], /solid, thick = thick
    plots, interpol(tmlt, tuts, facusputr), interpol(tilat, tuts, facusputr), $
        color = satcolors[j], thick = thick*5
    tmp = string((facusputr[0]-pocusputr[0])/3600d, format='(F4.1)')
    tmp = strtrim(tmp,2)
    xyouts, 0.6, 0.2, /normal, 'FAST: +, dT = '+tmp+' hr', color = satcolors[j]

    xyouts, 0.15, 0.9, /normal, 'Event ID: '+eventid+'    Footprints of Polar and FAST', charsize = 1.25
    sgclose


;---plot 3: omni, polar, fast overview.
    ofn = figdir+'/'+eventid+'/'+eventid+'_overview.pdf'
    sgopen, ofn, xsize = 6, ysize = 8, /inch
    device, decomposed = 0
    loadct2, ct
    titl = 'Event ID: '+eventid+'    Summary plot of OMNI, Polar, and FAST'
    plot_polar_fast_summary, stoepoch(alltr,'unix'), /no_delete, title = titl
    sgclose
    
end

id = '1998_0925_05'
id = '1999_1010_15'
cusp_save_1st_order_data, id, /save_data, /no_plot
end
