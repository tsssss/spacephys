;+
; Load Polar Hydra energy fluxes; Timas density, nflux, eflux, eflux>500eV;
; Load Fast ESA energy fluxes; Teams density;
; For both spacecraft, load electric and magnetic field,
; then calculate Poynting flux and filter it use the default filter.
;-

pro cusp_polar_save_data, eventid, test=test, no_plot=no_plot, reload=reload

    if n_elements(eventid) eq 0 then message, 'no event id ...'


;---Constants.
    re = 6378d & re1 = 1d/re
    va0 = 22d
    rad = !dpi/180d
    deg = 180d/!dpi
    secofday = 86400d

    omass = 16*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
    hmass = 1*(1.67e-27/1.6e-19)    ; E in eV, mass in kg.

    pres = ['po_']

;---Settings.
    ; rootdir to save the data and plots.
    rootdir = googledir()+'/works'
    if file_test(rootdir) eq 0 then message, 'Root directory not found ...'
    ; log file contains all the conjunc events.
    logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'


    ; dir to save plots and data.
    figdir = rootdir+'/works/cusp/cusp list polar 2-4Re'
    datdir = rootdir+'/data/cusp'
    if ~file_test(figdir,/directory) then file_mkdir, figdir
    if ~file_test(datdir,/directory) then file_mkdir, datdir
    datfn = datdir+'/'+eventid+'_polar_data.tplot'

    ; to prevent overwriting useful data without confirm.
    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif

    ; model used for mapping.
    model = 't89'


;---Plot settings.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    xyz = ['x','y','z']
    falabs = 'fa_'+['ilat','mlt','dis']
    polabs = 'po_'+['ilat','mlt','dis']
    rgb = [6,4,2]
    red = 6
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



;---Load basic variables.
    ; common variables.
    vars = ['ilat','mlt','dis','fpt_mlat','fpt_mlt','map_coef','pos_gsm', $
        'de_fac','db_fac','pf_fac','pfb_mor_spec', $
        'ele_eflux','ion_eflux','h_density','o_density','e_density','ion_nflux']
    common_vars = ['po_'+vars, 'fa_'+vars]

    ; polar only variables.
    po_vars = ['po_'+['h_nflux','o_nflux','h_eflux','o_eflux'],tnames('po_*_500eV')]

    ; info variable.
    info_var = 'event_info'

    ; basic vars.
    basic_vars = [info_var, common_vars, po_vars]

    load = 0
    if file_test(datfn) eq 0 then load = 1
    if keyword_set(reload) then load = 1
    if load eq 0 then begin
        if tnames(info_var) ne '' then begin
            get_data, info_var, tmp, eventinfo
            if eventinfo.id ne eventid then tplot_restore, filename=datfn
        endif else tplot_restore, filename=datfn
    endif else begin
        if file_test(datfn) eq 1 then file_delete, datfn
        store_data, basic_vars, /delete

        eventinfo = cusp_event_info()
        loginfo = cusp_polar_read_event_info(logfile, event=eventid)
        if size(loginfo,/type) ne 8 then message, 'no id found ...'

        poinfo = loginfo.polar
        poplotutr = poinfo.plot_time
        pocusputr = poinfo.cusp_time

        eventinfo.id = eventid
        eventinfo.model = model
        eventinfo.polar.plot_time = poplotutr
        eventinfo.polar.cusp_time = pocusputr

    ;---Load AE and Dst to init event info.
        tutr = minmax([pocusputr])
        omni = sread_omni(tutr)
        uts = sfmepoch(omni.epoch, 'unix')
        idx = where(uts ge tutr[0] and uts le tutr[1])
        eventinfo.ae = max((omni.ae_index)[idx])
        eventinfo.dst = min((omni.sym_h)[idx])


    ;---Read Polar data.
        ; polar ion and electron energy fluxes.
        ; po_[ion,ele]_eflux.
        polar_load_hydra_eflux, poplotutr, /ion
;        polar_load_hydra_eflux, poplotutr, /electron
        fn = rootdir+'/works/cusp/cusp list conjun/'+eventid+'/'+eventid+'_ke_special.svg'
        fn1 = rootdir+'/data/cusp/'+eventid+'_ke_special.svg'
        if file_test(fn1) eq 0 then file_copy, fn, fn1
        yrange = poinfo.keflux  ; [ele,ion].
        type = poinfo.ketype    ; etp,p2p,ptp.
        case type of
            'etp':ex = {epstopdf:1}
            'ptp':ex = {pstopdf:1}
            'p2p':ex = {ps2pdf:1}
        endcase
        polar_read_ke_flux, fn1, poplotutr, yrange[0]*[-1,1], yrange[1]*[-1,1], _extra = ex
        stplot_renew, 'po_ele_keflux', newname='po_ele_eflux', /delete
        stplot_renew, 'po_ion_keflux', newname='po_ion_eflux', /delete

        polar_read_hydra, poplotutr, 'moment', species='ion'


        hydroot = sdiskdir('Research')+'/sdata/opt_hydra/moment_data'
        pattern = time_string(pocusputr[0], tformat='YYYYMMDD')+'_hyd_mom_v*.cdf'
        fn = file_search(hydroot, pattern)
        cdf = scdfread(fn, ['epoch','n'])
        uts = sfmepoch(*cdf[0].value, 'unix')
        idx = where(uts ge poplotutr[0] and uts le poplotutr[1])
        dat = *cdf[1].value
        store_data, 'po_e_density', uts[idx], dat[idx], limits=$
            {ytitle:'(cm!U-3!N)', labels:'Hydra Ne', ylog:1}


        ; polar electric and magnetic fields.
        ; [po,fa]_[de,db]_fac, original field.
        pre0 = 'po_'
        fn = rootdir+'/data/cusp/po_sdt_fld_'+eventid+'.sdt'
        ; preprocess field, remove spike, separate model B and dB, etc.
        polar_sdt_prep_poynting_flux, fn, e56=poinfo.e56, tr=poplotutr, $
            eventid=eventid, cusptr=pocusputr, noplot=1, $
            orootdir=figdir+'/'+eventid, titpre='Event ID: '+eventid+'    '
        vars = pre0+['b0_spc','spc2fac','db_spc','de_spc','mlat']
        store_data, vars, /delete

        ; ion data.
        ; po_[h,o]_density, po_[h,o]_nflux, po_[h,o]_eflux.
        ; po_[h,o]_density_500eV, po_[h,o]_nflux_500eV, po_[h,o]_eflux_500eV.
        ; po_ion_nflux.
        polar_load_ion_moments, poplotutr
        polar_load_ion_moments, poplotutr, energy_range=[500,1e6]
        sys_add, 'po_o_nflux', 'po_h_nflux', to='po_ion_nflux'  ; overwrite hydra nflux with timas.


    ;---Mapping coefficient.
        ; [po,fa]_fpt_[mlat,mlt]
        ; [po,fa]_pos_gsm
        ; [po,fa]_map_coef
        foreach pre0, pres do $
            scalc_map_coef, pre0+'pos_gsm', pre0+'b', coord='gsm', /igrf, pre=pre0, model=model, altitude=eventinfo.h0

    ;---Fill in event info.
        foreach pre0, pres do begin
            if pre0 eq 'po_' then begin
                tinfo = eventinfo.polar
            endif else begin
                tinfo = eventinfo.fast
            endelse
            cusputr = tinfo.cusp_time
            plotutr = tinfo.plot_time

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
            get_data, pre0+'ion_eflux', uts
            tinfo.dr1 = sdatarate(uts)
            if tnames(pre0+'h_density') ne '' then begin
                get_data, pre0+'h_density', uts
                tinfo.dr2 = sdatarate(uts)
            endif

            ; density.
            get_data, pre0+'e_density', uts, tmp
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            tinfo.n_ele = max(tmp)

            if tnames(pre0+'h_density') eq '' then begin
                tinfo.n_ion = 0
                tinfo.r_no = 1
            endif else begin
                get_data, pre0+'h_density', uts, hden
                get_data, pre0+'o_density', uts, oden
                idx = where(uts ge cusputr[0] and uts le cusputr[1])
                tmp = hden[idx]+oden[idx]
                tinfo.n_ion = median(tmp)
                tinfo.r_no = median(oden[idx]/tmp)
            endelse

            ; Alfven speed derived from B and density.
            den = tinfo.n_ele
            tinfo.va = tinfo.b0/sqrt(den*(1+15*tinfo.r_no))*va0
            tinfo.va_h = tinfo.b0/sqrt(den)*va0
            tinfo.va_o = tinfo.b0/sqrt(den*16)*va0

            ; averaged s/c velocity and direction.
            tinfo.vsc = (tinfo.cusp_ilat[1]-tinfo.cusp_ilat[0])*$
                rad*tinfo.dis0/(tinfo.cusp_time[1]-tinfo.cusp_time[0])*re
            tinfo.vdir = tinfo.cusp_ilat[1] gt tinfo.cusp_ilat[0]


            ; correct for the east-west electric field.
            tinfo.pf_corr = cusp_correct_east_west_efield(pre0+'de_fac', pre0+'db_fac', cusputr)

            ; filters for Poynting flux.
            tinfo.cusp_duration = cusputr[1]-cusputr[0]
            get_data, pre0+'de_fac', uts
            tinfo.edge_duration = min([cusputr[0]-min(uts),max(uts)-cusputr[1]])
            s0 = 4*tinfo.dr0
            w0 = 6d
            s2t = 4*!dpi/(w0+sqrt(2+w0^2))
            tinfo.filters = [s0, min([tinfo.cusp_duration,tinfo.edge_duration])]*s2t    ; in sec, in period.

            ; scales for Poynting flux.
            tinfo.scaleinfo.s0 = s0     ; in sec, in scale.
            tinfo.scaleinfo.s1 = 0.5*(max(uts)-min(uts))
            tinfo.scaleinfo.dj = 1d/8
            ; the 'value' of pflux will be in period.
            stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scaleinfo, filter=tinfo.filters
            get_data, pre0+'pf_fac_mor_spec_'+string(eventinfo.pfidx+1,format='(I0)'), uts, dat, val
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            store_data, pre0+'pfb_mor_spec', uts[idx], dat[idx,*], 1d/val, $
                limits={spec:1, no_interp:1, ylog:1, ystyle:1, $
                yrange:minmax(1d/tinfo.filters), ytitle:'Freq (Hz)', $
                zrange:[-1,1]*max(abs(dat[idx,*]))*0.1}

            ; save values.
            if pre0 eq 'po_' then eventinfo.polar = tinfo else eventinfo.fast = tinfo
        endforeach

        eventinfo.hemisphere = (eventinfo.polar.ilat0 gt 0)? 1: -1

        store_data, info_var, eventid, eventinfo
        tplot_save, basic_vars, filename=datfn
    endelse


;---Map flux to ionosphere and convert parallel to earthward.
    vars = ['pfb','ele_eflux','ion_eflux','ion_nflux']
    mapped_vars = ['po_'+vars+'_map','fa_'+vars+'_map']
    load = 0
    foreach tvar, mapped_vars do begin
        if tnames(tvar) eq '' then begin
            print, tvar+' does not in tplot ...'
            load = 1
            break
        endif
    endforeach
    if keyword_set(reload) then load = 1
    if load eq 0 then begin
        if tnames(info_var) ne '' then begin
            get_data, info_var, tmp, eventinfo
            if eventinfo.id ne eventid then tplot_restore, filename=datfn
        endif else tplot_restore, filename=datfn
    endif else begin
        foreach pre0, pres do stplot_index, pre0+'pf_fac', eventinfo.pfidx, newname=pre0+'pfb'
        vars = ['pfb','ele_eflux','ion_eflux','ion_nflux']
        foreach pre0, pres do begin
            get_data, pre0+'map_coef', uts, c_map
            foreach tvar, vars do begin
                get_data, pre0+tvar, tuts, dat
                dat *= interpol(c_map, uts, tuts)
                if eventinfo.hemisphere eq -1 then dat = -dat
                store_data, pre0+tvar+'_map', tuts, dat
            endforeach

            tinfo = (pre0 eq 'po_')? eventinfo.polar: eventinfo.fast
            cusputr = tinfo.cusp_time

            ; Parallel Poynting flux vs total Poynting flux.
            get_data, pre0+'pf_fac', uts, dat
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            tinfo.r_pfb = total(dat[idx,2])/total(snorm(dat[idx,*]))

            ; ion ratio in terms of eflux and nflux.
            ; ion_eflux integration.
            get_data, pre0+'ion_eflux_map', uts, dat
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            dat = dat[idx]
            tinfo.r_eflux = total(dat)/total(abs(dat))
            tinfo.max_eflux = max(abs(dat))
            tinfo.kei_temporal = mean(dat)
            get_data, pre0+'fpt_mlat', tuts, mlat
            mlat = (interpol(mlat, tuts, uts))[idx]
            dlat = deriv(mlat)
            tinfo.kei_spatial = total(dlat*dat)/(mlat[-1]-mlat[0])

            get_data, pre0+'ion_nflux_map', uts, dat
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            dat = dat[idx]
            tinfo.r_nflux = total(dat)/total(abs(dat))
            tinfo.max_nflux = max(abs(dat))

            ; ele_eflux integration.
            get_data, pre0+'ele_eflux_map', uts, dat
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            dat = dat[idx]
            tinfo.kee_temporal = mean(dat)
            get_data, pre0+'fpt_mlat', tuts, mlat
            mlat = (interpol(mlat, tuts, uts))[idx]
            dlat = deriv(mlat)
            tinfo.kee_spatial = total(dlat*dat)/(mlat[-1]-mlat[0])

            ; pflux integration.
            get_data, pre0+'pfb_map', uts, dat
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            dat = dat[idx]
            tinfo.pfb_temporal = mean(dat)
            get_data, pre0+'fpt_mlat', tuts, mlat
            mlat = (interpol(mlat, tuts, uts))[idx]
            dlat = deriv(mlat)
            tinfo.pfb_spatial = total(dlat*dat)/(mlat[-1]-mlat[0])
;            if pre0 eq 'po_' then begin
;                filter = 1d/eventinfo.fast.filters
;                get_data, pre0+'pfb_mor_spec', uts, dat, val
;                idx = where(val ge min(filter) and val le max(filter), cnt)
;                dat = total(dat[*,idx], 2)
;                get_data, pre0+'map_coef', tuts, c_map
;                dat *= interpol(c_map, tuts, uts)
;                if eventinfo.hemisphere eq -1 then dat = -dat
;                idx = where(uts ge cusputr[0] and uts le cusputr[1])
;                dat = dat[idx]
;                tinfo.pfb_temporal = mean(dat)
;            endif

            ; the spatially integrated value.
            dlen = abs(tinfo.cusp_fpt_mlat[0]-tinfo.cusp_fpt_mlat[1])*rad*(eventinfo.h0+eventinfo.re)
            tinfo.pfb_int_spatial = tinfo.pfb_spatial*dlen
            tinfo.kei_int_spatial = tinfo.kei_spatial*dlen
            tinfo.kee_int_spatial = tinfo.kee_spatial*dlen
            dtime = abs(cusputr[0]-cusputr[1])
            tinfo.pfb_int_temporal = tinfo.pfb_temporal*dtime
            tinfo.kei_int_temporal = tinfo.kei_temporal*dtime
            tinfo.kee_int_temporal = tinfo.kee_temporal*dtime

            ; save values.
            if pre0 eq 'po_' then eventinfo.polar = tinfo else eventinfo.fast = tinfo
        endforeach

        store_data, 'event_info', eventid, eventinfo

        all_vars = [basic_vars, mapped_vars]
        tplot_save, all_vars, filename=datfn
    endelse

end

type = 'polar_2to4re_good'
ids = cusp_id_new(type)
;ids = '1999_0904_17'
;ids = '1998_0925_05'
;idx = where(ids eq '1998_1030_03')
;ids = ids[idx[0]:*]
foreach id, ids do cusp_polar_save_data, id, /reload
;cusp_combine_event_info, type
end