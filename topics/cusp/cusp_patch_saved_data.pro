
pro cusp_patch_saved_data, eventid, $
    patch = patch

    if n_elements(eventid) eq 0 then message, 'no event id ...'
    
    if n_elements(patch) eq 0 then patch = 0
    case patch of
        0: begin        ; apply all patches.
            patch1 = 1
            end
        else: ; do nothing.
    endcase
    

    ;---settings.
        yticklen = -0.01d
        enyrng = [20d,20000]
        enzrng = [1e4,1e8]
        pazrng = [1e5,1e9]
        nyrng = [1e-5,10]

    ;---constant.
        gamma = '!9'+string(71b)+'!X'
        rad = !dpi/180d
        deg = 180d/!dpi
        
        omass = 16*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
        hmass = 1*(1.67e-27/1.6e-19)    ; E in eV, mass in kg.


    ;---patch scinfo add r_oh1, r_oh2.
    if patch1 then begin
    ;---load pre-calced data.
        store_data, '*', /delete
        infofn = shomedir()+'/Google Drive/works/data/cusp/'+$
            eventid+'_scinfo.tplot'
        tplot_restore, filename = infofn

        get_data, 'scinfo', tutr, scinfo
        


    ;---polar.
        tinfo = scinfo.polar
        plotutr = tinfo.plot_time
        cusputr = tinfo.cusp_time
        pre0 = 'po_'

        dat = sread_polar_timas(plotutr, type='h0')
        if size(dat,/type) ne 8 then begin
            str_element, tinfo, 'r_oh_cusp', 0d, /add_replace
            str_element, tinfo, 'r_oh_pcap', 0d, /add_replace
            str_element, scinfo, 'polar', tinfo, /add_replace
        endif else begin
            uts = sfmepoch(dat.epoch_h, 'unix')
            utidx = where(uts ge plotutr[0] and uts le plotutr[1], nrec)
            if nrec eq 0 then begin
                str_element, tinfo, 'r_oh_cusp', 0d, /add_replace
                str_element, tinfo, 'r_oh_pcap', 0d, /add_replace
                str_element, scinfo, 'polar', tinfo, /add_replace
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
                store_data, pre1+'n', uts, ns, $
                    limits = {ytitle:'(cm!U-3!N)', ylog:1, yrange:nyrng}
                    
                    
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
                store_data, pre1+'n', uts, ns, $
                    limits = {ytitle:'(cm!U-3!N)', ylog:1, yrange:nyrng}
                    
                    
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
                
                
                str_element, tinfo, 'r_oh_cusp', double(no1/nh1), /add_replace
                str_element, tinfo, 'r_oh_pcap', double(no2/nh2), /add_replace
                str_element, scinfo, 'polar', tinfo, /add_replace
            endelse
        endelse
        
        
        
    ;---fast.
        tinfo = scinfo.fast
        plotutr = tinfo.plot_time
        cusputr = tinfo.cusp_time
        pre0 = 'fa_'
        dat = sread_fast_tms(plotutr)
        if size(dat,/type) ne 8 then begin
            str_element, tinfo, 'r_oh_cusp', 0d, /add_replace
            str_element, tinfo, 'r_oh_pcap', 0d, /add_replace
            str_element, scinfo, 'fast', tinfo, /add_replace
        endif else begin    ; have data in that day.
            uts = sfmepoch(dat.epoch, 'unix')
            utidx = where(uts ge plotutr[0] and uts le plotutr[1], nrec)
            if nrec eq 0 then begin
                str_element, tinfo, 'r_oh_cusp', 0d, /add_replace
                str_element, tinfo, 'r_oh_pcap', 0d, /add_replace
                str_element, scinfo, 'fast', tinfo, /add_replace
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
                store_data, pre1+'n', uts, ns, $
                    limits = {ytitle:'(cm!U-3!N)', ylog:1, yrange:nyrng}
                    
                    
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
                store_data, pre1+'n', uts, ns, $
                    limits = {ytitle:'(cm!U-3!N)', ylog:1, yrange:nyrng}
                    
                    
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
                
                
                str_element, tinfo, 'r_oh_cusp', double(no1/nh1), /add_replace
                str_element, tinfo, 'r_oh_pcap', double(no2/nh2), /add_replace
                str_element, scinfo, 'fast', tinfo, /add_replace
            endelse
        endelse
        
        
        ;---save scinfo.
        store_data, 'scinfo', tutr, scinfo
        tplot_save, 'scinfo', filename=infofn
    endif




end

;id = '1998_0925_05'
;id = '1998_1001_02'
;cusp_patch_saved_data, id

;ids = cusp_id('all')
;foreach id, ids do begin
;    print, '----processing '+id+' ...'
;    cusp_patch_saved_data, id
;endforeach
end
