;+
; Type: function.
; Purpose: Return info structure for cusp conjunction study.
; Parameters: none.
; Keywords: none.
; Return: struct.
; Notes: none.
; Dependence: none.
; History:
;   2015-01-22, Sheng Tian, create.
;-

function cusp_info_struct

    ; max # of filters.
    maxnfilter = 16

    ; unix time, ilat in deg, mlt in hr, dis in re.
    info_cusp = {ut:0d, ilat:0d, mlt:0d, dis:0d}

    ; integrated flux of low, mid, high freq bands, and all bands.
    ; fh: < faclim. Due to Alfven wave.
    ; fl: > faclim. Due to steady state and FAC.
    ; fs: total poynting flux at each band.
    info_pflux = {fh:0d, fl:0d, fs:dblarr(maxnfilter)}

    ; needed info for Polar and FAST.
    ; all energy fluxes are mapped and integrated, + is downward, - is upward.
    info_sc = {$
        cusp:{entry:info_cusp, exit:info_cusp},$; cusp entry/exit info.
        va:0d, $                                ; Va in km/s, rough estimation.
        hem:0d, $                               ; 1:north,-1:south hemisphere.
        kei:0d, $                               ; KEi.
        kee:0d, $                               ; KEe.
        nfilter:0d, $                           ; # of used filters.
        filters:dblarr(maxnfilter), $           ; filters.
        facdelim:0d, $                          ; limit for fac.
        noisedelim:dblarr(2), $                 ; limit [high,low]freq noise.
        ebratio:dblarr(maxnfilter,3), $         ; [E/B ratio, Emax, Bmax].
        sbratio:0, $                            ; sb/abs(sb+sv+sp).
        pfstar:dblarr(maxnfilter), $            ; poynting flux correct coef.
        sb:info_pflux, sv:info_pflux, sp:info_pflux, $ ; poynting flux [b,v,p].
        pflux:0d, $                             ; total poynting flux.
        eflux:0d}                               ; total energy flux.

    info = {$
        id:'', $            ; event id.
        conjtr:dblarr(2), $ ; conjunction time range.
        ae:0d, $            ; max AE, and min Dst in nT within conjunction,
        dst:0d, $           ; min[Polar,FAST] entry, max[Polar,FAST] exit.
        imfbz:dblarr(2), $  ; [min,max]Bz during conjunction.
        imfby:dblarr(2), $  ; [min,max]By during conjunction.
        dt:0d, $            ; in hour, FAST-Polar cusp entry time. 
        dmlt:0d, $          ; in hour, FAST-Polar middle MLT.
        dr:0d, $            ; in Re, Polar-FAST middle distance.
        max_kei:0d, $       ; polar max upward KEi mapped abs value.
        ratio_ion:0d, $     ; (up-down)/(up+down) KEi, -1: all down, 1: all up.
        ratio_eflux:0d, $   ; polar/fast eflux.
        ratio_pflux:0d, $   ; polar/fast pflux.
        polar:info_sc, $
        fast:info_sc}

    return, info
    
end    
