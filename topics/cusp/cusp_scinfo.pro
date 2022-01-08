function cusp_scinfo

    tinfo = {$
        plot_time: [0d,0], $    ; time range for plot.
        cusp_time: [0d,0], $    ; time range for cusp crossing.
        cusp_dis: [0d,0], $     ; distance range, in Re.
        cusp_mlt: [0d,0], $     ; MLT range, in hr.
        cusp_ilat: [0d,0], $    ; ilat range, in deg.
        cusp_fpt_mlt: [0d,0], $ ; footpoint MLT range, in deg.
        cusp_fpt_mlat: [0d,0], $; footpoint MLat range, in deg.
        b0: 0d, $               ; median background B, in nT.
        dis0: 0d, $             ; distance, in Re.
        mlt0: 0d, $             ; MLT, in hr.
        ilat0: 0d, $            ; ILat, in deg. Positive for northern hemisphere.
        dr0: 0d, $              ; field data rate, in sec.
        dr1: 0d, $              ; particle data rate, in sec.
        n0: 0d, $               ; max density, in cc.
        va: 0d, $               ; local Va in km/s, assume all H+.
        vsc: 0d, $              ; average Vsc in km/s.
        vdir: 0d, $             ; 1: low to high latitude, 0: otherwise.
        scaleinfo:$             ; info for Morlet wavelet.
            {s0:0d,s1:0d,dj:1d/8,ns:0d}, $
        filters: [0d,0], $      ; filters in sec.
        int_kei: 0d, $          ; integrated KEi, in W/m. Positive for parallel.
        int_kee: 0d, $          ; integrated KEe, in W/m. Positive for parallel.
        int_pfb: 0d, $          ; integrated S para, in W/m. Positive for parallel.
        int_pfs: [0d,0d,0d], $  ; integrated S FAC, in W/m.
        r_kei: 0d, $            ; ion ratio. In terms of parallel, not upward.
        r_kee: 0d, $            ; electron ratio. Same as above.
        r_pfb: 0d, $            ; S para ratio.
        r_ion: 0d, $            ; ion ratio. >0 for earthward, <0 for upward.
        max_kei: 0d, $          ; maximum absolute value of the mapped ion energy flux.
        r_oh_cusp: 0d, $        ; ratio of n(O/H) within cusp.
        r_oh_pcap: 0d}          ; ratio of n(O/H) within polar cap.
    scinfo = {$
        id: '', $
        model: '', $
        ae: 0d, $
        dst: 0d, $
        fast_orbit: 0d, $
        deidx: 0d, $            ; v: north-south.
        dbidx: 1d, $            ; p: east-west.
        pfidx: 2d, $            ; b: parallel.
        polar: tinfo, $
        fast: tinfo}

    return, scinfo
end
