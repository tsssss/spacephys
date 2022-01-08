;+
; Plot |B| of specified spacecraft and their position.
;-
;

    test = 0

;---Constant.
    deg = 180d/!dpi
    rad = !dpi/180d

;---Settings.
    utr0 = time_double(['2014-08-28/00:00','2014-08-29/00:00'])
    utr1 = time_double(['2014-08-28/04:30','2014-08-28/06:00'])     ; the time for plot.
    utr2 = time_double(['2014-08-28/04:55','2014-08-28/05:05'])     ; the time to count MLon.
    reload = 0

    all_probes = ['rbspa','rbspb','tha','thd','the','g13','g15']
    themis_probes = ['a','d']
    rbsp_probes = ['a']
    goes_probes = ['13','15']
    probes = ['rbspa','thd','tha','g13','g15']


;---Read R GSM.
    load = 0
    foreach tvar, probes+'_r_gsm' do if tnames(tvar) eq '' then begin
        load = 1
        break
    endif
    if keyword_set(reload) then load = 1
    if load then begin
        foreach probe, rbsp_probes do rbsp_read_orbit, utr0, probe=probe
        foreach probe, themis_probes do themis_read_orbit, utr0, probe=probe
        foreach probe, goes_probes do goes_read_orbit, utr0, probe=probe
    endif


;---Read B GSM.
    load = 0
    foreach tvar, probes+'_b_gsm' do if tnames(tvar) eq '' then begin
        load = 1
        break
    endif
    if keyword_set(reload) then load = 1
    if load then begin
        foreach probe, rbsp_probes do rbsp_read_bfield, utr0, probe=probe
        foreach probe, themis_probes do themis_read_bfield, utr0, probe=probe
        foreach probe, goes_probes do goes_read_bfield, utr0, probe=probe

        foreach tvar, probes+'_r_gsm' do read_geopack_info, tvar
        foreach probe, probes do begin
            bvar = probe+'_b_gsm'
            if tnames(bvar) eq '' then continue
            get_data, bvar, uts, bgsm
            bmag = snorm(bgsm)

            bmodvar = probe+'_bmod_gsm'
            get_data, bmodvar, tmp, bmodgsm
            bmodmag = snorm(bmodgsm)
            dbmag = bmag-interpol(bmodmag,tmp,uts)
            sdespike, uts, dbmag
            store_data, probe+'_dbmag', uts, dbmag
            add_setting, probe+'_dbmag', /smart, {$
                display_type: 'scalar', $
                unit: 'nT', $
                short_name: '|B|-|B!Dmod!N|'}
        endforeach
    endif




;---Calc MLon and MLT.
    load = 0
    foreach tvar, probes+'_mlon' do if tnames(tvar) eq '' then begin
        load = 1
        break
    endif
    if keyword_set(reload) then load = 1
    if load then begin
        foreach probe, probes do begin
            var = probe+'_r_gsm'
            get_data, var, times, rgsm
            ntime = n_elements(times)
            mlons = fltarr(ntime)
            mlts = fltarr(ntime)
            for ii=0, ntime-1 do begin
                ps = geopack_recalc(times[ii])
                geopack_conv_coord, rgsm[ii,0],rgsm[ii,1],rgsm[ii,2], /from_gsm, tx,ty,tz, /to_mag
                mlons[ii] = atan(ty,tx)*deg
                mlts[ii] = mlon2mlt(mlons[ii], times[ii])
            endfor
            store_data, probe+'_mlon', times, mlons, limits={$
                ytitle:'(deg)', labels:'MLon'}
            store_data, probe+'_mlt', times, mlts, limits={$
                ytitle:'(hr)', labels:'MLT'}
        endforeach
    endif



;---Plot |B|.
    pres = ['tha','g13','g15']+'_'
    npre = n_elements(pres)
    ; Sort by MLon.
    mlons = fltarr(npre)
    for ii=0, npre-1 do begin
        get_data, pres[ii]+'mlon', times, data
        mlons[ii] = median(data[where(times ge utr2[0] and times le utr2[1])])
        print, pres[ii]+'MLon (deg):', mlons[ii]
    endfor
    index = sort(mlons)
    vars = pres[index]+'dbmag'
    nvar = n_elements(vars)
    figlabels = ['a','b','c']+'.'

    tvar = 'g15_dbmag'
    options, tvar, 'labels', 'GOES-15'
    options, tvar, 'yrange', [-45,-5]
    options, tvar, 'ytickv', [-40,-25,-10]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    tvar = 'g13_dbmag'
    options, tvar, 'labels', 'GOES-13'
    options, tvar, 'yrange', [-5,50]
    options, tvar, 'ytickv', [0,25,50]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    tvar = 'tha_dbmag'
    options, tvar, 'labels', 'TH-A'
    options, tvar, 'yrange', [-15,45]
    options, tvar, 'ytickv', [0,20,40]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    ofn = sparentdir(srootdir())+'/plot/fig_case1_bmag_and_pos.pdf'
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=4, ysize=3

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    poss = sgcalcpos(nvar, lmargin=6, tmargin=1, bmargin=4, rmargin=8)
    tplot, vars, trange=utr1, position=poss

    tpos = poss[*,0]
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, '|B|-|B!Dmod!N|'

    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        xyouts, xchsz*1, tpos[3]-ychsz*0.8, /normal, figlabels[ii]
    endfor

    sgclose
end
