;+
; Test different methods to get MLT, either from MLon, or from GLon.
;
; Two methods are used here and shown consistent results:
;   1. A mapping methods used by Polar and Image spacecraft when dealing with auroral images
;   2. A calculation based on solar direction, encapsulated in slon2lt.
;
; The calculated results are checked against to Mlon of the all-sky imagers.
;-
;

;---Pick a random time.
    times = time_double('2014-08-28/08:00')

;---Read the mlon,mlat, glon,glat, midn of all cameras.
;    ; mlon,mlat are in deg. The position of the camera is known for sure.
;    ; glon,glat are in deg. They are also knonw for sure.
;    ; midn in hour. The UT in the day when the camera passes the midn, i.e., MLT=12 hr.
;    sites = themis_read_asi_sites()
;    nsite = n_elements(sites)
;
;    midns = dblarr(nsite)
;    mlons = dblarr(nsite)
;    mlats = dblarr(nsite)
;    glons = dblarr(nsite)
;    glats = dblarr(nsite)
;    foreach site, sites, i do begin
;        asc_var = 'thg_asc_'+site+'_'+['midn','mlon','mlat','glon','glat']
;        themis_read_asi, 0, id='asc', site=site, in_vars=asc_var, skip_index=1
;        str = get_var_data(asc_var[0])
;        if strpos(str,':') eq -1 then midns[i] = 0d else begin
;            tmp = double(strsplit(str,':',/extract))
;            midns[i] = tmp[0]+tmp[1]/60
;        endelse
;        mlons[i] = get_var_data(asc_var[1])
;        mlats[i] = get_var_data(asc_var[2])
;        glons[i] = get_var_data(asc_var[3])
;        glats[i] = get_var_data(asc_var[4])
;        store_data, asc_var, /delete
;    endforeach
;
;    ; For some sites, midn is not set, and midn=0 when that happens.
;    index = where(midns ne 0)
;    midns = midns[index]
;    mlons = mlons[index]
;    mlats = mlats[index]
;    glons = glons[index]
;    glats = glats[index]
;    sites = sites[index]
;    index = sort(mlons)
;    midns = midns[index]
;    mlons = mlons[index]
;    mlats = mlats[index]
;    glons = glons[index]
;    glats = glats[index]
;    sites = sites[index]

    secofday = 86400d
    ut = times[0]
    et = stoepoch(ut, 'unix')
    nsite = n_elements(mlons)

    ; test 1: use the Polar and Image way.
    apexfile = '/Users/Sheng/Projects/idl/slib/sread/support/mlatlon.1997a.xdr'
    geotoapex, glats, glons, apexfile, mlats_polar, mlons_polar
    get_local_time, et, glats, glons, apexfile, glts, mlts_polar
    index = where(mlts_polar gt 12, count)
    if count ne 0 then mlts_polar[index] -= 24

    ; test 2: use my conversion function.
    mlts_sheng = dblarr(nsite)
    for i=0, nsite-1 do mlts_sheng[i] = slon2lt(mlons[i], et, /mag, /degree)

    ; difference from the midnight UT to the current UT, gives the current MLT.
    mlts_asi = (ut mod 86400)/3600-midns


    ofn = srootdir()+'/test_calc_mlt_from_mlon_or_glon.pdf'
    figsz = 8
    sgopen, ofn, xsize=figsz, ysize=figsz, /inch
    poss = sgcalcpos(2,2, position=[0.1,0.1,0.95,0.95], rmargin=5, tmargin=5, xpad=8, ypad=8)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    symsz = 0.5
    ticklen = -0.02

    ; Show MLon agree.
    tpos = poss[*,0,0]
    mlon_range = [-150,50]
    plot, mlons_polar, mlons, /iso, /noerase, $
        xstyle=1, xrange=mlon_range, xticklen=ticklen, xtitle='MLon (deg) calc from the Polar method', $
        ystyle=1, yrange=mlon_range, yticklen=ticklen, ytitle='MLon (deg) read off from ASI', $
        position=tpos, psym=1, symsize=symsz
    plots, mlon_range, mlon_range, linestyle=1
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'a. MLon. ASI vs calc from Polar'

    ; Show MLat agree.
    tpos = poss[*,1,0]
    mlat_range = [55,85]
    plot, mlats_polar, mlats, /iso, /noerase, $
        xstyle=1, xrange=mlat_range, xticklen=ticklen, xtitle='MLat (deg) calc from the Polar method', $
        ystyle=1, yrange=mlat_range, yticklen=ticklen, ytitle='MLat (deg) read off from ASI', $
        position=tpos, psym=1, symsize=symsz
    plots, mlat_range, mlat_range, linestyle=1
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'b. MLat. ASI vs calc from Polar'

    ; Show MLT agree.
    tpos = poss[*,0,1]
    mlt_range = [-12,12]
    plot, mlts_polar, mlts_sheng, /iso, /noerase, $
        xstyle=1, xrange=mlt_range, xticklen=ticklen, xtitle='MLT (hr) calc from the Polar method', $
        ystyle=1, yrange=mlt_range, yticklen=ticklen, ytitle='MLT (hr) calc from slon2lt', $
        position=tpos, psym=1, symsize=symsz
    plots, mlt_range, mlt_range, linestyle=1
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'c. MLT. Sheng vs Polar (both from calc)'

    ; Show MLT disagree.
    tpos = poss[*,1,1]
    mlt_range = [-12,12]
    plot, mlts_polar, mlts_asi, /iso, /noerase, $
        xstyle=1, xrange=mlt_range, xticklen=ticklen, xtitle='MLT (hr) calc from the Polar method', $
        ystyle=1, yrange=mlt_range, yticklen=ticklen, ytitle='MLT (hr) calc from ASI midn', $
        position=tpos, psym=1, symsize=symsz
    plots, mlt_range, mlt_range, linestyle=1
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'd. MLT. ASI vs calc from Polar'

    sgclose
end
