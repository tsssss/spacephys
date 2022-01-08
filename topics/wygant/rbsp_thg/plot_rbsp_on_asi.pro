;+
; asivars: vars plotted on aurora images.
; pltvars: vars for tplot.
; vecvars: vectors plotted on aurora image.
; 
; ; run read_thm_asi_mosaic and rbsp_thg_prep_data first.
;-
pro plot_rbsp_on_asi, probe = probe, sites = sites, trange = rbtr, time = t0, $
    rootdir = rootdir, efwfn = efwfn, asifn = asifn, sg = sg, minlat = minlat


    ; preparing.
    if n_elements(probe) eq 0 then probe = 'a'
    rbx = 'rbsp'+strlowcase(probe[0])

    if n_elements(sites) eq 0 then message, 'asi sites undefined ...'
    if n_elements(minlat) eq 0 then minlat = 50

    ; time range to draw line plots.
    if n_elements(rbtr) ne 2 then message, 'time range undefined ...'
    rbetr = stoepoch(rbtr)
    utr = time_double(rbtr)
    uts = smkarthm(utr[0],utr[1],3,'dx')

    if n_elements(rootdir) eq 0 then $
        rootdir = shomedir()+'/Google Drive/works/data/rbsp_thg'

    if n_elements(sg) eq 0 then sg = 'w'
    ext = (sg eq 'ps')? 'pdf': 'png'


    ; rbsp e/b field, pos, fpoint, pflux.
    if n_elements(efwfn) eq 0 then begin
        efwfn = 'rbsp_efw_fld_'+sfmepoch(rbetr[0],'YYYY_MMDD_hh')+'.tplot'
        efwfn = rootdir+'/'+efwfn
    endif

    ; themis aurora.
    if n_elements(asifn) eq 0 then begin
        asifn = 'thg_asf_mosaic_'+sfmepoch(rbetr[0],'YYYY_MMDD_hhmm')+'_'+sfmepoch(rbetr[1],'hhmm')+'.cdf'
        asifn = rootdir+'/'+asifn
    endif

    if file_test(efwfn) eq 0 then rbsp_thg_prep_data, rbtr, probes = probe, sites = sites, mosaic = 0, ofn = efwfn
    if file_test(asifn) eq 0 then asifn = sread_thg_mosaic(rbtr, sites, ofn = asifn)


    ; load rbsp efw data.
    ; rbsp[a,b]_pos_gsm, [e,b]_gse, bmod_gse, [de,db]_fac,
    ; fpt_lonlat, pf_fac_mat, fpt_lonlat_thg_[sites]_cr.
    tplot_restore, filename = efwfn

    ; load themis asi data.
    asiid = cdf_open(asifn)
    vars = ['time','image_size','pixel_index','minlat','midn']
    asidat = scdfread(asiid, vars, skt=skt)
    asiuts = *(asidat[0].value)
    imgszs = *(asidat[1].value)
    pxidxs = *(asidat[2].value)
    minlat = *(asidat[3].value) & minlat = minlat[0]
    midns = *(asidat[4].value)


; **** plot settings.
    ; plot size, in inch.
    xsz = 7d
    ysz = 10d
    
    ; aurora image size and position.
    asixsz = 6d
    asiysz = 0.5*asixsz
    asipos = [(1-asixsz/xsz)*0.5,0,(1+asixsz/xsz)*0.5,0.95]
    asipos[1] = asipos[3]-double(asixsz)/xsz*imgszs[1]/imgszs[0]*xsz/ysz

    ; aurora range.
    xr = [-45,45]
    yr = [minlat,70]

    ; aurora index.
    isz = max(imgszs)
    xc = smkarthm(-1d,1,isz,'n') # (dblarr(isz)+1)
    yc = smkarthm(-1d,1,isz,'n') ##(dblarr(isz)+1)
    rc = sqrt(xc^2+yc^2)
    tc = atan(yc, xc)
    lat = 90-rc*(90-minlat)             ; in deg.
    lon = (tc/!dpi*180+360) mod 360     ; in deg.
    lat = (rotate(lat,1))[*,0:isz/2]
    lon = (rotate(lon,1)-180)[*,0:isz/2]
    
    idx = where(lon ge xr[0] and lon le xr[1] and lat ge yr[0] and lat le yr[1])
    idx = array_indices(lon, idx)
    xidx = minmax(idx[0,*])+[-1,1]
    yidx = minmax(idx[1,*])+[-1,1]
    xidx[0] >= 0 & xidx[1] <= isz-1
    yidx[0] >= 0 & yidx[1] <= isz/2
    
    
    !p.charsize = 1d

    tplot_options, 'labflag', 1
    tplot_options, 'num_lab_min', 10

    red = sgcolor('red')
    grn = sgcolor('lime')
    blu = sgcolor('blue')
    wht = sgcolor('white')
    blk = sgcolor('black')

    ; pflux.
    var = rbx+'_pf_fac_mat'
    options, var, 'labels', ['x','y','z']
    options, var, 'ytitle', 'S!C(mW/m!U2!N)'
    options, var, 'colors', [red,grn,blu]
    ;get_data, var, tmp, dat
    ;store_data, var+'_para', tmp, dat[*,0], limits = {ytitle:'S para!C(mW/m!U2!N)'}
    ;ylim, var+'_para', 1e-4, 1, 1

    ; e/b fields.
    vars = rbx+['_de_fac','_db_fac']
    options, vars, 'colors', [red,grn,blu]

    ; count rate along s/c footpoint.
    vars = rbx+['_fpt_lonlat_thg_'+sites+'_cr']
    ylim, vars, 1, 100, 1

    ; vector plotted on the aurora image.
    scle = 0.05     ; 1 mV/m ~ 10 deg.
    scle = 0.1      ; for log field.
    get_data, 'rbspb_de_fac', tmp, dat
    store_data, 'rbspb_de_vec', tmp, dat[*,0:1], scle, limits = {colors:red, linestyle:0}
    store_data, 'rbspb_fpt_lonlat', limits = {colors:wht, linestyle:1}


    ; vars plotted on the aurora image using plots.
    asivars = 'rbspb_fpt_lonlat'    ; s/c footprint.
    ; vars plotted on the aurora image as special vector.
    vecvars = [['rbspb_fpt_lonlat'],['rbspb_de_vec']]   ; [origin,length] pair.
    ; vars plotted using tplot.
    pltvars = ['rbspb_'+['de_fac','db_fac','vexb','pf_fac_mat','fpt_lonlat_thg_'+sites+'_cr']]
    
    ; log scale for vector vars.
    log = 1

    if sg eq 'ps' then ext = 'pdf' else ext = 'png'
    ofnptn = shomedir()+'/'+rbx+'_on_asi_YYYY_MMDD_hhmm_ss.'+ext


    ; loop through the times.
    ; time to draw time bar.
    if n_elements(t0) ne 0 then utbars = time_double(t0) else utbars = uts
    for i = 0, n_elements(utbars)-1 do begin
        ; reconstruct the aurora image at the time.
        ; tut, timg, tmidn.
        tut = utbars[i]
        tet = stoepoch(tut,'unix')
        print, 'processing '+time_string(tut)+' ...'
        idx = where(tut eq asiuts, cnt)
        if cnt eq 0 then continue
        timg = bytarr(imgszs)
        timg[pxidxs] = *((scdfread(asiid,'thg_mosaic',idx,skt=skt))[0].value)
        tmidn = midns[idx]

        ; create a canvas.
        ofn = (sg eq 'w')? 0: sfmepoch(tet,ofnptn)
        sgopen, ofn, xsize = xsz, ysize = ysz, /inch
        sgtruecolor
        loadct, 1
        ; trim the image to the lat/lon range then plot.
        timg = timg[xidx[0]:xidx[1],yidx[0]:yidx[1]]
        sgtv, timg, position = asipos, /nan
        sgset_map, xrange = xr, yrange = yr, pos = asipos, $
            ytickv = smkarthm(yr[0],yr[1],5,'dx'), ytickpos = 30, $
            xtickv = smkarthm(xr[0],xr[1],3,'n'), xtickpos = minlat+1.2, $
            xtickname = ['21','00','03'], xminor = 2, $
            xticknudge = [[0,-0.5],[0,-0.5],[0,-0.5]], yticknudge = [0.75,0], $
            color = wht, charsize = 1, majorgrid = 1, minorgrid = 1
        xyouts, asipos[0], asipos[1], color = blk, $
            'ASI: '+time_string(tut,tformat='hh:mm:ss'), /normal

        ; plot directly using plots.
        !p.color = blk
        !p.background = wht
        vars = asivars
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin  ; need to account for midnight longitude.
            get_data, vars[i], t0, dat, limits = lim
            plots, dat[*,0]-tmidn[0], dat[*,1], color = lim.colors, linestyle = lim.linestyle
        endfor

        ; plot vectors using plots after conversion.
        vars = vecvars
        nvar = n_elements(vars)/2
        for i = 0, nvar-1 do begin
            get_data, vars[i,0], t0, dat
            p0 = sinterpol(dat,t0, tut) & p0[0] -= tmidn[0]
            get_data, vars[i,1], t0, dat, scl, limits = lim
            v0 = sinterpol(dat,t0, tut)
            p1 = svec2map(v0, p0, scl, log = log)
            plots, [p0[0],p1[0]], [p0[1],p1[1]], color = lim.colors, linestyle = lim.linestyle
        endfor

        ; plot using tplot.
        vars = pltvars
        nvar = n_elements(vars)
        varpos = sgcalcpos(nvar, region = [0d,0,1,asipos[1]])
        tplot, pos = varpos, vars, /noerase, trange = tr
        timebar, tut, color = red

        sgclose

    endfor
    
    cdf_close, asiid
end



rootdir = shomedir()+'/Google Drive/works/data/rbsp_thg'
sites = ['atha','tpas']     ; the sites to load photon count and aurora image.
probe = 'b'
rbtr = ['2013-05-01/04:00','2013-05-01/10:00']  ; the time range to load rbsp data and photon count, and tplot these line plots.
ut0 = ['2013-05-01/07:38']                      ; the time for naming outfile and draw time bar.
sg = 'ps'

plot_rbsp_on_asi, probe = probe, sites = sites, trange = rbtr, time = ut0, $
    rootdir = rootdir, efwfn = efwfn, asifn = asifn, sg = sg

end
