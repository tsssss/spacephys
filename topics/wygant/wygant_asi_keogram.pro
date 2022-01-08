
pro wygant_asi_keogram, tr

    ; settings.
    range_mlat = [60,68]
    range_mlon = [-65,-45]
    ; for keogram.
    dmlon = 0.1 ; deg.
    dmlat = 0.1 ; deg.
    ; for photon count.
    delmlon = .1 ; deg.
    delmlat = .1 ; deg.
    h0 = 110 ; km. altitude.
    ; load asf images.
    sites = ['atha']
    vars = ['mlat','mlon','alti']
    asf = sread_thg_asi(tr, sites, type = 'asf')
    asc = sread_thg_asc(tr, sites, vars = vars, type = 'asf')
    asi = asf.img
    t0 = asf.utsec
    nrec = n_elements(t0)
    idx = where(asc.(0).alti eq h0*1e3)
    asi_mlons = reform(asc.(0).mlon[idx,*,*])
    asi_mlats = reform(asc.(0).mlat[idx,*,*])
    ; convert corner position to center position.
    npx = (size(asi_mlons, /dimensions))[0]
    asi_mlons = 0.5*(asi_mlons[1:npx-1,*]+asi_mlons[0:npx-2,*])
    asi_mlons = 0.5*(asi_mlons[*,1:npx-1]+asi_mlons[*,0:npx-2])
    asi_mlats = 0.5*(asi_mlats[1:npx-1,*]+asi_mlats[0:npx-2,*])
    asi_mlats = 0.5*(asi_mlats[*,1:npx-1]+asi_mlats[*,0:npx-2])
    asi_mlons = asi_mlons[*]
    asi_mlats = asi_mlats[*]
    
    ; load rbsp pos, fpoint, etc.
    fn = shomedir()+'/Google Drive/works/data/rbsp_thg/rbsp_efw_fld_2013_0501_04.tplot'
    tplot_restore, filename = fn
    
    ; rbspb_fpt_lonlat contains footpoint mlon, mlat.
    get_data, 'rbspb_fpt_lonlat', tmp, dat
    rbspb_mlons = interpol(dat[*,0], tmp, t0)
    rbspb_mlats = interpol(dat[*,1], tmp, t0)
    
    ; extract points around rbsp mlon, and between mlat range.
    lats = smkarthm(range_mlat[0], range_mlat[1], dmlat, 'dx')
    lons = smkarthm(range_mlon[0], range_mlon[1], dmlon, 'dx')
    nlat = n_elements(lats)
    nlon = n_elements(lons)
    keo = dblarr(nrec,nlat)
    ewo = dblarr(nrec,nlon)
    crs = dblarr(nrec)
    for i = 0, nrec-1 do begin
        timg = reform(asf.img[i,*,*]) & timg = timg[*]
        
        ; extract keogram.
        tmlon = rbspb_mlons[i]
        idx = where(abs(asi_mlons-tmlon) le dmlon, cnt)
        tlats = asi_mlats[idx]
        keo[i,*] = interpol(timg[idx],tlats,lats)
        ; remove extrapolation.
        idx = where(lats lt min(tlats) or lats gt max(tlats), cnt)
        if cnt ne 0 then keo[i,idx] = !values.d_nan
        
        ; extract ewogram.
        tmlat = rbspb_mlats[i]
        idx = where(abs(asi_mlats-tmlat) le dmlat, cnt)
        tlons = asi_mlons[idx]
        tmp = sort(tlons)
        ewo[i,*] = interpol((timg[idx])[tmp],tlons[tmp],lons)
        ; remove extrapolation.
        idx = where(lons lt min(tlons) or lons gt max(tlons), cnt)
        if cnt ne 0 then ewo[i,idx] = !values.d_nan
        
        ; count rate.
        idx = where(abs(asi_mlons-tmlon) le delmlon and $
            abs(asi_mlats-tmlat) le delmlat, cnt)
        crs[i] = cnt? max(timg[idx]): !values.d_nan
    endfor

    ; save keogram and ewogram.
    var = 'thg_atha_keo'
    store_data, var, t0, keo, lats
    options, var, 'zrange', minmax(keo)
    var = 'thg_atha_ewo'
    store_data, var, t0, ewo, lons
    options, var, 'zrange', minmax(ewo)
    vars = 'thg_atha_'+['keo','ewo']
    options, vars, 'spec', 1
    options, vars, 'ztitle', 'Photon Count'
    options, vars, 'zlog', 0
    options, vars, 'no_interp', 1
    
    ; save rbspb footpoint.
    store_data, 'rbspb_mlat', t0, rbspb_mlats, limits = {colors:6, ytitle:'MLat'}
    store_data, 'rbspb_mlon', t0, rbspb_mlons, limits = {colors:6, ytitle:'MLon'}

    ; save photon count.
    var = 'rbspb_atha_count'
    idx = where(~finite(crs,/nan))
    crs = interpol(crs[idx], t0[idx], t0)
    store_data, var, t0, crs
    options, var, 'labels', snum2str(delmlon,/shortfloat)+'x'+snum2str(delmlat,/shortfloat)+' deg'
;    options, var, 'psym', 0
    
    var = 'comb_keo'
    store_data, var, data = ['thg_atha_keo','rbspb_mlat']
    options, var, 'yrange', range_mlat
    options, var, 'ytitle', 'mlat (deg)'
    
    var = 'comb_ewo'
    store_data, var, data = ['thg_atha_ewo','rbspb_mlon']
    options, var, 'yrange', range_mlon
    options, var, 'ytitle', 'mlon (deg)'
        
    vars = ['comb_keo','comb_ewo','rbspb_pf_fac_mat','rbspb_atha_count']
    labs = ['rbspb_mlat','rbspb_mlon']
    titl = 'RBSP-B footpoint over Themis/ASI keogram and ewogram (site ATHA)'
    tmp = time_string(tr[0],tformat='hhmm')+'-'+time_string(tr[1],tformat='hhmm')
    ofn = shomedir()+'/tmp/2013_0501_keogram_'+tmp+'.pdf'
    sgopen, ofn
    sgindexcolor, 43
    tplot, vars, trange = tr, var_label = labs, title = titl
    sgclose
end

;tr = time_double(['2013-05-01/07:30','2013-05-01/07:50'])
;wygant_asi_keogram, tr
;tr = time_double(['2013-05-01/05:00','2013-05-01/10:00'])
;wygant_asi_keogram, tr
;tr = time_double(['2013-05-01/07:15','2013-05-01/08:15'])
;wygant_asi_keogram, tr
;tr = time_double(['2013-05-01/05:30','2013-05-01/06:10'])
;wygant_asi_keogram, tr

tr = time_double(['2013-05-01/07:25','2013-05-01/07:55'])
wygant_asi_keogram, tr


end