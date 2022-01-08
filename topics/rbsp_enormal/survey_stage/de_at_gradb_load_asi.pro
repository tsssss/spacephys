
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
minlat = 50                 ; min lat for aurora image.
tpos = [0.1,0.1,0.9,0.9]    ; position for aurora image.
asixr = [-90,90]            ; longitude/local time range for aurora image.
asiyr = [50,90]             ; latitude range for aurora image.
rgb = [6,4,2]


logfn = shomedir()+'/Google Drive/works/works/rbsp_de/de_at_gradb/list_de_at_gradb_asi_round2.log'
nheader = 3
headers = strarr(nheader)
nline = file_lines(logfn)-nheader
lines = strarr(nline)
openr, lun, logfn, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun

asc = sread_thg_asc(0, '*', vars = ['mlon','midn'])
nsite = n_elements(tag_names(asc))
cmidn = fltarr(nsite)
cmlon = fltarr(nsite)
for i = 0, nsite-1 do begin
    cmlon[i] = asc.(i).mlon
    tmp = asc.(i).midn
    if tmp eq '     ' then cmidn[i] = !values.d_nan $
    else cmidn[i] = float(strmid(tmp,0,2))+float(strmid(tmp,3,2))/60.0
endfor
idx = sort(cmidn) & cmlon = cmlon[idx] & cmidn = cmidn[idx]


for i = 0, nline-1 do begin
    tfn = file_basename(lines[i])
    tprobe = strlowcase(strmid(tfn,19,1))
    id = strmid(tfn,0,14)
    ut0 = time_double(id, tformat='YYYY_MMDD_hhmm')
    utr = ut0+[-2,4]*60
    sites = strtrim(strmid(tfn,59),2)
    pre0 = 'rbsp'+tprobe
    
    ; check site.
    if sites[0] eq 'na' then continue
    sites = strsplit(sites,', ',/extract)
    print, sites
    
    
    ; load s/c position.
    efwl3 = sread_rbsp_efw_l3(utr, probes = tprobe)
    if size(efwl3,/type) ne 8 then continue
    uts = sfmepoch(efwl3.epoch,'unix',/epoch16)
    
    store_data, pre0+'mlt', uts, (efwl3.mlt_lshell_mlat)[*,0], $
        limits = {ytitle:'MLT (hr)'}
    store_data, pre0+'lshell', uts, (efwl3.mlt_lshell_mlat)[*,1], $
        limits = {ytitle:'L (Re)'}
    store_data, pre0+'ilat', uts, acos(1/sqrt((efwl3.mlt_lshell_mlat)[*,1]))*deg, $
        limits = {ytitle:'ILat (deg)'}
    store_data, pre0+'mlat', uts, (efwl3.mlt_lshell_mlat)[*,2], $
        limits = {ytitle:'MLat (deg)'}
    store_data, pre0+'pos_gse', uts, efwl3.pos_gse, $
        limits = {colors:rgb, labels:['x','y','z']}


    ; prepare mapping.
    model = 't89'
    dir = -1        ; always north hem, b/c conjugate to thm_asi.
    r0 = 1+110*re1  ; 110 km altitude.
    sgeopack_par, utr, model, /delete  ; get tplot var <model>_par.
    t89 = 0 & t96 = 0 & t01 = 0
    case model of
        't89': t89 = 1
        't96': t96 = 1
        't01': t01 = 1
    endcase
    ; map pos to fpt, convert to mag.
    get_data, pre0+'pos_gse', data = tmp
    uts = tmp.x & ets = 1000D*uts+62167219200000D
    pos0 = tmp.y*re1 & pos1 = pos0      ; in re.
    ; interpolate par.
    if model ne 't89' then begin
        get_data, model+'_par', data = tmp
        pars = sinterpol(tmp.y, tmp.x, uts)
    endif
    ; loop for each time.
    for j = 0, n_elements(uts)-1 do begin
        ; set geopack.
        geopack_epoch, ets[j], yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001D, /date
        ; pos in gse, which is the mapping coord.
        x0 = pos0[j,0] & y0 = pos0[j,1] & z0 = pos0[j,2]
        ; convert from gse to gsm.
        geopack_conv_coord, x0, y0, z0, /from_gse, $
            x1, y1, z1, /to_gsm
        dir = (z1 gt 0)? -1: 1
        if model ne 't89' then par = reform(pars[j,*]) else par = 2
        geopack_trace, x1, y1, z1, dir, par, xf, yf, zf, $
            epoch = ets[j], /refine, /ionosphere, $
            t89 = t89, t96 = t96, t01 = t01
        ; convert from gse to mag.
        geopack_conv_coord, xf, yf, zf, /from_gsm, $
            x1, y1, z1, /to_mag
        pos1[j,*] = [x1,y1,z1]
    endfor
    mlat = asin(pos1[*,2]*(1/r0))*deg
    mlon = atan(pos1[*,1],pos1[*,0])*deg
    store_data, pre0+'fpt_mag', data = {x:uts, y:pos1}
    store_data, pre0+'fpt_mlat', data = {x:uts, y:mlat}, $
        limits = {ytitle:'MLat/fpt (deg)'}
    store_data, pre0+'fpt_mlon', data = {x:uts, y:mlon}, $
        limits = {ytitle:'MLon/fpt (deg)'}
    
    
    ; check hemisphere.
    get_data, pre0+'fpt_mlat', data = tmp
    if tmp.y[0] lt 0 then continue
    
    
    ; load thg asi.
    asi = sread_thg_mlt(utr, sites, type = 'asf', /half)
    nrec = n_elements(asi.epoch)
    szs = size(reform(asi.mltimg[0,*,*]),/dimensions)
    
    
    ; output.
    get_data, pre0+'fpt_mlon', tmp, mlon
    mlon = sinterpol(mlon, tmp, sfmepoch(asi.epoch,'unix'))
    get_data, pre0+'fpt_mlat', tmp, mlat
    mlat = sinterpol(mlat, tmp, sfmepoch(asi.epoch,'unix'))
    
    rootdir = shomedir()+'/de_at_gradb/asi/'+id
    for j = 0, nrec-1 do begin
        print, 'processing '+sfmepoch(asi.epoch[j])+' ...'
        ofn = rootdir+'/thg_asf_'+sfmepoch(asi.epoch[j],'YYYY_MMDD_hhmm_ss')+'.png'
        sgopen, ofn, xsize = szs[0], ysize = szs[1]
        sgtruecolor
        sgtv, reform(asi.mltimg[j,*,*]), position = tpos, ct = 1
        ; set up coord.
        sgset_map, xrange = asixr, yrange = asiyr, position = tpos, $
            ytickv = smkarthm(asiyr[0],asiyr[1],10,'dx'), xtickv = smkarthm(asixr[0],asixr[1],3,'n'),$
            ytickpos = 30, xtickpos = minlat+1.2, xtickname = ['21','00','03'], xminor = 2, $
            xticknudge = [[0,-0.5],[0,0],[0,-0.5]], yticknudge = [0.75,0], $
            color = sgcolor('white'), charsize = 1, majorgrid = 1, minorgrid = 1
        xyouts, tpos[0], tpos[1], /normal, 'ASF: '+sfmepoch(asi.epoch[j],'YYYY-MM-DD/hh:mm:ss'), color = sgcolor('white')
        
        midn = interpol(cmlon, cmidn, (sfmepoch(asi.epoch[j],'unix')/86400d mod 1)*24, /nan)
        plots, mlon[j]-midn, mlat[j], psym = 1, color = sgcolor('red')        
        sgclose
    endfor
endfor

end
