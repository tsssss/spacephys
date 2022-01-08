;+
; asivar: aurora image.
; asivars: vars plotted on aurora images.
; pltvars: vars for tplot.
; vecvars: vectors plotted on aurora image.
;-
pro vars_on_asi, asivar, asivars = asivars, pltvars = pltvars, $
    vecvars = vecvars, trange = tr, sg = sg, minlat = minlat

    xsz = 5d & ysz = 7d     ; inch.
    xsz = 7d & ysz = 10d    ; inch.
    asixsz = 6d & asiysz = 3d   ; inch.
    log = 1
    
    get_data, asivar, ut0, asiimg, midn
    
    if ~keyword_set(sg) then sg = 'w'
    case sg of
        'ps': begin
            ofn = shomedir()+'/rbspb_on_thg_'+$
                time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.eps'
            sgpsopen, ofn, xsize = xsz, ysize = ysz, /inch
        end
        'z': begin
            ofn = shomedir()+'/rbspb_on_thg/rbspb_on_thg_'+$
                time_string(ut0,tformat='YYYY_MMDD_hhmm_ss')+'.png'
            sgzopen, ofn, xsize = xsz, ysize = ysz, /inch
        end
        'w': begin
            ofn = 0
            sgwopen, ofn, xsize = xsz, ysize = ysz, /inch
        end
    endcase
    
    !p.charsize = 1d
    if n_elements(minlat) eq 0 then minlat = 50
    
    loadct2, 43
    sgtruecolor
    erase, color = sgcolor('white')
    
    timg = asiimg
    
    xr = [-45,45]
    yr = [minlat,70]
    
    isz = max(size(timg,/dimensions))
    xc = smkarthm(-1d,1,isz,'n') # (dblarr(isz)+1)
    yc = smkarthm(-1d,1,isz,'n') ##(dblarr(isz)+1)
    r = sqrt(xc^2+yc^2)
    t = atan(yc, xc)
    lat = 90-r*(90-minlat)              ; in deg.
    lon = (t/!dpi*180+360) mod 360      ; in deg.
    lat = (rotate(lat,1))[*,0:isz/2]
    lon = (rotate(lon,1)-180)[*,0:isz/2]
    
    idx = where(lon ge xr[0] and lon le xr[1] and lat ge yr[0] and lat le yr[1])
    idx = array_indices(lon, idx)
    xidx = minmax(idx[0,*])+[-1,1]
    yidx = minmax(idx[1,*])+[-1,1]
    xidx[0] >= 0 & xidx[1] <= isz-1
    yidx[0] >= 0 & yidx[1] <= isz/2
    
    timg = timg[xidx[0]:xidx[1],yidx[0]:yidx[1]]
    
    imgsz = size(timg,/dimensions)
    asipos = [(1-asixsz/xsz)*0.5,0,(1+asixsz/xsz)*0.5,0.95]
    asipos[1] = asipos[3]-double(asixsz)/xsz*imgsz[1]/imgsz[0]*xsz/ysz
    
;    tv, timg, asipos[0], asipos[1], /normal, /nan
    sgtv, timg, position = asipos, /nan
    xyouts, asipos[0], asipos[1], color = sgcolor('white'), $
        'ASI: '+time_string(ut0,tformat='hh:mm:ss'), /normal
    
    sgset_map, xrange = xr, yrange = yr, pos = asipos, $
        ytickv = smkarthm(yr[0],yr[1],5,'dx'), xtickv = smkarthm(xr[0],xr[1],3,'n'),$
        ytickpos = 30, xtickpos = minlat+1.2, xtickname = ['21','00','03'], xminor = 2, $
        xticknudge = [[0,-0.5],[0,-0.5],[0,-0.5]], yticknudge = [0.75,0], $
        color = sgcolor('white'), charsize = 1, majorgrid = 1, minorgrid = 1
    
    ; plot directly using plots.
    loadct2, 43
    sgtruecolor
    !p.color = sgcolor('black')
    !p.background = sgcolor('white')
    vars = asivars
    nvar = n_elements(vars)
    for i = 0, nvar-1 do begin      ; need to account for midnight longitude.
        get_data, vars[i], t0, dat, limits = lim
        plots, dat[*,0]-midn[0], dat[*,1], color = lim.colors, linestyle = lim.linestyle
    endfor
    
    ; plot using plots after conversion.
    vars = vecvars
    nvar = n_elements(vars)/2
    for i = 0, nvar-1 do begin
        get_data, vars[i,0], t0, dat
        p0 = sinterpol(dat,t0, ut0) & p0[0] -= midn[0]
        get_data, vars[i,1], t0, dat, scl, limits = lim
        v0 = sinterpol(dat,t0, ut0)
        p1 = svec2map(v0, p0, scl, log = log)
        plots, [p0[0],p1[0]], [p0[1],p1[1]], color = lim.colors, linestyle = lim.linestyle
    endfor
    
    ; plot using tplot.
    !p.charsize = 1d
    vars = pltvars
    nvar = n_elements(vars)
    varpos = sgcalcpos(nvar, region = [0d,0,1,asipos[1]])
    tplot, pos = varpos, vars, /noerase, trange = tr
    timebar, ut0, color = sgcolor('red')
    
    case sg of
        'ps': sgpsclose, /pdf
        'z': sgzclose
        'w': sgwclose
    endcase
end

tplot_options, 'labflag', 1
tplot_options, 'num_lab_min', 10
sites = ['atha','tpas']

r = sgcolor('red')
g = sgcolor('lime')
b = sgcolor('blue')
white = sgcolor('white')

etr = stoepoch(['2013-05-01/04:00','2013-05-01/10:00'])

; restore, v_exb, pf.
restore, filename = sdiskdir('Works')+$
    '/works/wygant/pflux_thm_asi/rbspb_2013_0501_pflux_scott/RBSPb_ExB_vels_FAC_11_sec_ave_fields_2013-05-01.dat'
store_data, 'rbspb_vexb', etimes, [[v_azim],[v_rad]], limits = $
    {ytitle:'V_exb!C(km/s)', colors:[r,g], labels:['azim','rad']}
fn = sdiskdir('Works')+'/works/wygant/pflux_thm_asi/'+ $
    '/rbspb_2013_0501_pflux_scott/RBSPb_mapped_pflux_2013-05-01_11-to-55-sec.dat'
restore, filename = fn
store_data, 'rbspb_pf1', etimes, s_para_mapped
options, 'rbspb_pf1', 'ytitle', 'S para!C(mW/m!U2!N)'

; load rbsp e/b field, pos, fpoint, etc.
fn = sdiskdir('Works')+'/data/rbsp_thg/rbsp_efw_fld_2013_0501_04.tplot'
tplot_restore, filename = fn

var = 'rbspb_pf_fac_mat'
options, var, 'labels', ['x','y','z']
options, var, 'ytitle', 'S!C(mW/m!U2!N)'
options, var, 'colors', [r,g,b]
;get_data, var, tmp, dat
;store_data, var+'_para', tmp, dat[*,0], limits = {ytitle:'S para!C(mW/m!U2!N)'}
;ylim, var+'_para', 1e-4, 1, 1

vars = ['rbspb_de_fac','rbspb_db_fac']
options, vars, 'colors', [r,g,b]

; load count rate along s/c footpoint.
vars = ['rbspb_fpt_lonlat_thg_'+sites+'_cr']
ylim, vars, 1, 100, 1

scle = 0.05     ; 1 mV/m ~ 10 deg.
scle = 0.1      ; for log field.
get_data, 'rbspb_de_fac', t0, dat
store_data, 'rbspb_de_vec', t0, dat[*,0:1], scle, limits = {colors:r, linestyle:0}
vecvars = [['rbspb_fpt_lonlat'],['rbspb_de_vec']]

asivar = 'asi'
asivars = 'rbspb_fpt_lonlat'
store_data, 'rbspb_fpt_lonlat', limits = {colors:white, linestyle:1}
pltvars = ['rbspb_'+['de_fac','db_fac','vexb','pf_fac_mat','fpt_lonlat_thg_'+sites+'_cr']]
tr = sfmepoch(etr,'unix')

; read asi image.
mofn = sdiskdir('Research')+'/sdata/themis/thg/mosaic/thg_asf_mosaic_2013_0501_0400_to_1000.cdf'
moid = cdf_open(mofn)
mouts = *((scdfread(moid,'time',skt=skt))[0].value)
mosz = *((scdfread(moid,'image_size',skt=skt))[0].value)
moids = *((scdfread(moid,'pixel_index',skt=skt))[0].value)
minlat = *((scdfread(moid,'minlat',skt=skt))[0].value) & minlat = minlat[0]
midn = *((scdfread(moid,'midn',skt=skt))[0].value)

uts = time_double(['2013-05-01/07:38','2013-05-01/07:38:03'])
uts = smkarthm(uts[0],uts[1],3,'dx')
sg = 'ps'
for i = 0, n_elements(uts)-1 do begin
    tut = uts[i]
    print, 'processing '+time_string(tut)+' ...'
    idx = where(tut eq mouts, cnt)
    if cnt eq 0 then continue
    tmo = *((scdfread(moid,'thg_mosaic',idx,skt=skt))[0].value)
    img = bytarr(mosz)
    img[moids] = tmo
    store_data, 'asi', tut, img, midn[idx]
    vars_on_asi, asivar, asivars = asivars, pltvars = pltvars, vecvars = vecvars, trange = tr, sg = sg, minlat = minlat
endfor
cdf_close, moid
end