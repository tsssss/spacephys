

utr = time_double(['2013-06-07/04:30','2013-06-07/05:30'])
;utr = time_double(['2013-06-07/05:00','2013-06-07/05:30'])
sites = ['pina','kapu','chbg']
probes = ['a','b']
syms = [1,1]
colors = [6,255]

deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1

reload = 1


; load high res B field.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    tvar = pre0+'b'
    load = tnames(tvar) eq ''
    if reload then load
    if load then begin
        emfisis = sread_rbsp_emfisis_l3(utr, type = 'hires', probe = tprobe)
        store_data, tvar, sfmepoch(emfisis.epoch,'unix'), emfisis.mag, $
            limit = {ytitle:'B GSM!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
    endif
endforeach

; load high res E field.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    tvar = pre0+'e'
    if tnames(tvar) eq '' then begin
        efw = sread_rbsp_efw_l2(utr, type = 'euvw', probe = tprobe)
        efw.efield_uvw[*,2] = 0
        store_data, tvar, sfmepoch(efw.epoch,'unix'), efw.efield_uvw, $
            limit = {ytitle:'E UVW!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
    endif
endforeach

stop


; load pitch 2d.
types = ['electron','proton']
log = 1
unit = 'energy'

uts = smkarthm(utr[0],utr[1],12,'dx')

foreach tprobe, probes do begin
    hopel3 = sread_rbsp_hope_l3(utr, probes = tprobe)
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l3_pitch2d, uts[i], types[j], unit = unit, $
                log = log, hopel3 = hopel3, probe = tprobe
endforeach

stop

; load pina asf.
tvar = 'pina_asf'
if tnames(tvar) eq '' then begin
    asi = sread_thg_asi(utr, 'pina', type = 'asf')
    store_data, tvar, asi.utsec, asi.img
endif


xsz = 800
ysz = xsz
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

get_data, tvar, uts, imgs
nrec = n_elements(uts)
for i = 0, nrec-1 do begin

    ofn = shomedir()+'/thg_asi/pina/thg_asf_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
    ;    ofn = 0
    sgopen, ofn, xsize = xsz, ysize = ysz
    
    device, decomposed = 0
    loadct, 39
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    timg = reform(imgs[i,*,*])
    timg *= 64d/(median(timg) > 1)
    timg = bytscl(timg)
    timg = congrid(timg, xsz, ysz)
    
    sgtv, timg, position = tpos
    
    xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
    
    sgclose
endfor

stop

; load mltimg.
tvar = 'mltimg'
if tnames(tvar) eq '' then begin
    asi = sread_thg_mlt(utr, sites, /half, type = 'asf')
    store_data, tvar, sfmepoch(asi.epoch, 'unix'), asi.mltimg
endif


; load sc position.
tvar = 'rbspa_fpt_mlt'
if tnames(tvar) eq '' then begin
    rbsp_load_spice_kernels, trange = utr
    uts = smkarthm(utr[0], utr[1], 60, 'dx')    ; 1 min resolution.
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        rbsp_load_spice_state, probe = tprobe, coord = 'gsm', times = uts, /no_spice_load
        get_data, pre0+'state_pos_gsm', tmp, pos0
        pos0 = pos0*re1
        store_data, pre0+'state_*', /delete
        
        ets = stoepoch(uts, 'unix')
        pos1 = pos0
        ; loop for each time.
        for j = 0, n_elements(uts)-1 do begin
            ; set geopack.
            geopack_epoch, ets[j], yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001D, /date
            ; pos in gse, which is the mapping coord.
            x1 = pos0[j,0] & y1 = pos0[j,1] & z1 = pos0[j,2]
            dir = (z1 gt 0)? -1: 1
            par = 2
            geopack_trace, x1, y1, z1, dir, par, xf, yf, zf, $
                epoch = ets[j], /refine, /ionosphere, /t89
            ; convert from gse to mag.
            geopack_conv_coord, xf, yf, zf, /from_gsm, $
                x1, y1, z1, /to_mag
            pos1[j,*] = [x1,y1,z1]
        endfor
        mlat = asin(pos1[*,2]*(1/r0))*deg
        mlon = atan(pos1[*,1],pos1[*,0])*deg
        mlt = slon2lt(mlon, stoepoch(uts,'unix'), /mag, /deg)/15    ; in hour.
        mlt = (mlt+24) mod 24
        store_data, pre0+'fpt_mlon', uts, mlon, $
            limits = {ytitle:'MLon/fpt (deg)'}
        store_data, pre0+'fpt_mlat', uts, mlat, $
            limits = {ytitle:'MLat/fpt (deg)'}
        store_data, pre0+'fpt_mlt', uts, mlt, $
            limits = {ytitle:'MLT/fpt (hr)'}
            
    endforeach
    rbsp_load_spice_kernels, trange = utr, /unload
endif



xsz = 800
ysz = xsz*0.5
tpos = [0.1,0.1,0.9,0.9]
white = 255
xr = [-90,90]
chsz = 1.5

get_data, 'mltimg', uts, imgs
nrec = n_elements(uts)
for i = 0, nrec-1 do begin
    
    ofn = shomedir()+'/thg_asi/thg_rb_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
;    ofn = 0
    sgopen, ofn, xsize = xsz, ysize = ysz
    
    device, decomposed = 0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    timg = reform(imgs[i,*,*])
    timg = congrid(timg, xsz, ysz)
    
    sgtv, timg, position = tpos
    sgset_map, position = tpos, color = white, xrange = xr
    
    xyouts, tpos[0], tpos[1], /normal, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz

    for j = 0, n_elements(probes)-1 do begin
        pre0 = 'rbsp'+probes[j]+'_'
        get_data, pre0+'fpt_mlat', tuts, mlat
        get_data, pre0+'fpt_mlt', tuts, mlt
        mlat = interpol(mlat, tuts, uts[i], /quadratic)
        mlt = interpol(mlt, tuts, uts[i], /quadratic)
        plots, mlt*15, mlat, psym = syms[j], color = colors[j]
        xyouts, tpos[0], tpos[1]+ychsz*(j+1), /normal, 'RBSP-'+strupcase(probes[j]), color = colors[j], charsize = chsz

    endfor
        
    sgclose
endfor

end