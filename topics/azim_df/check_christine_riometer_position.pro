;+
; Check the locations of riometers.
;-

rio_info = dictionary($
    'daws', [64.05,220.89], $
    'fsim', [61.76,238.77], $
    'fsmi', [60.02,248.05], $
    'rabb', [58.22,256.32], $
    'gill', [56.38,265.36], $
    'mcmu', [56.66,248.79], $
    'isll', [53.86,265.34])
time_range = time_double(['2014-08-28/10:10','2014-08-28/10:20'])
times = mean(time_range)
ntime = n_elements(times)
ndim = 3

rios = rio_info.keys()
nrio = n_elements(rios)
rad = constant('rad')
deg = constant('deg')
mlts = fltarr(nrio,ntime)
mlats = fltarr(nrio,ntime)
foreach rio, rios, rio_id do begin
    glatlon = rio_info[rio]
    glat = glatlon[0]*rad
    glon = glatlon[1]*rad
    r_geo = [cos(glat)*cos(glon),cos(glat)*sin(glon),sin(glat)]
    r_mag = dblarr(ntime,ndim)
    foreach time, times, time_id do begin
        r_mag[time_id,*] = cotran(r_geo, times, 'geo2mag')
    endforeach
    mlon = atan(r_mag[*,1], r_mag[*,0])*deg
    mlt = mlon2mlt(mlon, times)
    mlts[rio_id,*] = mlt
    mlats[rio_id,*] = asin(r_mag[*,2])*deg
endforeach

index = sort(mlts[*,0])
mlts = mlts[index,*]
mlats = mlats[index,*]
rios = rios[index]
strlen = 10
foreach rio, rios, rio_id do begin
    print, extend_string(rio,length=strlen)+strjoin(string(reform(mlts[rio_id,*])))
endforeach


; Load SC MLT.
probes = ['thd','g15','1991-080','g13','LANL-01A']
if check_if_update(probes[0]+'_mlt') then _2014_0828_10_load_data
nprobe = n_elements(probes)
sc_mlts = fltarr(nprobe,ntime)
sc_mlats = fltarr(nprobe,ntime)
foreach probe, probes, probe_id do begin
    sc_mlts[probe_id,*] = get_var_data(probe+'_mlt', at=times)
    sc_mlats[probe_id,*] = get_var_data(probe+'_fmlat_t89', at=times)
endforeach

index = sort(sc_mlts[*,0])
sc_mlts = sc_mlts[index,*]
sc_mlats = sc_mlats[index,*]
probes = probes[index]

foreach probe, probes, probe_id do begin
    print, extend_string(probe,length=strlen)+strjoin(string(reform(sc_mlts[probe_id,*])))
endforeach


; Plot.
test = 1
plot_file = join_path([homedir(),'fig_2014_0828_sc_and_rio_mlt.pdf'])
if keyword_set(test) then plot_file = 0

margins = [0,0,8,0]
;poss = fig_init_simple(plot_file, pansize=[6,6], size=fig_size, margins=margins)
poss = panel_pos(plot_file, pansize=[6,6], fig_size=fig_size, margins=margins)


sgopen, plot_file, xsize=fig_size[0], ysize=fig_size[1], /inch, xchsz=xchsz, ychsz=ychsz

xticklen_chsz = -0.15
yticklen_chsz = -0.30


tpos = poss
xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

plot, [-1,1], [-1,1], $
    xstyle=5, ystyle=5, $
    position=tpos, nodata=1, noerase=1

    
    rr = 1
    tmp = smkarthm(0,2*!dpi,50,'n')
    xxs = rr*cos(tmp)
    yys = rr*sin(tmp)
    plots, xxs, yys, /data, linestyle=1

psym = 1
min_mlat = 55   ; deg.
xtickv = make_bins([min_mlat,90], 10, /inner)
foreach mlat, xtickv do begin
    rr = (90-mlat)/(90-min_mlat)
    tmp = smkarthm(0,2*!dpi,50,'n')
    xxs = rr*cos(tmp)
    yys = rr*sin(tmp)
    plots, xxs, yys, /data, linestyle=1
    ang = !dpi*0.5
    tx = rr*cos(ang)
    ty = rr*sin(ang)
    xyouts, tx,ty,/data, string(mlat,format='(I0)')
endforeach

foreach ang, make_bins([0,2*!dpi],!dpi*0.25) do begin
    rr = [0,1]
    xxs = rr*cos(ang)
    yys = rr*sin(ang)
    plots, xxs, yys, /data, linestyle=1
endforeach

; Add earth.
tmp = 50
tmp = findgen(tmp)/(tmp-1)*2*!dpi
rr = 0.1
xs = cos(tmp)*rr
ys = sin(tmp)*rr
polyfill, xs>0, ys, color=sgcolor('silver')
polyfill, xs<0, ys, color=sgcolor('white')
plots, xs, ys


rr = 0.95
foreach mlt, make_bins([0,6],3) do begin
    ang = mlt*15*rad
    tx = rr*cos(ang)
    ty = rr*sin(ang)
    xyouts, tx,ty,/data, string(mlt,format='(I02)')
endforeach

red = sgcolor('red')
foreach rio, rios, rio_id do begin
    mlt = mlts[rio_id]
    tt = (mlt)*15*rad
    rr = (90-mlats[rio_id])/(90-min_mlat)
    tx = cos(tt)*rr
    ty = sin(tt)*rr
    plots, tx,ty,/data, psym=psym, color=red
    tmp = convert_coord([tx,ty], /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx+xchsz*0.5,ty,/norm, strupcase(rio)
endforeach

foreach probe, probes, probe_id do begin
    mlt = sc_mlts[probe_id]
    tt = (mlt)*15*rad
    rr = (90-mlats[probe_id])/(90-min_mlat)
    tx = cos(tt)*rr
    ty = sin(tt)*rr
    plots, tx,ty,/data, psym=psym, color=red
    tmp = convert_coord([tx,ty], /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    if probe eq '1991-080' then begin
        tx = tmp[0]-xchsz*3
        ty = tmp[1]+ychsz*0.5
    endif
    xyouts, tx+xchsz*0.5,ty,/norm, strupcase(probe)
endforeach



    ; Add labels.
    tx = xchsz*2
    ty = tpos[3]-2*ychsz
    tab1 = 8
    tab2 = 16
    xyouts, tx,ty, /normal, 'Rio&SC'
    xyouts, tx+xchsz*tab1,ty, /normal, 'MLat (deg)'
    xyouts, tx+xchsz*tab2,ty, /normal, 'MLT (hr)'
        
    foreach rio, rios, rio_id do begin
        ty = tpos[3]-(rio_id+3)*ychsz
        xyouts, tx,ty, /normal, strupcase(rio)
        xyouts, tx+xchsz*tab1,ty, /normal, $
            string(mlats[rio_id],format='(F4.1)')
        xyouts, tx+xchsz*tab2,ty, /normal, $
            string(mlts[rio_id],format='(F4.1)')
    endforeach

    foreach probe, probes, probe_id do begin
        ty = tpos[3]-(probe_id+nrio+5)*ychsz
        xyouts, tx,ty, /normal, strupcase(probe)
        xyouts, tx+xchsz*tab1,ty, /normal, $
            string(sc_mlats[probe_id],format='(F4.1)')
        xyouts, tx+xchsz*tab2,ty, /normal, $
            string(sc_mlts[probe_id],format='(F4.1)')
    endforeach


if keyword_set(test) then stop
sgclose


end
