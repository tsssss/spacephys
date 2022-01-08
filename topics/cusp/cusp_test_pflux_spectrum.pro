;+
; study on pflux spectrum on frequency.
;-



rad = !dpi/180
deg = 180/!dpi
re = 6378d

infos = []
       
ids = cusp_calc_id(cusp_id('south_imf'),cusp_id('polar_south_imf'))
foreach id, ids do begin
    tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 1)
    infos = [infos,tinfo]
endforeach

ids = cusp_id('south_imf')
foreach id, ids do begin
    tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 0)
    infos = [infos,tinfo]
endforeach


;; Plot 1: shows v_s/c as a function of distance.
;sgopen, shomedir()+'/cusp_fig_vsc_vs_r.pdf', xsize = 7, ysize = 5, /inch
;plot, podiss, povscs, psym = 1, xtitle = 'R (Re)', ytitle = 'Vs/c (km/s)', $
;    title = 'Vs/c vs R, Polar-black, FAST-red'
;oplot, fadiss, favscs, psym = 1, color = sgcolor('red')
;sgclose

;; Plot 2: shows the lower limit converted to spatial scale in km.
;sgopen, shomedir()+'/cusp_fig_min_scale_vs_r.pdf', xsize = 7, ysize = 5, /inch
;plot, podiss, pohfreqlims, psym = 1, yrange = [1,100], /ylog, $
;    xtitle = 'R (Re)', ytitle = 'Min Scale (km)', $
;    title = 'Min Scale vs R, Polar-black, FAST-red'
;oplot, fadiss, fahfreqlims, psym = 1, color = sgcolor('red')
;sgclose
;
;; Plot 3: shows the upper limit converted to deg in ILat.
;sgopen, shomedir()+'/cusp_fig_max_scale_vs_r.pdf', xsize = 7, ysize = 5, /inch
;plot, podiss, polfreqlims, psym = 1, yrange = [1,100], /ylog, $
;    xtitle = 'R (Re)', ytitle = 'Max Scale (deg)', $
;    title = 'Max Scale vs R, Polar-black, FAST-red'
;oplot, fadiss, falfreqlims, psym = 1, color = sgcolor('red')
;sgclose
;
;; Plot 4: shows the fac limit converted to deg in ILat.
;sgopen, shomedir()+'/cusp_fig_fac_scale_vs_r.pdf', xsize = 7, ysize = 5, /inch
;plot, podiss, pofaclims, psym = 1, yrange = [1,100], /ylog, $
;    xtitle = 'R (Re)', ytitle = 'FAC Scale (deg)', $
;    title = 'FAC Scale vs R, Polar-black, FAST-red'
;oplot, fadiss, fafaclims, psym = 1, color = sgcolor('red')
;sgclose

;; Plot 5: this is the pflux within original filters.
;sgopen, shomedir()+'/cusp_fig_pflux0_vs_r.pdf', xsize = 7, ysize = 5, /inch
;plot, podiss, popflux0, psym = 1, yrange = [1,1e4], /ylog, $
;    xtitle = 'R (Re)', ytitle = 'FAC Scale (deg)', $
;    title = 'FAC Scale vs R, Polar-black, FAST-red'
;oplot, fadiss, fapflux0, psym = 1, color = sgcolor('red')
;sgclose

;ofn = shomedir()+'/cusp_fig_scale_vs_r.pdf'
;;ofn = 0
;sgopen, ofn, xsize = 7, ysize = 5, /inch
;
;tpos = [0.1,0.1,0.8,0.9]
;plot, infos.dis, infos.lfreqlim, psym = 1, /ylog, yrange = [.005,50], $
;    position = tpos, title = 'Scale vs R, Polar and FAST orbits', symsize = 0.8, $
;    xrange = [1,6], xstyle = 1, ystyle = 1, xtitle = 'R (Re)', ytitle = 'Scale (deg)'
;oplot, infos.dis, infos.faclim, psym = 4, color = sgcolor('red'), symsize = 0.8
;oplot, infos.dis, infos.hfreqlim, psym = 4, color = sgcolor('blue'), symsize = 0.8
;
;xchsz = double(!d.x_ch_size)/!d.x_size
;ychsz = double(!d.y_ch_size)/!d.y_size
;
;tx = tpos[2]+2*xchsz & ty = tpos[3]-ychsz*2
;ttext = 'Cusp Scale' & tsym = 1 & tcolor = sgcolor('black')
;xyouts, tx, ty, /normal, ttext
;plots, tx+strlen(ttext)*xchsz, ty+0.2*ychsz, /normal, psym = tsym, color = tcolor
;
;tx = tpos[2]+2*xchsz & ty = tpos[3]-ychsz*4
;ttext = 'Max Scale' & tsym = 4 & tcolor = sgcolor('red')
;xyouts, tx, ty, /normal, ttext
;plots, tx+strlen(ttext)*xchsz, ty+0.2*ychsz, /normal, psym = tsym, color = tcolor
;
;tx = tpos[2]+2*xchsz & ty = tpos[3]-ychsz*6
;ttext = 'Min Scale' & tsym = 4 & tcolor = sgcolor('blue')
;xyouts, tx, ty, /normal, ttext
;plots, tx+strlen(ttext)*xchsz, ty+0.2*ychsz, /normal, psym = tsym, color = tcolor
;
;sgclose


;; Plot 6: pflux on dis, scale plane.
;ofn = 0
;ofn = shomedir()+'/cusp_fig_pflux_alt_spectrogram.pdf'
;sgopen, ofn, xsize = 7, ysize = 5, /inch
;tpos = [0.1,0.1,0.8,0.9]
;
;ninfo = n_elements(infos)
;plot, [1,6], [0.05,20], /ylog, /nodata, position = tpos, xstyle = 1, ystyle = 1, $
;    xtitle = 'R (Re)', ytitle = 'Scale (deg)'
;    
;zr = [0.01,5e2]
;ncolors = 255
;ct = 40
;
;for i = 0, ninfo-1 do begin
;    tdis = infos[i].dis
;    nfilter = infos[i].nfilter
;    filters = infos[i].filters[0:nfilter-1]*infos.vsc/tdis/re*deg
;    pfluxes = infos[i].pfluxes[0:nfilter-1]
;    tcolors = bytscl(pfluxes, min = zr[0], max = zr[1], top = ncolors-1)
;    tradius = sqrt(tcolors/256d)
;
;    for j = 0, nfilter-1 do begin
;        txs = cos(findgen(11)*2*!dpi/10)*tradius[j]
;        tys = sin(findgen(11)*2*!dpi/10)*tradius[j]
;        usersym, txs, tys, color = sgcolor(tcolors[j], ct = ct), /fill
;        plots, tdis, filters[j], psym = 8
;    endfor
;endfor
;
;tpos[0] = 0.8+2*xchsz & tpos[2] = tpos[0]+1*xchsz
;sgcolorbar, indgen(ncolors), position = tpos, zrange = zr, ct = ct, $
;    ztitle = 'Line-integrated S!D||!N(W/m)'
;
;sgclose

;; Plot 7: pflux overlap on each other.
;ofn = 0
;ofn = shomedir()+'/cusp_fig_pflux_spectrum.pdf'
;sgopen, ofn, xsize = 7, ysize = 5, /inch
;tpos = [0.15,0.1,0.85,0.9]
;ninfo = n_elements(infos)
;plot, [1,1e4], [0.01,1e4], /nodata, xstyle = 1, ystyle = 1, position = tpos, $
;    xlog = 1, ylog = 1, ytitle = 'Line-integrated S!D||!N(W/m)', xtitle = 'Scale (sec)'
;for i = 0, ninfo-1 do begin
;    nfilter = infos[i].nfilter
;    filters = infos[i].filters[0:nfilter-1]
;    pfluxes = infos[i].pfluxes[0:nfilter-1]
;    plots, filters, pfluxes, psym = 1
;;    plots, filters, pfluxes, psym = -1, linestyle = 1
;endfor
;sgclose

;; Plot 8: pflux vs dis.
;ofn = 0
;ofn = shomedir()+'/cusp_fig_pflux_vs_r.pdf'
;sgopen, ofn, xsize = 7, ysize = 5, /inch
;tpos = [0.1,0.1,0.8,0.9]
;plot, infos.dis, infos.pflux, psym = 1, ylog = 1, yrange = [2,2e4], xrange = [1,6], $
;    xtitle = 'R (Re)', ytitle = 'Line-integrated S!D||!N(W/m)', xstyle = 1, ystyle = 1
;sgclose

end
