

infofn = shomedir()+'/Google Drive/works/data/rbsp_de/32hz_event_info.tplot'
tplot_restore, filename = infofn
tvar = 'psbl_de_32hz_info'
get_data, tvar, tmp, infos

ninfo = n_elements(infos)


mlats = infos.fptmlat
mlts = infos.fptmlt
mlts[where(mlts gt 12)]-= 24


; the conditions.
maxdes = snorm(transpose(infos.maxde.de_fac))
idx = where(maxdes gt 0)


if cnt eq 0 then message, 'no event found ...'
mlats = mlats[idx]
mlts = mlts[idx]
infos = infos[idx]
ninfo = n_elements(infos)



types = ['e']
maxv = 3d/10   ; deg/mV/m.
maxv = 4d/10   ; deg/mV/m.
pos = [0.15,0.15,0.9,0.9]
xr = [-6,6]*15
yr = [-1,1]*90


foreach type, types do begin
    
    udcnt = intarr(2,2,2)   ; mlt, mlat, n/s.
    
    ofn = shomedir()+'/fig_e_equadrant.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    
    ; set the box and labels.
    plot, xr/15, yr, xrange = xr/15, yrange = yr, /nodata, position = pos, $
        noerase = noerase, $
        xstyle = 1, ystyle = 1, xtitle = 'MLT (hr)', ytitle = 'Mapped MLat (deg)', $
        title = 'RBSP dE at PSBL'
        
    ; set the data coord.
    plot, xr, yr, xrange = xr, yrange = yr, /nodata, /noerase, position = pos, $
        xstyle = 5, ystyle = 5
    plots, mean(xr)+[0,0], yr, linestyle = 1
    plots, xr, 65+[0,0], linestyle = 1
    plots, xr,-65+[0,0], linestyle = 1
    
    for i = 0, ninfo-1 do begin
        v0 = [mlts[i]*15,mlats[i]]
        v1 = infos[i].maxde.de_fac
        print, infos[i].id, snorm(v1)
        v1 = v1[1:2]
        udcnt[mlts[i] gt 0,mlats[i] gt 0,v1[1] gt 0]++
        v1 = v0+v1*maxv
        arrow, v0[0],v0[1], v1[0],v1[1], /data, color = sgcolor('blue'), hsize = 50, /solid
    endfor
    
    plots, mean(xr)+maxv*0.5*[-1,1]*50, mean(yr)+[0,0], color = sgcolor('blue')
    xyouts, mean(xr), mean(yr)-10, /data, '50 mV/m', alignment = 0.5
    
    ; lower left.
    xyouts,-3*15,-40, /data, alignment = 0.5, $
        sgnum2str(udcnt[0,0,0])+'D:'+sgnum2str(udcnt[0,0,1])+'U'
    ; lower right.
    xyouts, 3*15,-40, /data, alignment = 0.5, $
        sgnum2str(udcnt[1,0,0])+'D:'+sgnum2str(udcnt[1,0,1])+'U'
    ; upper left.
    xyouts,-3*15, 40, /data, alignment = 0.5, $
        sgnum2str(udcnt[0,1,0])+'D:'+sgnum2str(udcnt[0,1,1])+'U'
    ; upper right.
    xyouts, 3*15, 40, /data, alignment = 0.5, $
        sgnum2str(udcnt[1,1,0])+'D:'+sgnum2str(udcnt[1,1,1])+'U'
        
    sgclose

endforeach

end
