;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro mat_eg_step_fourier, p

    pistr = '!9'+string(112b)+'!X'
    
    nrec = 5000
    np = 4
    colors = ['black','magenta','blue','cyan','green','yellow','red']

    ps = nrec/smkarthm(1,2,np,'x0')         ; wave periods in # of record.
    xr = [0,nrec]
    yr = [min(ps)/3>5,max(ps)*1.2]
    zr = 4/!dpi*1.1*[-1,1]
    zr = 1.5*[-1,1]
    nztick = 20
    zticks = smkarthm(zr[0],zr[1],nztick,'n')
    ytickv = [1,0.2,0.05]*nrec
    ytickn = reverse(['2','0.4','0.1'])+pistr
    xticks = 4
    xminor = 10
    

    nx = nrec
    x = (dindgen(nx)/(nx-1)-0.5)*2*!dpi     ; [-pi,pi].
    

    lims = findgen(np)*2+1        ; odd.
    z = replicate(0,nrec)
    
    ; prepare spectrogram limit.
    zm = 4/!dpi*1.1
    fm = 0.5
;    ylim = [2*!dpi/max(lims)/2,10]    ; to show first few terms better.
    zlim = [-zm,zm]
    
    ct = 66
    
    p = plot3d(xr, yr, zr, axis_style = 2, $
        /nodata, /overplot, /current, /perspective, /ylog, $
        depth_cue = [2,2], $
        xtitle = 'Time', ytitle = 'Period', ztitle = 'Signal', $
        xtickformat = '(A1)', ztickformat = '(A1)', $
;        xy_shadow = 1, yz_shadow = 1, xz_shadow = 1, $
        xminor = xminor, xmajor = xticks+1, yminor = 0, ymajor = np, zmajor = 5, zminor = 5, $
        xrange = xr, yrange = yr, zrange = zr)
    
    p.rotate, -8, /zaxis
    p.rotate, 25, /xaxis
    
    
    for i = 0, np-1 do begin
        f = 4/!dpi/(2*i+1)*sin((2*i+1)*x)
        w = replicate(ps[i],nrec)
        p = plot3d(findgen(nrec), w, f, 'o', /sym_filled, sym_size = 0.9, $
            /overplot, /current, $
;            xy_shadow = 1, yz_shadow = 1, xz_shadow = 1, $
            rgb_table = ct, vert_colors = ((-fm>f<fm)/fm+1)*127)
    endfor
    
    for i = 0, np-1 do begin
        p = plot3d(xr, ps[i]+[0,0], zr[0]+[0,0], linestyle = 2, thick = 2, $
            color=(colors[i]), /overplot, /current)
        p = plot3d(xr[1]+[0,0], ps[i]+[0,0], zr, linestyle = 2, thick = 2, $
            color=(colors[i]), /overplot, /current)
    endfor
    
    
    ax = p.axes

    idx = [1,4,7,10]
    for i = 0, n_elements(idx)-1 do begin
        ax[idx[i]].xrange = xr
        ax[idx[i]].major = 2
        ax[idx[i]].tickvalue = ytickv
        ax[idx[i]].minor = 10
    endfor
    ax[2].showtext = 0
    ax[8].showtext = 8
    ax[1].tickname = ytickn
    
    idx = [2,6,7]
    for i = 0, n_elements(idx)-1 do ax[idx[i]].hide = 1

end


pro fig_step_2d3d

    mat_eg_step_fourier, p
    ofn = shomedir()+'/step_3d.eps'
    p.save, ofn
    sps2pdf, ofn, /rm
    
    stop
    
    ;p.translate, 0, 0, -0.1, /normal
    p.perspective = 0
    p.rotate, reset = 1
    p.aspect_ratio = 1200
    p.font_size = 10
    
    scale = 1.6
    p.scale, scale, scale, scale
    
    ax = p.axes
    for i = 0, n_elements(ax)-1 do ax[i].ticklen = 0.02
    
    idx = [2,6,7]
    for i = 0, n_elements(idx)-1 do ax[idx[i]].hide = 0
;    ax[6].ticklen = 0.015   ; up x-axis.
;    ax[9].ticklen = 0.015   ; down x-axis.
;    ax[7].ticklen = 0.015   ; left y-axis.
;    ax[10].ticklen = 0.015  ; right y-axis.

    ofn = shomedir()+'/step_2d.eps'
    p.save, ofn
    sps2pdf, ofn, /rm
    
    obj_destroy, p
    


end
