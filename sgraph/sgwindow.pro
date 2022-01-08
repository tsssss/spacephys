;+
; Type: procedure.
; Purpose: Create graphics area to plot.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro sgwindow, sgidw, xsize = xsize, ysize = ysize, $
    cm = cm, inch = inch, magnify = magcoef, nowindow = nowindow, _extra = extra
    
    ; original device and plot.
    d0 = !d

    ; window size.
    if n_elements(xsize) eq 0 then xsz = !d.x_size else xsz = double(xsize)
    if n_elements(ysize) eq 0 then ysz = !d.y_size else ysz = double(ysize)
    xsz = double(xsz) & ysz = double(ysz)
;    pchsz = [225d,375] & zchsz = [9d,15] & wchsz = [9d,15]
    pchsz = [216d,360] & zchsz = [9d,15] & wchsz = [9d,15]   
        
    if n_elements(magcoef) ne 0 then begin
        magcoef = double(magcoef)
        xsz *= magcoef & ysz *= magcoef
        pchsz *= magcoef & zchsz *= magcoef & wchsz *= magcoef
    endif
    
    cm2px = 40d & inch2cm = 2.54d & inch2px = 101.6d ; inch2cm*cm2px.
    if keyword_set(cm) then begin
        if keyword_set(inch) then message, 'cannot set inch and cm both ...'
        xsz *= cm2px & ysz *= cm2px
    endif else if keyword_set(inch) then begin
        xsz *= inch2px & ysz *= inch2px
    endif
    
    ; draw window.
    devw = (!version.os_family eq 'unix')?'x':'win'
    set_plot, devw
    if n_elements(sgidw) eq 0 then begin
        window, xsize = xsz, ysize = ysz, /free, /pixmap, _extra = extra
        sgidw = !d.window
    endif else window, sgidw, xsize = xsz, ysize = ysz, _extra = extra
    erase, color = sgcolor('white')
    if keyword_set(nowindow) then wdelete, sgidw
    
    ; settings that smears out device differences.
    !p.font = 1         ; font system, use truetype.
;    !p.charsize = 0.9   ; overall scale so characters are around 10-12pt.
    
    ; os coef. mac: 1.045; linux: 1.412; windows: 0.925.
    ; !p.charsize *= oscoef, pmode.area *= oscoef.
    
    ; cm to pixel.
    px2cm = 0.025d
;    set_plot, 'ps' & px2cm = 1d/!d.x_px_cm
    set_plot, d0.name
    pmode = {sgid: filepath('idl.ps', root_dir = shomedir()), device: 'ps', $
        area: [xsz,ysz]*px2cm, char: pchsz, thick: 2d}
    zmode = {sgid: filepath('idl.png', root_dir = shomedir()), device: 'z', $
        area: [xsz,ysz], char: zchsz, thick: 1d}
    wmode = {sgid: fix(sgidw), device: devw, $
        area: [xsz,ysz], char: wchsz, thick: 1d}

    defsysv, '!sgraph', exists = flag
    if flag eq 0 then begin
        sgraph = {d0:d0, d1:!d, p0:!p, x0:!x, y0:!y, z0:!z, $
            pmode:pmode, zmode:zmode, wmode:wmode}
        defsysv, '!sgraph', sgraph
    endif else begin
        !sgraph.pmode = pmode
        !sgraph.zmode = zmode
        !sgraph.wmode = wmode
    endelse
end

sgwindow, 0, xsize = 600, ysize = 400
device, set_character_size = !sgraph.wmode.char
!p.background = !p.color & !p.color = 0     ; this is dangerous.
plot, loaddata(1), xtitle = 'latitude', ytitle = 'speed', title = 'Hello'

set_plot, !sgraph.pmode.device
device, filename = !sgraph.pmode.sgid, $
    set_character_size = !sgraph.pmode.char, $
    xsize = !sgraph.pmode.area[0], ysize = !sgraph.pmode.area[1]
plot, loaddata(1), xtitle = 'latitude', ytitle = 'speed', title = 'Hello'
device, /close
set_plot, !sgraph.d0.name
end
