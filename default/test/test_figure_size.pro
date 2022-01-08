;+
; Type: procedure.
; Purpose: Create similar size figure in both postscrip and window modes.
; Parameters: none.
; Keywords: none.
; Notes: none.z
; Dependence: none.
; History:
;   2014-01-23, Sheng Tian, create.
;-

pro test_figure_size

    ; common settings.
    !p.font = 1
    px2cm = 1d/40d
    xsize = 1900/2 & ysize = 1000/2     ; in pixel.
    figaspect = double(ysize)/xsize
    dev0 = (!version.os_family eq 'Windows')? 'win': 'x'

    set_plot, 'ps'
    print, 'ps'
    print, !d.x_px_cm, !d.y_px_cm
    print, !d.x_size, !d.y_size, !d.x_size/!d.x_px_cm, !d.y_size/!d.y_px_cm, double(!d.y_size)/!d.x_size
    print, !d.x_ch_size, !d.y_ch_size, double(!d.y_ch_size)/!d.x_ch_size, double(!d.y_size)/!d.y_ch_size
    set_plot, dev0
    print, dev0
    print, !d.x_px_cm, !d.y_px_cm
    print, !d.x_size, !d.y_size, !d.x_size/!d.x_px_cm, !d.y_size/!d.y_px_cm, double(!d.y_size)/!d.x_size
    print, !d.x_ch_size, !d.y_ch_size, double(!d.y_ch_size)/!d.x_ch_size, double(!d.y_size)/!d.y_ch_size
    set_plot, 'z'
    print, 'z'
    print, !d.x_px_cm, !d.y_px_cm
    print, !d.x_size, !d.y_size, !d.x_size/!d.x_px_cm, !d.y_size/!d.y_px_cm, double(!d.y_size)/!d.x_size
    print, !d.x_ch_size, !d.y_ch_size, double(!d.y_ch_size)/!d.x_ch_size, double(!d.y_size)/!d.y_ch_size
    
    ; window mode.
    set_plot, dev0
    figid = 2
    figxsz = xsize
    figysz = ysize
    charysz = 14.08
    window, figid, xsize = figxsz, ysize = figysz
    for i = 0, figysz/charysz do xyouts, 0, i*charysz, string(i), /device
    
    ; postscript mode.
    set_plot, 'ps'
    figid = shomedir()+'/test_figure_size_'+dev0+'.ps'
    figxsz = xsize*px2cm
    figysz = ysize*px2cm
    charysz = !d.y_ch_size
    device, filename = figid, xsize = figxsz, ysize = figysz, /cm
    for i = 0, figysz*!d.y_px_cm/charysz do xyouts, 0, i*charysz, string(i), /device
    device, /close & set_plot, dev0
    
end