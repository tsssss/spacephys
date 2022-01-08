;+
; Type: procedure.
; Purpose: Test uniform display between PostScript and window.
; Parameters: none.
; Keywords: none.
; Notes: none.
; Dependence: none.
; History:
;   2014-01-23, Sheng Tian, create.
;-

pro test_uniform_display
    ; aspect ratio is useful.
    ; letter paper size: 8.5x11 inch.
    ; !d system variable.
    ; relavent units are pixel, cm, inch.
    
    ; os.
    case !version.os of
        'linux': os = 'linux'
        'darwin': os = 'mac'
        'Win32': os = 'win'
    endcase

    devp = 'ps'
    devw = (os eq 'win')? 'win': 'x'

    ; character size, aspect ratio 0.6.
    charp = [225,375]
    charw = [9,15]      ; charp*1000/40

    ; handler.
    sgidp = shomedir()+'/test_uniform_display_'+os+'_'+devp+'.ps'
    sgidw = 1

    ; default size.
    px2cm = 0.025d
    wsize = [400,500]
    psize = wsize*px2cm
    
    ; test output string.
    str0 = 'MMMMMM'

    ; postscript mode.
    dev0 = !d.name & set_plot, devp
    !p.font = 1 & !p.charsize = 1
    device, set_character_size = charp
    device, filename = sgidp, xsize = psize[0], ysize = psize[1], /cm
    ;for i = 0, !d.y_size, charp[1] do xyouts, 0, i, str0, /device
    plot, loaddata(1), xtitle = 'latitude', ytitle = 'speed', title = 'Hello'
    device, /close & set_plot, dev0
    
    ; window mode.
    dev0 = !d.name & set_plot, devw
    !p.font = 1 & !p.charsize = 1
    device, set_character_size = charw
    window, sgidw, xsize = wsize[0], ysize = wsize[1]
    ;for i = 0, !d.y_size, charw[1] do xyouts, 0, i, str0, /device
    plot, loaddata(1), xtitle = 'latitude', ytitle = 'speed', title = 'Hello'
    set_plot, dev0

end
