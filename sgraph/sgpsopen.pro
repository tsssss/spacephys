
pro sgpsopen, fn, _extra = extra
    compile_opt idl2
    
    dir = file_dirname(fn)
    if ~file_test(dir,/directory) then file_mkdir, dir
    
    defsysv, 'sgraph', exist = flag
    if flag eq 0 then sgwindow, _extra = extra, /nowindow
    !sgraph.pmode.sgid = fn
    mode = !sgraph.pmode
    
    set_plot, mode.device, /copy, /interpolate
    !sgraph.d1 = !d
    
    opt = {filename: mode.sgid, $
        set_character_size: mode.char, $
        encapsulate: 1, copy:1, interpolate:1, $
        xsize: mode.area[0], ysize: mode.area[1], cm: 1}
    device, _extra = opt
    
    !p.thick = mode.thick
    !p.charthick = mode.thick
    !x.thick = mode.thick
    !y.thick = mode.thick
    !z.thick = mode.thick

end