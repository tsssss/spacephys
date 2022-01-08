
pro sgzopen, fn, _extra = extra
    compile_opt idl2
    
    dir = file_dirname(fn)
    if ~file_test(dir,/directory) then file_mkdir, dir
    
    defsysv, 'sgraph', exist = flag
    if flag eq 0 then sgwindow, _extra = extra, /nowindow
    !sgraph.zmode.sgid = fn
    mode = !sgraph.zmode
    
    set_plot, mode.device
    !sgraph.d1 = !d
    
    opt = {set_character_size: mode.char, $
        set_resolution: mode.area}
    device, _extra = opt
    
    !p.thick = mode.thick
    !p.charthick = mode.thick
    !x.thick = mode.thick
    !y.thick = mode.thick
    !z.thick = mode.thick
    
    !p.background = sgcolor('white')
    !p.color = sgcolor('black')
    
    erase, color = sgcolor('white')
end