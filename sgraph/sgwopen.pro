
pro sgwopen, wid, _extra = extra
    compile_opt idl2
    
    defsysv, 'sgraph', exist = flag
    if flag eq 0 then sgwindow, wid, _extra = extra
    !sgraph.wmode.sgid = wid
    mode = !sgraph.wmode
    
    set_plot, mode.device
    !sgraph.d1 = !d
    
    opt = {set_character_size: mode.char}
    device, _extra = opt
    
    !p.thick = mode.thick
    !p.charthick = mode.thick
    !x.thick = mode.thick
    !y.thick = mode.thick
    !z.thick = mode.thick
    !p.color = sgcolor('black')
    !p.background = sgcolor('white')
    
    erase, color = sgcolor('white')
end