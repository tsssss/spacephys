
pro sgzclose
    fn = !sgraph.zmode.sgid
    device, get_pixel_depth = true
    true = (true eq 8)? 0: 1
    case strlowcase((reverse(strsplit(fn,'.',/extract)))[0]) of
        'png': begin
            if true then write_png, fn, tvrd(/true) else begin
                tvlct, r, g, b, /get
                write_png, fn, tvrd(), r, g, b
            endelse
        end
        'jpg': begin tvlct, r, g, b, /get
            write_jpeg, fn, tvrd(true=3), true=3 & end
        'jpeg': begin tvlct, r, g, b, /get
            write_jpeg, fn, tvrd(true=3), true=3 & end
    endcase
    
    ; restore !d for 'z'.
    d0 = !sgraph.d1
    opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size], $
        set_resolution: [d0.x_size,d0.y_size]}
    device, _extra = opt
    
    ; restore color.
    if d0.n_colors eq 256 then sgindexcolor else sgtruecolor
    
    ; restore, !p, ![xyz].
    !p = !sgraph.p0
    !x = !sgraph.x0
    !y = !sgraph.y0
    !z = !sgraph.z0
    
    d0 = !sgraph.d0
    set_plot, d0.name
    opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size]}
    device, _extra = opt
    if d0.n_colors eq 256 then sgindexcolor else sgtruecolor
end