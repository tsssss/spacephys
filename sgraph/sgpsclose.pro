
pro sgpsclose, pdf = pdf
    device, /close
    
    ; restore !d for 'ps'.
    d0 = !sgraph.d1
    opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size], $
        xsize: d0.x_size, ysize: d0.y_size}
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
    
    if keyword_set(pdf) then begin
        psfn = !sgraph.pmode.sgid
        sps2pdf, psfn, /rm
    endif
end