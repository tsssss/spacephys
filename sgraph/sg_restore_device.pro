;+
; restore the current device with the given device structure.
; reset [xy]_size, [xy]_vsize, [xy]_ch_size.
;-
pro sg_restore_device, d0

    ; device names.
    devw = (!version.os_family eq 'unix')? 'x': 'win'
    devp = 'ps' & devz = 'z'
    
    px2cm = 0.001d  ; 1/!d.x_px_cm.
    
    dev0 = strlowcase(d0.name)
    if dev0 eq devp then $  ; no way to reset inch, encapsulate.
        opt = {set_character_size:[d0.x_ch_size,d0.y_ch_size], $
        xsize:d0.x_size*px2cm, ysize:d0.y_size*px2cm, inch:0}
    if dev0 eq devz then $
        opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size], $
        set_resolution: [d0.x_size,d0.y_size]}
    if dev0 eq devw then $
        opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size]}
    set_plot, dev0
    device, _extra = opt

    device, decomposed = d0.decompose
    if d0.decompose eq 0 then tvlct, d0.r0, d0.g0, d0.b0
end