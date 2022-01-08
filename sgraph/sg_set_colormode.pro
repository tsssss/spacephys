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

pro sg_set_colormode, dev, decomposed = decomposed, depth = depth

    on_error, 2

    ; get the original device and target device.
    dev0 = !d.name
    dev1 = (n_elements(dev) eq 0)? dev0: dev

    dev0 = strlowcase(dev0)
    dev1 = strlowcase(dev1)

    ; set to target device.
    if dev1 ne dev0 then set_plot, dev1

    if n_elements(decomposed) eq 0 then decomposed = 0
    if n_elements(depth) eq 0 then depth = decomposed? 8: 24

    case dev1 of
        'ps': device, decomposed = decomposed, /color, bits_per_pixel = depth
        'z': device, decomposed = decomposed, set_pixel_depth = depth
        'x': device, decomposed = decomposed
        'win': device, decomposed = decomposed
        else: message, 'unknown device: '+dev1+' ...'
    endcase

    if dev1 ne dev0 then set_plot, dev0

end
