
function sgcolor_index2true, img, ct0
    
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

    case n_elements(ct0) of
        ; color table id.
        1: loadct, ct0, rgb_table = rgb
        ; use current table.
        0: rgb = [[r_curr],[g_curr],[b_curr]]
        ; ct0 is rgb.
        else: rgb = ct0
    endcase

    ncolor = n_elements(rgb)/3
    imgr = (rgb[*,0])[img mod ncolor]
    imgg = (rgb[*,1])[img mod ncolor]
    imgb = (rgb[*,2])[img mod ncolor]
    
    return, [[[imgr]],[[imgg]],[[imgb]]]
end

pro test_sg_colormap

xsz = 600
ysz = 400
ct = 30

; index color mode.
device, decomposed = 0
loadct, ct
img0 = dist(xsz,ysz)
window, 0, xsize = xsz, ysize = ysz
tv, img0

; true color mode.
device, decomposed = 1
window, 1, xsize = xsz, ysize = ysz
img1 = sgcolor_index2true(img0)
tv, img1, true = 3


end