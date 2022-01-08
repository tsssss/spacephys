;+
; Type: procedure.
; Purpose: Test true color vs index color imaging.
; Parameters: none.
; Keywords: none.
; Notes: none.
; Dependence: none.
; History:
;   2014-02-03, Sheng Tian, create.
;-

pro test_image_color
    on_error, 2

    img = cgdemodata(5)
    size = size(img, /dimensions)
    ctid = 1

    ; index color.
    device, decomposed = 0
    loadct, ctid, /silent
    window, 0, xsize = size[0], ysize = size[1]
    tv, img

    ; true color: works perfectly.
    device, decomposed = 0 & loadct, ctid, /silent
    common colors, r0, g0, b0, r1, g1, b1
    device, decomposed = 1
    img1 = [[[r1[img]]],[[g1[img]]],[[b1[img]]]]
    window, 1, xsize = size[0], ysize = size[1]
    tv, img1, true = 3
    
    ; true color 2: does not work.
    device, decomposed = 1
    img2 = r1[img]+g1[img]*256L+b1[img]*65536L
    window, 2, xsize = size[0], ysize = size[1]
    tv, img2

end
