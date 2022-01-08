;+
; if color mode is true color, using it.
; if color mode is index color, convert indexed colors to true color.
;-
pro test_sgtv

    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

    ; input.
    img0 = dist(500,300)
    ; true = true
    
    device, get_decomposed = cmode0
    if cmode0 eq 1 then cmode = 1
    
    imgsz = size(img0,/dimensions)
    case n_elements(imgsz) of
        2: cmode = 1
        3: cmode = 0
        else: message, 'invalid image dimension ...'
    endcase
    
    ; cmode is suggested by decompose, but determined by image size.
    if cmode eq 0 then img = img0 else begin
        img = bytarr([3,imgsz])
        true = 1
        ncolor = n_elements(r_curr)
        if ncolor eq 0 then loadct, 0, /silent
        ncolor = n_elements(r_curr)
        img[0,*,*] = r_curr[img0 mod ncolor]
        img[1,*,*] = g_curr[img0 mod ncolor]
        img[2,*,*] = b_curr[img0 mod ncolor]
    endelse
    
    ; determine true.
    imgsz = size(imgsz,/dimensions)
    if n_elements(true) eq 0 then true = where(imgsz eq 3)+1
    
    tv, img, true = true
    
    device, decomposed = cmode0

end