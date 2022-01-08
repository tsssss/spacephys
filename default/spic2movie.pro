;+
; dir: the directory containing the images (string) or the images (strarr[n]).
; vfn: movie filename.
; fps: frame per sec, default is 10.
;-
pro spic2movie, vfn, fps, plot_files=plot_files, plot_dir=dir, plot_ext=ext
    compile_opt idl2

    
    sep = path_sep()
    if n_elements(plot_files) eq 0 then begin
        fns = file_search(dir)+sep+'*'
        if n_elements(ext) ne 0 then fns += '.'+ext
        fns = file_search(fns)
    endif else begin
        fns = plot_files
    endelse
    nfn = n_elements(fns)
    
    if n_elements(ext) eq 0 then ext = fgetext(fns[0])
    
    ; get colors, image size.
    case ext of
        'gif': begin
            read_gif, fns[0], img, r, g, b
            sz = size(img, /dimensions)
            tvlct, r, g, b
            end
        'png': begin
            read_png, fns[0], img
            dims = size(img, /dimensions)
            sz = dims[where(dims ne 3)]
            end
    endcase
    
    ; create video.
    if n_elements(fps) eq 0 then fps = 10        ; frame per sec.
    ovid = idlffvideowrite(vfn)
    vidstream = ovid.AddVideoStream(sz[0],sz[1],fps)
    timg = bytarr(3,sz[0],sz[1])
    
    for i = 0, nfn-1 do begin
        if fns[i] eq '' then continue
        case ext of
            'gif': begin
                read_gif, fns[i], img
                timg[0,*,*] = r[img]
                timg[1,*,*] = g[img]
                timg[2,*,*] = b[img]
                end
            'png': read_png, fns[i], timg
        endcase
        time = ovid.put(vidstream, timg)
    endfor
    
    ovid.Cleanup

end

dir = shomedir()+'/Google Drive/works/works/aurora compare/2008_0202_clip1'
vfn = shomedir()+'/thm_asf_2008_0202_clip1.mp4'

dir = '/Users/Sheng/psbl_de_32hz/thg_asi_3site'
vfn = '/Users/Sheng/psbl_de_32hz/thg_2013_0607_3site.mp4'

;dir = '/Users/Sheng/psbl_de_32hz/thg_asi_quarter'
;vfn = '/Users/Sheng/psbl_de_32hz/thg_2013_0607_pina.mp4'

dir = '/Users/Sheng/psbl_de_32hz/thg_asi_insitu_t89'
vfn = '/Users/Sheng/psbl_de_32hz/thg_2013_0607_insitu.mp4'


dir = '/Users/Sheng/psbl_de_32hz/asi_iono2tail_v2'
vfn = '/Users/Sheng/psbl_de_32hz/thg_2013_0607_insitu_v2.mp4'

dir = '/Users/Sheng/psbl_de_32hz/asi_iono2tail_v3'
vfn = '/Users/Sheng/psbl_de_32hz/thg_2013_0607_insitu_v3.mp4'


ext = 'png'
spic2movie, dir, vfn, ext

end