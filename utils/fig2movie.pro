;+
; Make a movie from a list of figures.
; To replace spic2movie.
;
; movie_file. Required. The movie file name.
; fig_files=. Required. The filenames of the figures, or the directory of the figures.
; fig_extension=. The figure extension. By default is all extensions.
; fps=. The frame per sec, by default is 10.
;-

pro fig2movie, movie_file, fps=fps, $
    fig_files=fig_files, fig_extension=fig_extension, errmsg=errmsg

    compile_opt idl2

    if n_elements(movie_file) eq 0 then begin
        errmsg = 'No input movie_file ...'
        return
    endif
    if n_elements(fig_files) eq 1 then begin
        path = fig_files[0]
        if n_elements(fig_extension) eq 0 then begin
            msg = '*'
        endif else msg = '*'+fig_extension[0]
        files = file_search(join_path([path,msg]))
    endif else files = fig_files

    nfile = n_elements(files)
    if nfile eq 0 then begin
        errmsg = 'No input fig_files or path ...'
        return
    endif

    flags = fltarr(nfile)
    foreach file, files, file_id do begin
        flags[file_id] = file_test(file)
    endforeach
    index = where(flags eq 1, nfile)
    if nfile eq 0 then begin
        errmsg = 'No input fig_files or path ...'
        return
    endif else begin
        files = files[index]
    endelse



    ; get colors, image size.
    if n_elements(fig_extension) eq 0 then fig_extension = fgetext(files[0])
    case fig_extension of
        'gif': begin
            read_gif, files[0], img, r, g, b
            sz = size(img, /dimensions)
            tvlct, r, g, b
            end
        'png': begin
            read_png, files[0], img
            dims = size(img, /dimensions)
            sz = dims[where(dims ne 3)]
            end
    endcase

    ; create video.
    if n_elements(fps) eq 0 then fps = 10        ; frame per sec.
    ovid = idlffvideowrite(movie_file)
    vidstream = ovid.AddVideoStream(sz[0],sz[1],fps)
    timg = bytarr(3,sz[0],sz[1])

    foreach file, files do begin
        case fig_extension of
            'gif': begin
                read_gif, file, img
                timg[0,*,*] = r[img]
                timg[1,*,*] = g[img]
                timg[2,*,*] = b[img]
                end
            'png': read_png, file, timg
        endcase
        time = ovid.put(vidstream, timg)
    endforeach

    ovid.Cleanup

end
