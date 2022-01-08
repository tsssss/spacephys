pro play_thm_asi

device, decomposed = 0
loadct, 1
thm_init
rootdir = spreproot('themis')


show_time = '2013-04-14/07:00:00'
thm_asi_create_mosaic, show_time, /thumb, /no_color


; exclude = ['chbg','fykn','kian','kuuj','mcgr','yknf']

;thm_asi_merge_mosaic, show_time, /verbose, $
;    exclude = ['pgeo','yknf'], /no_color, /gif_out

;thm_asi_merge_mosaic, show_time, /verbose, $
;    show = show, exclude = exclude, /no_color, projection = 'Orthographic', $
;    /gif_out, rotation = 0, central_lat = 70, central_lon = 330, $
;    xsize = 1000, ysize = 500, /mag
;    
;show_time = '2013-06-01/02:01:00'
;thm_asi_merge_mosaic, show_time, /verbose, $
;    exclude = exclude, show = show, $
;    /no_color, /gif_out;, xsize = 321, ysize = 321, /mag

end