
; img, in, [npx,npx], input image.
; lat, in, [npx+1,npx+1], input lat at pixel corner.
; lon, in, [npx+1,npx+1], input lon at pixel corner.
; mask, in, [npx,npx], input mask at pixel.

function thm_asi_proj, img, lat, lon, mask

    fillval = !values.d_nan

    ; **** no standard lat/lon grid.
    minlat = min(lat, max = maxlat, /nan)
    minlon = min(lon, max = maxlon, /nan)
    npx = (size(img, /dimensions))[0]
    clat = npx/(maxlat-minlat)
    clon = npx/(maxlon-minlon)
    img1 = fltarr(npx,npx)
    cnt = intarr(npx,npx)
    
    ; bin the lat/lon corner according to uniform lat/lon grid.
    ilat = floor((lat-minlat)*clat+0.5)
    ilon = floor((lon-minlon)*clon+0.5)
    
    for i = 0, npx-1 do begin
        for j = 0, npx-1 do begin
            ib = max(ilat[i:i+1,j:j+1], min = ia)
            jb = max(ilon[i:i+1,j:j+1], min = ja)
            if max([ia,ib,ja,jb]) gt npx then continue   ; nan.
            ib = ib < (npx-1) > (ia+1)
            jb = jb < (npx-1) > (ja+1)
            img1[ia:ib,ja:jb] += img[i,j]
            cnt[ia:ib,ja:jb] += 1
        endfor
    endfor
    
    img1 /= cnt
;    loadct, 1
    tv, bytscl(img,/nan), 0
    tv, rotate(bytscl(img1,/nan),1),1
    stop

end

device, decomposed = 0
loadct, 1
; **** thumbnail data.
; thm_asi_create_mosaic.
show_time = '2008-03-09/06:50:00'
show = ['kapu']
thm_init, local_data_dir = spreproot(), /reset
thm_asi_create_mosaic, show_time, show = show, /thumb, /no_color
; read_thm_asi.
et = stoepoch('2008-03-09 06:50:00')
sites = 'kapu'
type = 'asf'
rootdir = spreproot()
asi = read_thm_asi(et, sites, type = type, $
    rootdir = rootdir, /debug)

end