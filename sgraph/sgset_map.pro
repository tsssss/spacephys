;+
; 0 degree is in downward direction, longitude increases counterclockwise.
; 
; [xy]major, in, int, opt. Set # of major ticks. Default x:4, y:5.
; [xy]minor, in, int, opt. Set # of minor ticks. Default x:1, y:0.
; [xy]tickv, in, intarr[n], opt. Set major tick values. Overwrite [xy]major.
; [xy]tickname, in, strarr[n], opt. Set major tick names.
; [xy]tickpos, in, int/intarr[n]/intarr[2,n]. Set position to pose tick name.
;   For xtick, int set the latitude, intarr[n] sets latitude for each tick,
;   intarr[n,2] sets [lon,lat] position for each tick. ytick is similar.
; [xy]ticknudge, in, intarr[2], opt. Nudge ytick position in [xy]charsize unix.
; [major,minor]grid, in, 0-6, opt. Set linestyle for major/minor grid. Major
;   grid has label while minor grid does not. Default is solid line for both.
; no[xy]tick, in, opt. Set to suppress ticknames.
; mlon, in, boolean, opt. Set ytick name to be longitude.
;-


pro sgset_map, $
    xrange = xrange, xmajor = xmajor, xminor = xminor, xtickv = xtickv, $
    xtickname = xtickname, xtickpos = xtickpos, xticknudge = xticknudge, $
    yrange = yrange, ymajor = ymajor, yminor = yminor, ytickv = ytickv, $
    ytickname = ytickname, ytickpos = ytickpos, yticknudge = yticknudge, $
    majorgrid = majorgrid, minorgrid = minorgrid, $
    noxtick = noxtick, noytick = noytick, $
    color = color, position = position, limit = limit, $
    flip = flip, south = south, mlon = mlon, yalign = yalign, $
    charsize = charsize
    
    ; flip = 1, means s-hem, 06 and 18 are flipped.
    if keyword_set(flip) then reverse = 1 else reverse = 0
    
    on_error, 0
    
    if n_elements(position) eq 0 then position = [0.1,0.1,0.9,0.9]
    if n_elements(xrange) eq 0 then xrange = [0,360]
    if n_elements(yrange) eq 0 then yrange = [50,90]
    if keyword_set(south) then begin
        yrange = -yrange
        reverse += 2
    endif
    limit = [yrange[0],xrange[0],yrange[1],xrange[1]]
    p0lat = limit[2] & p0lon = 0 & rot = 0
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    chsz0 = !p.charsize
    if keyword_set(charsize) then !p.charsize = charsize
    
    map_set, name = 'AzimuthalEquidistant', p0lat, p0lon, rot, $
        /noborder, /noerase, /isotropic, position = position, limit = limit, $
        reverse = reverse
    
    ; x:lon, y:lat.
    if n_elements(xtickv) eq 0 then begin
        if n_elements(xmajor) eq 0 then $
            xmajor = floor(abs(xrange[1]-xrange[0])/90)+1
        xmajor = xmajor>2
        if xmajor eq 0 then xtickv = xrange else $
            xtickv = xrange[0]+findgen(xmajor)*(xrange[1]-xrange[0])/(xmajor-1d)
    endif else xmajor = n_elements(xtickv)
    if n_elements(xminor) eq 0 then xminor = 1
    if n_elements(xtickname) eq 0 then begin
        if ~keyword_set(mlon) then begin
            xtickname = string((xtickv mod 360)/15,format='(I02)')
        endif else xtickname = string(xtickv mod 360,format='(I3)')
    endif
    
    if n_elements(ytickv) eq 0 then begin
        if n_elements(ymajor) eq 0 then $
            ymajor = floor(abs(yrange[1]-yrange[0])/10)+1
        ymajor = ymajor>2
        if ymajor eq 0 then ytickv = yrange else $
            ytickv = yrange[0]+findgen(ymajor)*(yrange[1]-yrange[0])/(ymajor-1d)
    endif else ymajor = n_elements(ytickv)
    if n_elements(yminor) eq 0 then yminor = 0
    if n_elements(ytickname) eq 0 then ytickname = snum2str(ytickv,0)
    
    if n_elements(majorgrid) eq 0 then majorgrid = 0
    if n_elements(minorgrid) eq 0 then minorgrid = 0
    
    ; draw minor grid.
;    if xminor ne 0 then $
;        map_grid, lats = yrange, glinestyle = minorgrid, color = color, $
;            lons = (xtickv[0:-2]##(fltarr(xminor)+1)+$
;                (findgen(xminor)+1)#(xtickv[1:*]-xtickv[0:-2])/(xminor+1d))[*]
    if xminor ne 0 then begin
        lons = (xtickv[0:-2]##(fltarr(xminor)+1)+$
            (findgen(xminor)+1)#(xtickv[1:*]-xtickv[0:-2])/(xminor+1d))[*]
        for i = 0, n_elements(lons)-1 do $
            plots, lons[[i,i]], yrange, linestyle = minorgrid, color = color
    endif
    if yminor ne 0 then $
        map_grid, lons = xrange, glinestyle = minorgrid, color = color, $
            lats = reform((ytickv[0:-2]##(fltarr(yminor)+1)+$
                (findgen(yminor)+1)#(ytickv[1:*]-ytickv[0:-2])/(yminor+1d))[*])
    
    ; draw major grid.
    map_grid, lons = xtickv, lats = yrange, glinestyle = majorgrid, color = color
    map_grid, lons = xrange, lats = ytickv, glinestyle = majorgrid, color = color
        
    ; draw ticks.
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    case n_elements(xtickpos) of    ; set latitude mostly.
        0: xtkpos = transpose([[xtickv],[dblarr(xmajor)+yrange[0]]])
        1: xtkpos = transpose([[xtickv],[dblarr(xmajor)+xtickpos]])
        xmajor-1: xtkpos = transpose([[xtickv],[(xtickpos)[*],xtickpos[0]]])
        xmajor: xtkpos = transpose([[xtickv],[(xtickpos)[*]]])
        else: ; do nothing.
    endcase
    if n_elements(yalign) eq 0 then yalign = 0.5
    tmp = convert_coord(/data,/to_norm, xtkpos)
    tmp[1,*] -= yalign*ychsz
    xtkpos = (convert_coord(/norm,/to_data, tmp))[0:1,*]
    if n_elements(xticknudge) ge 2 then begin
        if n_elements(xticknudge) eq 2 then $
            xticknudge = xticknudge[*]#(dblarr(xmajor)+1)
        if keyword_set(flip) then xticknudge[0,*] *= -1
        tmp = convert_coord(/data,/to_norm, xtkpos)
        tmp[0,*] += xchsz*xticknudge[0,*] & tmp[1,*] += ychsz*xticknudge[1,*]
        xtkpos = (convert_coord(/norm,/to_data, tmp))[0:1,*]
    endif
    if ~keyword_set(noxtick) then begin
        tmp = convert_coord(xtkpos[0,*], xtkpos[1,*], /data, /to_normal)
        xyouts, tmp[0,*], tmp[1,*]+ychsz*0.2*chsz0, xtickname, alignment=0.5, /normal, color=color
    endif
    
    case n_elements(ytickpos) of
        0: ytkpos = transpose([[dblarr(ymajor)+xrange[0]+45],[ytickv]])
        1: ytkpos = transpose([[dblarr(ymajor)+ytickpos],[ytickv]])
        ymajor: ytkpos = transpose([[ytickv],[reform(ytickpos)]])
    endcase
    if n_elements(yticknudge) ge 2 then begin
        if n_elements(yticknudge) eq 2 then $
            yticknudge = yticknudge[*]#(dblarr(ymajor)+1)
            if keyword_set(flip) then yticknudge[0,*] *= -1
        tmp = convert_coord(/data,/to_norm, ytkpos)
        tmp[0,*] += xchsz*yticknudge[0,*] & tmp[1,*] += ychsz*yticknudge[1,*]
        ytkpos = (convert_coord(/norm,/to_data, tmp))[0:1,*]
    endif
    if ~keyword_set(noytick) then begin
        tmp = convert_coord(ytkpos[0,*], ytkpos[1,*], /data, /to_normal)
        xyouts, tmp[0,*], tmp[1,*], ytickname, alignment=0.5, /normal, color=color
    endif
    
    !p.charsize = chsz0
end

chsz0 = !p.charsize
!p.charsize = 1.5
window, 1, xsize = 500, ysize = 500
;map_set, 90, 0, 0, name = 'AzimuthalEquidistant', /noborder, $
;    limit = [50,0,90,360], /isotropic, $
;    latdel = 10, londel = 45, glinestyle = 1
;plots, [0,90,180,270], [90,80,70,60,50], color = sgcolor('blue')
;wait, 1
;sgset_map, ytickpos = 135, xtickpos = 48, xrange = [90,270]
;sgset_map, ytickpos = 135, xtickpos = [48,49,48], xrange = [90,270]
;sgset_map, ytickpos = 135, xtickpos = [48,49,48], xrange = [90,270], yticknudge = [2.2,-0.3]
;sgset_map, xrange = [90,270], xtickpos = 48, ytickpos = 90, yticknudge = [[0,-1.2],[0,-1.2],[0,-1.5],[0,-1.7],[0,-2]]
sgset_map, xtickpos = 48, ytickpos = 45, yticknudge = [2.2,-0.4]
plots, [0,90,180,270,360], [90,80,70,60,50], color = sgcolor('red')
!p.charsize = chsz0
end
