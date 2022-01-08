;+
; Generate pics for given time, on Polar and FAST footpoint tracks
; mapped to the ILat-MLT plane.
;   It's possible to extend the function to plot footpoint of other s/c.
;-
pro footpoint_movie, ets, dt = dt, ps = ps, jpg = jpg
    compile_opt idl2
    
    ; constants.
    xsize = 800 & ysize = 400
    poss = [[0.05,0.10,0.45,0.90],[0.55,0.10,0.95,0.90]]
    
    ; polar:b,*, fast:r,+.
    satnames = ['polar (*,b)','fast (+,r)']
    satcolors = [sgcolor('blue'),sgcolor('red')]
    satsyms = [2,1]     ; [*,+].
    
    if n_elements(ets) ne 2 then message, 'wrong epoch input ...'
    if n_elements(dt) eq 0 then dt = 10  ; 10 min.
    
    ; read satellite orbit.
    vars = ['Epoch','mlt','ilat','dis']
    pos_po = sread_polar_pos(t = ets, vars = vars, dt = dt)
    pos_fa = sread_fast_pos(t = ets, vars = vars, dt = dt)
    pos = [pos_po,pos_fa]   ; {epoch,ilat,mlt,dis}.
    nsat = n_elements(pos)
    nrec = n_elements(pos[0].epoch)
    
    ; rootdir module.
    if n_elements(rootdir) eq 0 then begin
        case susrhost() of
            'Sheng@Xps': rootdir = sdiskdir('Research')
            'sheng@XpsMintv': rootdir = sdiskdir('Research')
            'sheng@dhcp-131.spa.umn.edu': rootdir = sdiskdir('Tian')
            else: message, 'cannot recognize the usr@host ...'
        endcase
    endif
    dataptn = rootdir+'/sdata/fpmovie/yyyy/MM/fpmovie_yyyy_MMdd_hhmm.'
    tmp = 'png'
    if keyword_set(ps) then tmp = 'ps' & if keyword_set(jpg) then tmp = 'jpg'
    dataptn += tmp
    
    ; prepare parameters.
    insouth = bytarr(nsat)
    innorth = bytarr(nsat)
    for j = 0, nsat-1 do begin
        tmp = pos[j].ilat[0] ge 0    ; 1 for N-hem, 0 for S-hem.
        innorth[j] = tmp eq 1 & insouth[j] = tmp ne 1
    endfor
    
    bufidxs = lonarr(nsat,2)    ; start from [0,0].
    bufidxn = lonarr(nsat,2)
    
    for i = 1, nrec-1 do begin
        ; update position buffer.
        for j = 0, nsat-1 do begin
            tmp = pos[j].ilat[i] ge 0   ; 1 for N-hem, 0 for S-hem.
            flip = tmp xor innorth[j]   ; 1 for flip hem.
            innorth[j] = tmp eq 1 & insouth[j] = tmp ne 1
            if flip then begin  ; flip hemisphere.
                if innorth[j] then bufidxn[j,*] = [i-1,i] $
                else bufidxs[j,*] = [i-1,i]
            endif else $        ; remain in same hemisphere.
                if innorth[j] then bufidxn[j,1] += 1 else bufidxs[j,1] += 1
        endfor
        tet = pos[0].epoch[i]
        
        ; draw window, set true color.
        sgwindow, 0, xsize = xsize, ysize = ysize, /nowindow
        !p.font = 0 & !sgraph.zmode.char *= 0.9
        if keyword_set(ps) then sgpsopen, sptn2fn(dataptn, tet) $
        else sgzopen, sptn2fn(dataptn, tet)
        lineskip = 1.5d*!d.y_ch_size/ysize
        sgtruecolor
        
        ; add label.
        xyouts, 0.5, 0.95, sfmepoch(tet,'yyyy-MM-dd hh:mm')+' UT', $
            /normal, align = 0.5, color = sgcolor('black')
        for j = 0, nsat-1 do begin
            tmp = satnames[j]+': ('+snum2str(pos[j].ilat[i])+', '+$
                string(pos[j].mlt[i],format='(F4.1)')+', '+$
                string(pos[j].dis[i],format='(F3.1)')+')'
            xyouts, 0.5, 0.95-(j+1)*lineskip, tmp, /normal, align = 0.5, $
                color = sgcolor('black')
        endfor
        
        ; draw coord and footpoint.
        sgset_map, pos = poss[*,0], color = sgcolor('black'), $
            flip = 0, south = 0, xtickpos = 47, ytickpos = 45, $
            yticknudge = [0.1,-1.3], yalign = 0.25, ytickv = [50,60,70,80]
        for j = 0, nsat-1 do begin
            idx0 = bufidxn[j,0] & idx1 = bufidxn[j,1]
            if idx0 eq idx1 then continue
            plots, pos[j].mlt[idx0:idx1]*15, pos[j].ilat[idx0:idx1], $
                color = satcolors[j]
            plots, pos[j].mlt[idx0:idx1]*15, pos[j].ilat[idx0:idx1], $
                color = satcolors[j], psym = satsyms[j], symsize = 0.8
        endfor
        ; south.
        sgset_map, pos = poss[*,1], color = sgcolor('black'), $
            flip = 1, south = 1, xtickpos =-47, ytickpos =-45, $
            yticknudge = [0.1,-1.3], yalign = 0.25, ytickv = -[50,60,70,80]
        for j = 0, nsat-1 do begin
            idx0 = bufidxs[j,0] & idx1 = bufidxs[j,1]
            if idx0 eq idx1 then continue
            plots, pos[j].mlt[idx0:idx1]*15, pos[j].ilat[idx0:idx1], $
                color = satcolors[j]
            plots, pos[j].mlt[idx0:idx1]*15, pos[j].ilat[idx0:idx1], $
                color = satcolors[j], psym = satsyms[j], symsize = 0.8
        endfor
        ; save and finish.
        print, 'saving '+sptn2fn(dataptn, tet)
        if keyword_set(ps) then sgpsclose else sgzclose
    endfor
end

;footpoint_movie, ['2000-12-28','2006-01-01']
;footpoint_movie, ['1999-10-06','1999-10-07']
footpoint_movie, ['2005-12-31','2008-06-15']
end