
pro fig_fpoint, ps = ps, png = png

    fn = '~/fpoint.'
    if keyword_set(ps) then fn += 'eps'
    if keyword_set(png) then fn += 'png'

    sgwindow, 0, xsize = 3.5, ysize = 2, /inch
    if keyword_set(ps) then sgpsopen, fn
    if keyword_set(png) then sgzopen, fn
    sgtruecolor
    
    pos = [0.1,0.1,0.9,0.9]
    sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
        ytickv = [50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
        xtickpos = 47
    
    vars = ['Epoch','mlt','ilat','dis']
    
    ; polar:b,*, fast:r,+.
    satnames = ['Polar (blue)','Fast (red)']
    satcolors = [sgcolor('blue'),sgcolor('red')]
    satsyms = [2,1]     ; [*,+].
    
    ; polar.
    ets = stoepoch(['1998-09-25/01:00','1998-09-25/06:10']) & dt = 30 & j = 0
    pos_po = sread_polar_pos(t = ets, vars = vars, dt = dt)
    mlt = pos_po.mlt*15 & ilat = pos_po.ilat & et = pos_po.epoch
    plots, mlt[0:-2], ilat[0:-2], color = satcolors[j]
    plots, mlt[0:-2], ilat[0:-2], color = satcolors[j], $
        psym = satsyms[j], symsize = 0.4
    arrow, mlt[-2], ilat[-2], mlt[-1], ilat[-1], /data, $
        color = satcolors[j], /solid, thick = 2
    i0 = -1 & i1 = -2 & chsz = 0.7
    xyouts, mlt[i0], ilat[i0], sfmepoch(et[i0],'hh:mm'), $
        alignment = 1.1, color = satcolors[j], charsize = chsz
    xyouts, mlt[i1], ilat[i1], sfmepoch(et[i1],'hh:mm'), $
        alignment = -0.3, color = satcolors[j], charsize = chsz
    xyouts, 0.8, 0.30, /norm, satnames[j], color = satcolors[j], align = 1

    ; fast.
    ets = stoepoch(['1998-09-25/04:20','1998-09-25/04:45']) & dt = 5 & j = 1
    pos_fa = sread_fast_pos(t = ets, vars = vars, dt = dt)
    mlt = pos_fa.mlt*15 & ilat = pos_fa.ilat & et = pos_fa.epoch
    plots, mlt[0:-2], ilat[0:-2], color = satcolors[j]
    plots, mlt[0:-2], ilat[0:-2], color = satcolors[j], $
        psym = satsyms[j], symsize = 0.4
    arrow, mlt[-2], ilat[-2], mlt[-1], ilat[-1], /data, $
        color = satcolors[j], /solid, thick = 2
    i0 = 2 & i1 = 3 & chsz = 0.7
    xyouts, mlt[i0], ilat[i0], sfmepoch(et[i0],'hh:mm'), $
        alignment = -0.3, color = satcolors[j], charsize = chsz
    xyouts, mlt[i1], ilat[i1], sfmepoch(et[i1],'hh:mm'), $
        alignment = 1.3, color = satcolors[j], charsize = chsz
    xyouts, 0.8, 0.24, /norm, satnames[j], color = satcolors[j], align = 1
        
    if keyword_set(ps) then begin sgpsclose & wdelete, 0 & endif
    if keyword_set(png) then begin sgzclose & wdelete, 0 & endif
end

fig_fpoint, /ps
fig_fpoint
end