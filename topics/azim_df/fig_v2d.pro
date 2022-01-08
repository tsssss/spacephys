;+
; For the v_triad determined for all events:
;   Plot the direction and magnitude.
;-


test = 1
    magnify = (keyword_set(test))? 2: 1

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)
    timing_key = '_obs_time'


    fig_xsize = 5
    fig_ysize = 3
    file = join_path([project.plot_dir,'fig_v2d.pdf'])
    if keyword_set(test) then file = 0
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch, magnify=magnify
    margins = [7,4,2,2]
    poss = sgcalcpos(1, margins=margins, xchsz=xchsz, ychsz=ychsz)

    label_size = constant('label_size')
    xticklen_chsz = -0.3
    yticklen_chsz = -0.5

;---Panel a. Direction.
    tpos = poss
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xrange = [-1,1]*9
    xstep = 3
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = 3
    xtitle = 'MLT (hr)'

    yrange = [-1,1]*210
    yrange = [-180,270]+[-30,0]
    ystep = 90
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 3
    ytitle = tex2str('phi')+' (deg)'

;---Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtitle=xtitle, $
        ystyle=5, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase


;---4 schemes.
    schemes = orderedhash($
        'eastward', dictionary($
            'label', 'Eastward', $
            'pre_midn', [-12.,0], $
            'pre_phi', [-270.,-90], $
            'post_midn', [-0.,12], $
            'post_phi', [-90.,90], $
            'color', sgcolor('blue')), $
         'westward', dictionary($
            'label', 'Westward', $
            'pre_midn', [0.,12], $
            'pre_phi', [90.,270], $
            'post_midn', [-12,0], $
            'post_phi', [-90.,90], $
            'color', sgcolor('green')), $
        'earthward', dictionary($
            'label', 'Radial.In', $
            'pre_midn', [-12.,0], $
            'pre_phi', [-180.,0], $
            'post_midn', [0.,12], $
            'post_phi', [0.,180], $
            'color', sgcolor('red')), $
        'outward', dictionary($
            'label', 'Radial.Out', $
            'pre_midn', [-12.,0], $
            'pre_phi', [0.,180], $
            'post_midn', [0.,12], $
            'post_phi', [-180.,0], $
            'post_phi', [180.,360], $
            'color', sgcolor('purple')))


    linestyle = 2
    thick = (keyword_set(test))? 2: 8
    foreach scheme, schemes do begin
        color = scheme.color
        foreach key, ['pre','post'] do begin
            txs = scheme[key+'_midn']
            tys = scheme[key+'_phi']
            oplot, txs,tys, linestyle=linestyle, color=color
        endforeach
    endforeach


;---Plot all data.
    deg = constant('deg')
    foreach event, events do begin
        direction = (strsplit(event.region,'%',/extract))[1]
        scheme = schemes[direction]
        color = scheme.color

        v_phi = list()
        v_mlt = list()

        triad_list = event.triad_list
        foreach triad, triad_list do begin
            center_r_sm = triad.center_r_sm

            the_mlt = pseudo_mlt(center_r_sm)
            v_mlt.add, the_mlt

            vhat = triad['vhat'+timing_key]
            the_phi = atan(vhat[1],vhat[0])*deg

            mlt_range = (the_mlt lt 0)? scheme.pre_midn: scheme.post_midn
            phi_range = (the_mlt lt 0)? scheme.pre_phi: scheme.post_phi
            target_phi = interpol(phi_range,mlt_range,the_mlt)
            if the_phi-target_phi lt 180 then the_phi += 360
            if the_phi-target_phi gt 180 then the_phi -= 360
            v_phi.add, the_phi
        endforeach

        v_phi = v_phi.toarray()
        v_mlt = v_mlt.toarray()
        index = sort(v_mlt)
        v_phi = v_phi[index]
        v_mlt = v_mlt[index]
        grid_color = sgcolor('silver')
        plots, v_mlt, v_phi, color=grid_color
        plots, v_mlt, v_phi, psym=1, symsize=label_size*0.5, color=color
    endforeach

;---Add labels.
    full_ychsz = constant('full_ychsz')
    tx_label = tpos[0]+xchsz*1
    ty_label = tpos[3]-ychsz*full_ychsz
    ii = 0
    foreach scheme, schemes do begin
        color = scheme.color
        tx = tx_label
        ty = ty_label-ychsz*full_ychsz*ii & ii += 1
        xyouts, tx,ty,/normal, scheme.label, color=color, charsize=label_size
    endforeach



;---Add axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase

    grid_color = sgcolor('silver')
    constants = [-180,-90,0,90,180]
    foreach val, constants do plots, xrange, val+[0,0], linestyle=1, color=grid_color
    constants = [-6,0,6]
    foreach val, constants do plots, val+[0,0], yrange, linestyle=1, color=grid_color

    tx = tpos[0]-xchsz*0
    ty = tpos[3]+ychsz*0.6
    xyouts, tx,ty,/normal, 'd. v!D2D!N direction vs MLT'



    if keyword_set(test) then stop
    sgclose

end
