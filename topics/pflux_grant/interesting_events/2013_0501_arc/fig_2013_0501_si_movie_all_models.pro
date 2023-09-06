;+
; Generate the movie for supplemental information.
;
; Adopted from fig_si_movie
;-

function fig_2013_0501_si_movie_all_models, movie_file, event_info=event_info


;---Load data and settings.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
test = 1

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = event_info['asi_time_range']
    model_setting = event_info['model_setting']
    internal_model = model_setting['internal_model']
    models = model_setting['models']
    asi_setting = event_info['asi_setting']
    site = (asi_setting['sites'])[0]

    if n_elements(movie_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_all_models.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [7,4,3,1]

;---ASI.
    mlt_range = [-2,0.5]
    mlat_range = [58,68]
    asi_ct = 49
    top_color = 254

    xtitle = 'MLT (h)'
    xrange = mlt_range+24
    xstep = 0.5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtickn = string(xtickv,format='(F4.1)')
    foreach tx, xtickv, ii do begin
        if tx lt 24 then continue
        xtickn[ii] = strtrim(string(xtickv[ii]-24,format='(F4.1)'),2)
    endforeach

    ytitle = 'MLat (deg)'
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 5


    zrange = [0,1e4]
    ztitle = 'ASI Count (#)'

    mlt_image_var = 'thg_asf_mlt_image_rect'
    get_data, mlt_image_var, times, mlt_images
    index = where_pro(times, '[]', time_range)
    times = times[index]
    mlt_images = mlt_images[index,*,*]

    plot_files = []
    foreach time, times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_rect',$
            'thg_asf_mlt_image_rect_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, xsize=4, ysize=2.4, magnify=magn
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, margins=margins)

        mlt_image = reform(mlt_images[time_id,*,*])
        mlt_bins = get_setting(mlt_image_var, 'mlt_bins')
        mlt_index = where_pro(mlt_bins, '[]', mlt_range)
        mlt_bins = mlt_bins[mlt_index]
        mlt_image = mlt_image[mlt_index,*]
        mlat_bins = get_setting(mlt_image_var, 'mlat_bins')
        mlat_index = where_pro(mlat_bins, '[]', mlat_range)
        mlat_bins = mlat_bins[mlat_index]
        mlt_image = mlt_image[*,mlat_index]

        tpos[3] = tpos[3]-ychsz*2.5
        cbpos = tpos
        cbpos[1] = tpos[3]+ychsz*0.2
        cbpos[3] = cbpos[1]+ychsz*0.5

        zzs = bytscl(mlt_image, min=zrange[0],max=zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos, resize=1
        sgcolorbar, findgen(top_color), horizontal=1, ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Add grid.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle='', xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xminor=1, xticklen=1, xgridstyle=1, xtickname=xtickn, $
            ystyle=1, ylog=0, yrange=yrange, ytitle='', ytickformat='(A1)', ytickv=ytickv, yticks=yticks, yminor=1, yticklen=1, ygridstyle=1, $
            position=tpos, nodata=1, noerase=1, ynozero=1, color=sgcolor('gray')

        ; Draw axes.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
            ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
            position=tpos, nodata=1, noerase=1, ynozero=1

    ;---Add footpoint.
        fmlts = get_var_data(prefix+'fmlt_'+internal_model, at=time)+24
        fmlats = get_var_data(prefix+'fmlat_'+internal_model, at=time)
        model_colors = sgcolor(['red','green','blue','black'])
        foreach external_model, models, model_index do begin
            sc_color = event_info['sc_color']
            model_color = model_colors[model_index]
            fmlt = fmlts[model_index]
            fmlat = fmlats[model_index]

            tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=0.5, color=model_color

            ; Add model name.
            tx = tpos[0]+xchsz*(1+model_index*5)
            ty = tpos[3]-ychsz*1.5
            xyouts, tx,ty,normal=1, strupcase(external_model), color=model_color
        endforeach



    ;---Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = strupcase(site)+' '+time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('black')


;    ;---Add FAC.
;        j_ver = get_var_data('thg_j_ver', at=time, limits=lim)
;        pixel_mlons = lim.mlon_grids
;        pixel_mlats = lim.mlat_grids
;        pixel_mlts = mlon2mlt(pixel_mlons, time)+24
;        index = where(pixel_mlts gt xrange[0] and pixel_mlts lt xrange[1] and $
;            pixel_mlats gt yrange[0] and pixel_mlats lt yrange[1], npixel)
;        xxs = pixel_mlts[index]
;        yys = pixel_mlats[index]
;        zzs = j_ver[index]
;        fac_zrange = [-1,1]*1e5
;        ccs = bytscl(zzs, min=fac_zrange[0], max=fac_zrange[1])
;        dc = 80
;        ccs[where(zzs ge 0)] = 128+dc
;        ccs[where(zzs lt 0)] = 128-dc
;        ct = 70
;        tmp = smkarthm(0,2*!dpi,30,'n')
;        txs = cos(tmp)
;        tys = sin(tmp)
;        usersym, txs, tys
;        
;        symszs = (abs(zzs/fac_zrange[1]))^0.25*0.5
;        for ii=0,npixel-1 do begin
;            cc = sgcolor(ccs[ii], ct=ct)
;            symsz = symszs[ii]
;            ;symsz = 0.5
;            ;cc = (zzs[ii] ge 0)? sgcolor('blue'): sgcolor('red')
;            plots, xxs[ii], yys[ii], color=cc, psym=8, symsize=symsz
;        endfor


        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files
    return, movie_file

end


print, fig_2013_0501_si_movie_all_models(event_info=event_info)
end