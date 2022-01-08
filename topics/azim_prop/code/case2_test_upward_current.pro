;+
; Test to create a ewogram for the upward current.
;-

    mlat_range = [60,80]
    time_range = time_double(['2014-08-28/10:00','2014-08-28/10:50'])
    mlon_range = [-120,-0]
    time_range = time_double(['2014-08-28/08:00','2014-08-28/24:00'])
    mlon_range = [50,-100]
    mlat_binsize = 2
    mlon_binsize = 4
    mlat_bins = make_bins(mlat_range,mlat_binsize)
    mlon_bins = make_bins(mlon_range,mlon_binsize)
    nmlon_bin = n_elements(mlon_bins)
    nmlat_bin = n_elements(mlat_bins)
    new_image_size = [nmlon_bin,nmlat_bin]


    j_var = 'thg_j_ver'
    if check_if_update(j_var, time_range) then themis_read_weygand, time_range

    ; Read the orig data, glon/glat bins.
    get_data, j_var, times, j_orig
    glonbins = get_setting(j_var, 'glonbins')
    glatbins = get_setting(j_var, 'glatbins')
    glonbinsize = total(glonbins[[0,1]]*[-1,1])
    glatbinsize = total(glatbins[[0,1]]*[-1,1])
    nglonbin = n_elements(glonbins)
    nglatbin = n_elements(glatbins)
    old_image_size = [nglonbin,nglatbin]


    ; Get the mlon/mlat for glon/glat bins.
    pixel_glons = fltarr(nglonbin,nglatbin)
    pixel_glats = fltarr(nglonbin,nglatbin)
    for ii=0,nglatbin-1 do pixel_glons[*,ii] = glonbins
    for ii=0,nglonbin-1 do pixel_glats[ii,*] = glatbins

;    corner_glons = fltarr(nglonbin,nglatbin,4)
;    corner_glats = fltarr(nglonbin,nglatbin,4)
;    for ii=0,nglatbin-1 do begin
;        corner_glons[*,ii,0] = glonbins-glonbinsize*0.5
;        corner_glons[*,ii,3] = glonbins-glonbinsize*0.5
;        corner_glons[*,ii,1] = glonbins+glonbinsize*0.5
;        corner_glons[*,ii,2] = glonbins+glonbinsize*0.5
;    endfor
;    for ii=0,nglonbin-1 do begin
;        corner_glats[ii,*,0] = glatbins-glatbinsize*0.5
;        corner_glats[ii,*,1] = glatbins-glatbinsize*0.5
;        corner_glats[ii,*,2] = glatbins+glatbinsize*0.5
;        corner_glats[ii,*,3] = glatbins+glatbinsize*0.5
;    endfor

    apexfile = join_path([homedir(),'Projects','idl','spacephys','aurora','image','support','mlatlon.1997a.xdr'])
    geotoapex, pixel_glats, pixel_glons, apexfile, pixel_mlats, pixel_mlons


    ; Map to uniform mlon/mlat bins.
    mlon_bin_min = mlon_range[0]
    mlat_bin_min = mlat_range[0]
    i0_bins = round((pixel_mlons-mlon_bin_min)/mlon_binsize)
    j0_bins = round((pixel_mlats-mlat_bin_min)/mlat_binsize)

    i1_range = [0,nmlon_bin-1]
    j1_range = [0,nmlat_bin-1]

    i_bins = make_bins(i1_range, 1)
    j_bins = make_bins(j1_range, 1)
    ni_bin = nmlon_bin
    nj_bin = nmlat_bin

    index_map_from_old = list()
    index_map_to_new = list()
    for ii=0, nmlon_bin-1 do begin
        mlon_range = mlon_bins[ii]+[-1,1]*mlon_binsize*0.5
        for jj=0, nmlat_bin-1 do begin
            mlat_range = mlat_bins[jj]+[-1,1]*mlat_binsize*0.5
            index = where($
                pixel_mlons ge mlon_range[0] and $
                pixel_mlons lt mlon_range[1] and $
                pixel_mlats ge mlat_range[0] and $
                pixel_mlats lt mlat_range[1], count)
            if count eq 0 then continue
            index_map_from_old.add, index
            index_map_to_new.add, ii+jj*nmlon_bin
        endfor
    endfor

    ntime = n_elements(times)
    j_new = fltarr([ntime,new_image_size])
    for ii=0,ntime-1 do begin
        img_old = reform(j_orig[ii,*,*])
        img_new = fltarr(new_image_size)
        foreach pixel_new, index_map_to_new, pixel_id do begin
            img_new[pixel_new] = mean(img_old[index_map_from_old[pixel_id]])
        endforeach
        j_new[ii,*,*] = img_new
    endfor


;---EWOgram.
    mlat_range = [60,70]
    mlat_index = lazy_where(mlat_bins, '[]', mlat_range)
    ewo = total(-j_new[*,*,mlat_index], 3)/n_elements(mlat_index)
    ewo = fltarr(ntime,nmlon_bin)
    foreach time, times, ii do begin
        foreach mlon, mlon_bins, jj do begin
            tmp = reform(-j_new[ii,jj,mlat_index])
            index = where(tmp lt 0, count)
            if count ne 0 then tmp[index] = 0
            ewo[ii,jj] = mean(tmp)   ; total -> zrange [0.2.5e5], max -> zrange [0,2.5e5], mean -> zrange [0,1.5e5], total -> zrange [3e5]
        endforeach
    endforeach
    ewo_var = 'thg_j_ver_ewo'
    store_data, ewo_var, times, ewo, mlon_bins, limits={$
        spec: 1, $
        no_interp: 1, $
        ytitle: 'MLon (deg)', $
        ystyle: 1, $
        yrange: reverse(minmax(mlon_bins)), $
        ztitle: 'Upward current (A)', $
        zlog: 0 , $
        zrange: [0,.5e5], $
        yticklen: -0.02, $
        xticklen: -0.02 }


test = 1 
    ofn = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot','case2_ewogram_of_upward_current_global.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=6, ysize=4
    margins = [10,4,12,2]
    tpos = sgcalcpos(1, margin=margins, xchsz=xchsz, ychsz=ychsz)
    loadct2, 40
    device, decomposed=0
    tplot, ewo_var, position=tpos, trange=time_range
    device, decomposed=1

    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    msg = 'EWOgram of upward vertical current in '+strjoin(string(mlat_range,format='(I0)'),'-')+' deg MLat'
    xyouts, tx,ty,/normal, msg

    mlt_times = time_double(['2014-08-28/10:10:40','2014-08-28/10:47:25'])
    mlts = [0.4,6.5]
    ;mlts = [0.1,5.5]
    mlons = [0,0]
    foreach mlt_time, mlt_times, ii do begin
        et = convert_time(mlt_time, from='unix', to='epoch')
        mlt_bins = slon2lt(mlon_bins, et, /deg, /mag)
        mlons[ii] = interpol(mlon_bins, mlt_bins, mlts[ii])
    endforeach

    mlon_plot_range = reverse(minmax(mlon_bins))
    plot, time_range, mlon_plot_range, $
        xstyle=5, ystyle=5, yrange=mlon_plot_range, $
        /nodata, /noerase, position=tpos
    ;oplot, mlt_times, mlons, color=sgcolor('white'), linestyle=1
    ;oplot, mlt_times, mlons-0.2*15, color=sgcolor('white'), linestyle=1
    ;oplot, mlt_times-2.5*60, mlons-0.2*15, color=sgcolor('white'), linestyle=1
    ;oplot, mlt_times+2.5*60, mlons-0.2*15, color=sgcolor('white'), linestyle=1

    asi_times = time_double(['2014-08-28/10:12:10','2014-08-28/10:12:48'])
    mlons = [-78,-88]
    oplot, asi_times, mlons, color=sgcolor('white'), linestyle=0
    mlons = [-78,-65]
    oplot, asi_times, mlons, color=sgcolor('white'), linestyle=0
    if keyword_set(test) then stop
    sgclose

    stop


;---Make plots.
    ;test = 0
    ct = 70
    max_j = 2.5e5
    margins = [8,4,10,1]
    xticklen = -0.015
    colors = reverse(findgen(255))
    fig_ysize = 3
    time_range = minmax(times)


 test = 1
 test_time = time_double('2014-08-28/10:10')
    ;---Mlonlat.
    fig_xsize = 6
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    tpos = sgcalcpos(1, margins=margins, xchsz=xchsz, ychsz=ychsz)
    yticklen = xticklen/(tpos[2]-tpos[0])*(tpos[3]-tpos[1])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    fig_dir = join_path([homedir(),'vertical_current_mlonlat'])
    fig_files = strarr(ntime)
    jdata_m = j_new
    for ii=0,ntime-1 do begin
        if keyword_set(test) then if times[ii] lt test_time then continue
        ofn = join_path([fig_dir,'vertical_current_mlonlat_'+time_string(times[ii],tformat='YYYY_MMDD_hhmm_ss')+'.png'])
        if keyword_set(test) then ofn = 0
        sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
        erase, sgcolor('white')
        sgtv, bytscl(reform(-jdata_m[ii,*,*]),min=-max_j,max=max_j), ct=ct, position=tpos, /resize
        plot, mlon_bins, mlat_bins, /nodata, /noerase, $
        xstyle=1, xticklen=xticklen, xtitle='MLon (deg)', $
        ystyle=1, yticklen=yticklen, ytitle='MLat (deg)', $
        position=tpos
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty, /normal, time_string(times[ii])+' UT'
        sgcolorbar, colors, ct=ct, zrange=[-1,1]*max_j, ztitle='Vertical J (A), red (>0) for upward', position=cbpos
        if keyword_set(test) then stop
        sgclose
        fig_files[ii] = ofn
    endfor

    mov_file = join_path([homedir(),'vertical_current_mlonlat_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_')+'.mp4'])
    spic2movie, fig_dir, mov_file, 'png'

end
