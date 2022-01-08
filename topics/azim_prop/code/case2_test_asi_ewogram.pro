;+
; Make a EWOgram for asi.
;-

    time_range = time_double(['2014-08-28/10:00','2014-08-28/10:50'])
    plot_mlon_range = [0,-120]
    mlat_range = [60,70]
    time_range = time_double(['2014-08-28/10:05','2014-08-28/10:20'])
    plot_mlon_range = [-60,-100]


    asi_var = 'thg_mlonimg'
    get_data, asi_var, times, mlon_image, limits=lim
    mlat_bins = lim.mlat_bins
    mlon_bins = lim.mlon_bins

    mlat_index = lazy_where(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    ntime = n_elements(times)
    nmlon_bin = n_elements(mlon_bins)
    ewo = fltarr(ntime,nmlon_bin)
    foreach time, times, ii do begin
        foreach mlon, mlon_bins, jj do begin
            tmp = reform(mlon_image[ii,jj,mlat_index])
            ewo[ii,jj] = mean(tmp)
        endforeach
    endforeach
    ewo_var = 'thg_mlonimg_ewo'
    store_data, ewo_var, times, ewo, mlon_bins, limits={$
        spec: 1, $
        no_interp: 1, $
        ytitle: 'MLon (deg)', $
        ystyle: 1, $
        yrange: plot_mlon_range, $
        ztitle: 'Photon count (#)', $
        zlog: 0 , $
        yticklen: -0.02, $
        xticklen: -0.02 }


test = 0
    ofn = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot','case2_ewogram_of_asf.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=6, ysize=3



    margins = [10,4,12,2]
    tpos = sgcalcpos(1, margin=margins, xchsz=xchsz, ychsz=ychsz)
    device, decomposed=0
    loadct2, 40
    tplot, ewo_var, position=tpos, trange=time_range
    device, decomposed=1
    
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    msg = 'EWOgram of auroral brightness in '+strjoin(string(mlat_range,format='(I0)'),'-')+' deg MLat'
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
    plot, time_range, plot_mlon_range, $
        xstyle=5, ystyle=5, yrange=plot_mlon_range, $
        /nodata, /noerase, position=tpos
    ;oplot, mlt_times, mlons, color=sgcolor('white'), linestyle=1
;    oplot, mlt_times, mlons-0.2*15, color=sgcolor('white'), linestyle=1
;    oplot, mlt_times-2.5*60, mlons-0.2*15, color=sgcolor('white'), linestyle=1
;    oplot, mlt_times+2.5*60, mlons-0.2*15, color=sgcolor('white'), linestyle=1

    asi_times = time_double(['2014-08-28/10:12:10','2014-08-28/10:12:48'])
    mlons = [-78,-88]
    oplot, asi_times, mlons, color=sgcolor('white'), linestyle=0
    mlons = [-78,-65]
    oplot, asi_times, mlons, color=sgcolor('white'), linestyle=0
    
    
    if keyword_set(test) then stop
    sgclose
    

end
