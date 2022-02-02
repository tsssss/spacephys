;+
; Plot spatial distribution of cusp search results.
;-

test = 0
data_path = join_path([googledir(),'data','cusp_polar'])
data_file = join_path([data_path,'cusp_polar_cusp_sections.tplot'])
tplot_restore, filename=data_file

cusp_infos = get_var_data('section_with_cusp')
ncusp = n_elements(cusp_infos)

mlt_range = [-1d,1]*12+12
mlt_bin_size = 1d
ilat_range = [60d,90]
ilat_bin_size = 1.5d
dis_range = [1d,10]
dis_bin_size = 0.5d

mlt_bins = smkarthm(0.5,23.5, mlt_bin_size, 'dx')
nmlt_bin = n_elements(mlt_bins)-1
mlt_vals = mlt_bins[0:nmlt_bin-1]+mlt_bin_size*0.5

ilat_bins = make_bins(ilat_range, ilat_bin_size)
nilat_bin = n_elements(ilat_bins)-1
ilat_vals = ilat_bins[0:nilat_bin-1]+ilat_bin_size*0.5

dis_bins = make_bins(dis_range, dis_bin_size)
ndis_bin = n_elements(dis_bins)-1
dis_vals = dis_bins[0:ndis_bin-1]+dis_bin_size*0.5

cusp_counts = fltarr(nmlt_bin,nilat_bin,ndis_bin)
foreach cusp_info, cusp_infos do begin
    if product(cusp_info.mlt-mlt_range) gt 0 then continue
    if product(cusp_info.ilat-ilat_range) gt 0 then continue
    if product(cusp_info.dis-dis_range) gt 0 then continue

    xids = floor(minmax((cusp_info.mlt-min(mlt_range))/mlt_bin_size))
    yids = floor(minmax((cusp_info.ilat-min(ilat_range))/ilat_bin_size))
    zids = floor(minmax((cusp_info.dis-min(dis_range))/dis_bin_size))

    cusp_counts[xids[0]:xids[1],yids[0]:yids[1],zids[0]:zids[1]] += 1
endforeach


plot_file = join_path([srootdir(),'fig_cusp_search_plot_spatial_distr.pdf'])
if keyword_set(test) then plot_file = 0
xsize = 6
ysize = 2
sgopen, plot_file, xsize=xsize, ysize=ysize, xchsz=xchsz, ychsz=ychsz

xticklen_chsz = -0.30
yticklen_chsz = -0.40

;;---Plot a slice at noon.
;    tpos = sgcalcpos(1)
;    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
;    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
;
;    xid = nmlt_bin/2-1
;    the_slice = reform(cusp_counts[xid,*,*])
;    xrange = [5,0d]
;    xtitle = 'X (Re)'
;    xminor = 2
;    xtickv = make_bins(xrange, xminor, inner=1)
;    xticks = n_elements(xtickv)-1
;
;    yrange = [0d,10]
;    ytitle = 'Z (Re)'
;    yminor = 2
;    ytickv = make_bins(yrange, yminor, inner=1)
;    yticks = n_elements(ytickv)-1
;
;    zrange = [1,1e3]
;    ztitle = 'Count (#)'
;    zrange_log = alog10(zrange)
;    color_range = [10,250]
;    color_table = 50
;
;    plot, xrange, yrange, $
;        xstyle=5, xrange=xrange, xtitle=xtitle, $
;        ystyle=5, yrange=yrange, ytitle=ytitle, $
;        position=tpos, nodata=1, noerase=1, iso=1
;
;
;    for xid=0,nilat_bin-1 do begin
;        the_ilat = ilat_bins[xid:xid+1]*constant('rad')
;        for yid=0,ndis_bin-1 do begin
;            the_dis = dis_bins[yid:yid+1]
;            p1 = the_dis[0]*[cos(the_ilat[0]),sin(the_ilat[0])]
;            p2 = the_dis[0]*[cos(the_ilat[1]),sin(the_ilat[1])]
;            p3 = the_dis[1]*[cos(the_ilat[1]),sin(the_ilat[1])]
;            p4 = the_dis[1]*[cos(the_ilat[0]),sin(the_ilat[0])]
;            color = interpol(color_range, zrange_log, alog10(the_slice[xid,yid])>0)
;            color = sgcolor(color, ct=color_table)
;            if the_slice[xid,yid] eq 0 then continue
;            polyfill, data=1, color=color, $
;                [p1[0],p2[0],p3[0],p4[0]], $
;                [p1[1],p2[1],p3[1],p4[1]]
;        endfor
;    endfor
;
;    plot, xrange, yrange, $
;        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
;        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
;        position=tpos, nodata=1, noerase=1, iso=1
;
;;---Draw earth.
;    tmp = smkarthm(0,2*!dpi,40,'n')
;    txs = cos(tmp)
;    tys = sin(tmp)
;    oplot, txs, tys


    count_yrange = [0,1.3e4]
    margins = [8,4,2,1]
    npanel = 3
    poss = sgcalcpos(1,npanel, xpad=1, margins=margins)
    fig_labels = letters(npanel)+'.'
    for ii=0,npanel-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty,normal=1, fig_labels[ii]
    endfor

;---Histogram of ILat.
    plot_mlt_range = [11,13]
    plot_mlt_range = mlt_range
    xids = lazy_where(mlt_vals, '[]', plot_mlt_range)
    hist = total(total(cusp_counts[xids,*,*],1),2)
    bins = ilat_vals
    nbin = n_elements(bins)
    binsize = ilat_bin_size

    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xrange = ilat_range
    xtitle = 'ILat (deg)'
    xminor = 10
    xtickv = make_bins(xrange, xminor, inner=1)
    xticks = n_elements(xtickv)-1

    yrange = [0,max(hist)]
    yrange = count_yrange
    ytitle = 'Count (#)'
    ystep = 5e3
    yminor = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, nodata=1, noerase=1
    for ii=0,nbin-1 do begin
        vals = bins[ii]+[-1,1]*0.5*binsize
        oplot, vals, hist[ii]+[0,0]
        foreach val, vals do begin
            oplot, val+[0,0], [yrange[0],hist[ii]]
        endforeach
    endfor

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        position=tpos, nodata=1, noerase=1

;---Histogram of MLT.
    plot_dis_range = [1,3]
    plot_dis_range = dis_range
    zids = lazy_where(dis_vals, '[]', plot_dis_range)
    hist = total(total(cusp_counts[*,*,zids],3),2)
    bins = mlt_vals
    nbin = n_elements(bins)
    binsize = mlt_bin_size

    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xrange = [0.5,23.5]
    xtitle = 'MLT (h)'
    xstep = 6
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep

    yrange = [0,max(hist)]
    ytitle = 'Count (#)'
    yrange = count_yrange
    ytitle = ' '
    ytickformat = '(A1)'
    ystep = 5e3
    yminor = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1


    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, nodata=1, noerase=1
    for ii=0,nbin-1 do begin
        vals = bins[ii]+[-1,1]*0.5*binsize
        oplot, vals, hist[ii]+[0,0]
        foreach val, vals do begin
            oplot, val+[0,0], [yrange[0],hist[ii]]
        endforeach
    endfor
    oplot, 12+[0,0], yrange, linestyle=1

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickformat=ytickformat, $
        position=tpos, nodata=1, noerase=1

;---Histogram of dis.
    hist = total(total(cusp_counts,1),1)
    bins = dis_vals
    nbin = n_elements(bins)
    binsize = dis_bin_size

    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xrange = dis_range
    xtitle = '|R| (Re)'
    xstep = 3
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep

    yrange = [0,max(hist)]
    ytitle = 'Count (#)'
    yrange = count_yrange
    ytitle = ' '
    ytickformat = '(A1)'
    ystep = 5e3
    yminor = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1


    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, nodata=1, noerase=1
    for ii=0,nbin-1 do begin
        vals = bins[ii]+[-1,1]*0.5*binsize
        oplot, vals, hist[ii]+[0,0]
        foreach val, vals do begin
            oplot, val+[0,0], [yrange[0],hist[ii]]
        endforeach
    endfor

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickformat=ytickformat, $
        position=tpos, nodata=1, noerase=1

if keyword_set(test) then stop
sgclose

end
