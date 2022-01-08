;+
; Check binned B data.
;-


test = 0
    fn = '/Volumes/GoogleDrive/My Drive/works/works/pflux_survey/data/binned_data/pflux_survey_spherical_spatial.cdf'
    probes = ['a','b']
    rgb = constant('rgb')
    
    plot_file = join_path([homedir(),'fig_check_binned_b0_sm.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, xsize=8, ysize=4, xchsz=xchsz, ychsz=ychsz
    
    xtitle = 'Mean Value per bin (nT)'
    xrange = [-1,1]*400
    yrange = [1,100]
    xticklen_chsz = -0.20
    yticklen_chsz = -0.40

    margins = [10,4,3,2]
    poss = sgcalcpos(1,2, margins=margins, xpad=3)
    
    tpos = poss[*,0]
    tx = tpos[0]+xchsz*1
    ty = tpos[1]+ychsz*3.2
    xyz = constant('xyz')
    for ii=0,2 do begin
        msg = 'SM B'+xyz[ii]
        xyouts, tx, ty-ychsz*(ii+0.2), /normal, msg, color=rgb[ii]
    endfor

    foreach probe, probes, probe_id do begin
        prefix = 'rbsp'+probe+'_'
        mean = cdf_read_var(prefix+'b0_sm_mean', filename=fn)
        stddev = cdf_read_var(prefix+'b0_sm_stddev', filename=fn)
        
        tpos = poss[*,probe_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytitle = (probe_id eq 0)? 'Stddev Value per bin (nT)': ''
        ytickformat = (probe_id eq 0)? '': '(A1)'
        
        ; Set up coord and box.
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=1, xlog=0, xtitle=xtitle, xrange=xrange, xticklen=xticklen, $
            ystyle=1, ylog=1, ytitle=ytitle, yrange=yrange, yticklen=yticklen, ytickformat=ytickformat, $
            position=tpos
        
        ; Plot per component.
        for ii=0,2 do begin
            color = rgb[ii]
            plots, mean[*,*,*,ii], stddev[*,*,*,ii], $
                psym=1, symsize=0.25, color=color
        endfor
        
        ; Figure label.
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*1.2
        xyouts, tx,ty,/normal, probe+'. RBSP-'+strupcase(probe)
        
        
        ; Fit a trend.
        txs = abs(mean)
        tys = stddev
        index = where(finite(txs))
    endforeach
    
    
    tpos = poss[*,0]
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    msg = 'Bin "model" B field (<0.4 mHz) using 3-year data (Oct 2012-2015) >4 Re'
    xyouts, tx,ty,msg, /normal
    
    if keyword_set(test) then stop
    sgclose

end