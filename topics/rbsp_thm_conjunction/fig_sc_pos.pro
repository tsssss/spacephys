
    ;_2014_1222_load_data
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr2 = time_double(['2014-12-22/05:25'])
    pres = ['rba','rbb','tha','the','thd']+'_'
    
    ofn = shomedir()+'/fig_sc_pos.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=5, ysize=5, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    x0rng = [2,-8]
    y0rng = [-5,15]
    z0rng = [-4,2]
    
    x0title = 'X GSM (Re)'
    y0title = 'Y GSM (Re)'
    z0title = 'Z GSM (Re)'
    
    pos0 = [.2,.15,.6,.95]
    pos1 = [.6,.15,.9,.95]+[1,0,1,0]*xchsz*2
    
    colors = [1,2,4,0,6]
    psym = 6
    hsize = (size(ofn,/type) eq 7)? 100: 5
    
    
;---X-Y plane.
    xrng = x0rng
    yrng = y0rng
    xtitle = x0title
    ytitle = y0title
    idx0 = 0
    idx1 = 1
    tpos = pos0
    plot, xrng, yrng, position=tpos, /nodata, $
        xstyle=1, xtitle=xtitle, xticklen=-0.01, $
        xticks=5, xminor=2, $
        ystyle=1, ytitle=ytitle, yticklen=-0.02, $
        /noerase, /isotropic
    tmp = findgen(51)/50*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    oplot, txs, tys
    polyfill, txs<0, tys, linestyle=0, orientation=45
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1.2, /normal, 'a. GSM X-Y plane'

    for i=0,1 do oplot, txs*(i+1)*5, tys*(i+1)*5, linestyle=1

        
    foreach pre0, pres, i do begin
        get_data, pre0+'r_gsm', uts, rgsm
        tmp = min(uts-utr2[0], /absolute, idx)
        plots, rgsm[idx,idx0], rgsm[idx,idx1], color=colors[i], psym=psym
        tmp = convert_coord(rgsm[idx,idx0], rgsm[idx,idx1],/data,/to_normal)
        xyouts, tmp[0], tmp[1]-ychsz*1.2, /normal, alignment=0.5, strupcase(strmid(pre0,0,3)), color=colors[i]
        if pre0 eq 'the_' then begin
            get_data, pre0+'vbulk', uts, dat
            tutr = time_double(['2014-12-22/04:40','2014-12-22/05:40'])
            tidx = where(uts ge tutr[0] and uts le tutr[1])
            velx = mean(dat[tidx,0])
            vely = mean(dat[tidx,1])
            velz = mean(dat[tidx,2])
            print, velx,vely,velz
            vmag = sqrt(velx^2+vely^2)
            tr = 0.2
            tx = rgsm[idx,idx0]
            ty = rgsm[idx,idx1]
            arrow, tx,ty, tx+velx*tr, ty+vely*tr, /data, /solid, hsize=hsize
            tmp = convert_coord(tx+velx*tr, ty+vely*tr,/data,/to_normal)
            xyouts, tmp[0], tmp[1]+ychsz*0.2, /normal, alignment=0.5, $
                '<V>='+sgnum2str(vmag,nsgn=2)+' km/s'
        endif
    endforeach
    
    
;---Y-Z plane.
    xrng = z0rng
    yrng = y0rng
    xtitle = z0title
    ytitle = ''
    idx0 = 2
    idx1 = 1
    tpos = pos1
    plot, xrng, yrng, position=tpos, /nodata, $
        xstyle=1, xtitle=xtitle, xticklen=-0.01, $
        xticks=3, xminor=2, $
        ystyle=1, ytitle=ytitle, yticklen=-0.02, $
        ytickformat='(A1)', $
        /noerase, /isotropic
    tmp = findgen(51)/50*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    oplot, txs, tys
    polyfill, txs, tys, linestyle=0, orientation=45
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1.2, /normal, 'b. GSM Y-Z plane'
    
    foreach pre0, pres, i do begin
        get_data, pre0+'r_gsm', uts, rgsm
        tmp = min(uts-utr2[0], /absolute, idx)
        plots, rgsm[idx,idx0], rgsm[idx,idx1], color=colors[i], psym=psym
        tmp = convert_coord(rgsm[idx,idx0], rgsm[idx,idx1],/data,/to_normal)
        ;xyouts, tmp[0], tmp[1]+ychsz*0.8, /normal, alignment=0.5, strupcase(strmid(pre0,0,3)), color=colors[i]
    endforeach
    
        
    sgclose

end
