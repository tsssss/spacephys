;+
; Event distribution.
;-

;---Constants.
    secofday = 86400d   ; sec.
    cns = -1

    rad = !dpi/180d
    deg = 180d/!dpi

;---Settings.
    xrng = [0,-20]      ; in Re.
    yrng = [-1,1]*20    ; in Re.
    zrng = [4,6]
    xtitl = 'X (Re)'
    ytitl = 'Y (Re)'
    top = 254
    ct0 = 40
    
    ztitl = 'Count per bin'
    ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')
    
    bintype = 'dis_mlt_dns'
    
    ; file names.    
    rootdir = shomedir()+'/global_e_distr'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
    stryrs = string(smkarthm(2013,2015,1,'dx'),format='(I4)')
    probes = ['a','d','e']
    foreach tprobe, probes do foreach stryr, stryrs do begin
        
        datadir = rootdir+'/data/'+stryr
        if file_test(datadir,/directory) eq 0 then file_mkdir, datadir

        logfn = rootdir+'/global_e_distr_bin_th'+tprobe+'_'+bintype+'_'+stryr+'.log'
        if file_test(logfn) then file_delete, logfn & stouch, logfn

        binfn = datadir+'/global_e_distr_bin_th'+tprobe+'_'+bintype+'_'+stryr+'.sav'
        utfn = datadir+ '/global_e_distr_bin_th'+tprobe+'_'+stryr+'.sav'
        posfn = datadir+'/global_e_distr_bin_th'+tprobe+'_'+stryr+'_'+bintype+'.sav'

        restore, binfn
        
        
        tmp = size(bininfos,/dimensions)
        nbin1 = tmp[0]
        nbin2 = tmp[1]
        nbin3 = tmp[2]
        bincnts = dblarr(nbin1,nbin2,nbin3)
        for i=0,nbin1-1 do for j=0,nbin2-1 do for k=0,nbin3-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            bincnts[i,j,k] = n_elements(*bininfos[i,j,k].uts)
        endfor
        
        tcnt = total(bincnts,3)
        tdat = alog10(tcnt)
        tdat = bytscl(tdat, min=zrng[0],max=zrng[1], top=top)
        
        ofn = 0
        ofn = rootdir+'/fig/th'+tprobe+'_'+stryr+'_data_coverage.pdf'
        sgopen, ofn, xsize=4, ysize=8, /inch
        
        pos0 = [0.2,0.1,0.9,0.8]
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        device, decomposed=0
        loadct, ct0
        plot, xrng, yrng, /nodata, position=pos0, $
            xstyle=5, xrange=xrng, xtitle=xtitl, $
            ystyle=5, yrange=yrng, ytitle=ytitl, $
            /noerase, /isotropic
        
        bin2s = bin2s*15*rad
        for i=0, nbin1-1 do for j=0, nbin2-1 do begin
            if tcnt[i,j] eq 0 then continue
            txs =-[bin1s[i]*cos(bin2s[j]),bin1s[i]*cos(bin2s[j+1]),bin1s[i+1]*cos(bin2s[j+1]),bin1s[i+1]*cos(bin2s[j]),bin1s[i]*cos(bin2s[j])]
            tys = [bin1s[i]*sin(bin2s[j]),bin1s[i]*sin(bin2s[j+1]),bin1s[i+1]*sin(bin2s[j+1]),bin1s[i+1]*sin(bin2s[j]),bin1s[i]*sin(bin2s[j])]
            polyfill, txs, tys, /data, color=tdat[i,j]
        endfor
        
;        ; draw earth.
;        tmp = findgen(51)/50*2*!dpi
;        txs = cos(tmp)
;        tys = sin(tmp)
;        idx = where(txs le 0)
;        polyfill, txs[idx], tys[idx], /line_fill, orientation=45
;        polyfill, txs[idx], tys[idx], /line_fill, orientation=-45
;        plots, txs, tys, /data, color=tcolor, thick=thick
        
        loadct, 0
        titl = 'TH-'+strupcase(tprobe)+', Year='+stryr
        plot, xrng, yrng, /nodata, position=pos0, $
            xstyle=1, xrange=xrng, xtitle=xtitl, $
            ystyle=1, yrange=yrng, ytitle=ytitl, $
            /noerase, /isotropic

        pos1 = pos0 & pos1[1]=pos0[3]+ychsz*0.5 & pos1[3]=pos1[1]+ychsz*0.5
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle=ztitl, $
            ztickv=ztickv, ztickname=ztickn
        
        xyouts, (pos0[0]+pos0[2])*0.5, pos1[3]+ychsz*3, /normal, alignment=0.5, titl
        
        sgclose
    endforeach

end
