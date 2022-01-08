;+
; Figure for plotting Dst in a wider range and zoom in to Dst and AE.
;-
;


    utr0 = time_double(['2014-08-25','2014-09-05'])
    utr1 = time_double(['2014-08-28','2014-08-28/15:00'])
    times = time_double('2014-08-28/'+['05:00','06:00','10:00','11:00'])


    vars = ['dst','ae']
    load = 0
    foreach tvar, vars do if tnames(tvar) eq '' then load = 1
    if load then omni_read_index, utr0, index=vars
    
    var = 'dst'
    options, var, 'constant', [0,-50]
    options, var, 'yrange', [-100,20]
    options, var, 'yticks', 2
    options, var, 'ytickv', [-2,-1,0]*50
    options, var, 'yminor', 5
    
    var = 'ae'
    options, var, 'constant', 1000
    options, var, 'yrange', [0,1500]
    options, var, 'yticks', 3
    options, var, 'ytickv', [3,2,1]*500
    options, var, 'yminor', 5
    

    ofn = sparentdir(srootdir())+'/plot/fig_omni_overview.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=5, ysize=4, /inch
    device, decomposed=0
    loadct2, 43
    
    red = 6
    thick = (size(ofn,/type) eq 7)? 12: 4
    
    poss = sgcalcpos(3,tmargin=2,bmargin=5,rmargin=5)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    pos1 = poss[*,0]
    pos2 = pos1 & pos2[1] = poss[1,2] & pos2[3] = pos1[1]-ychsz*4
    poss = sgcalcpos(2,position=pos2)
    tplot, 'dst', trange=utr0, /novtitle, position=pos1, /noerase
    xyouts, xchsz*2, pos1[3]-ychsz*0.8, /normal, 'a. Dst'
    timebar, utr1, color=red
    txs = (utr1-utr0[0])/(utr0[1]-utr0[0])*(pos1[2]-pos1[0])+pos1[0]
    plots, [txs[0],poss[0,0]], [pos1[1],poss[3,0]], color=red
    plots, [txs[1],poss[2,0]], [pos1[1],poss[3,0]], color=red
    
    tplot, vars, trange=utr1, position=poss, /novtitle, /noerase
    xyouts, xchsz*2, poss[3,0]-ychsz*0.8, /normal, 'b. Dst zoom-in'
    xyouts, xchsz*2, poss[3,1]-ychsz*0.8, /normal, 'c. AE'
    
    txs = (times-utr1[0])/(utr1[1]-utr1[0])*(poss[2,0]-poss[0,0])+poss[0,0]
    ty = poss[3,1]
    plots, txs[0:1], ty+[0,0], /normal, thick=thick, color=blue
    plots, txs[2:3], ty+[0,0], /normal, thick=thick, color=blue
    
    chsz = 0.8
    xyouts, mean(txs[0:1]), ty+ychsz*0.7, /normal, 'Case 1', alignment=0.5, charsize=chsz
    xyouts, mean(txs[2:3]), ty+ychsz*0.7, /normal, 'Case 2', alignment=0.5, charsize=chsz

    sgclose

end