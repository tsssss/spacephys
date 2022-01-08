
pro test_sgraph_axis

    nrec = 101
    x = findgen(nrec)*2*!const.pi/(nrec-1)
    y1 = 3*sin(x)
    y2 = sin(x)+1
    
    !x.style = 1
    !y.style = 1
    plot, x, y1, position = [0.1,0.1,0.9,0.9], ystyle = 9
    plots, replicate(x[nrec-1],2), [0,-3], linestyle = 5
    print, !y.range
    print, !y.crange
    stop
    
    plot, x, y2, position = [0.1,0.5,0.9,0.9], /noerase, $
        xstyle = 5, ystyle = 5
    axis, x[nrec-1], yaxis = 1, yrange = [0,2], ystyle = 1, /save
    oplot, x, y2

end