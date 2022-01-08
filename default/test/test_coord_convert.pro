
pro test_coord_convert

    !x.s = [0,1] & !y.s = [0,1]
    
    ; ![xy].window = [0,1].
    plot, [0,1], /nodata, position = pos, xstyle = 5, ystyle = 5

    !x.s = [0.5,0.5] & !y.s = [0.5,0.5]
    !x.type = 3 & !y.type = 3
    print, convert_coord(0,50,/data,/to_norm)
end