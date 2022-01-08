;+
; show fonts are the same size in all modes.
; show nominal character size is real character size plus kerns.
; show font size is 9pt in 8in wide paper.
;-
pro test_sgraph_font

    ofn = shomedir()+'/sg_font.'
    str = sgfont('2014WorldCup','times')
    str2 = sgfont('2014 World Cup is coming, I plan to watch the 1st game, '+$
        'which is Brazil vs Croatia. Here is a piece of news from the '+$
        'Internet.!C!CThe hosts are almost unbeatable at home, won the '+$
        'Confederations Cup last summer, and are loaded with talent, '+$
        'as usual.!CNothing less than a World Cup title a month from now '+$
        'will satisfy the fans, at least those who are not actively '+$
        'protesting!Cthe tournament.','times')
    black = sgcolor('black')
    xsize = 8 & ysize = 2    ; inch.

    sgwindow, 0, xsize = xsize, ysize = ysize, /inch
    
    x0 = !d.x_ch_size*0.7
    y0 = !d.y_ch_size*0.6
    xyouts, 0,0,/device, str, color = black
    xyouts, x0,y0,/device, str, color = black

    sgpsopen, ofn+'eps'
    sgtruecolor
    x0 = !d.x_ch_size*0.7
    y0 = !d.y_ch_size*0.6
    xyouts, 0,0,/device, str, color = black
    xyouts, x0,y0,/device, str, color = black
    xyouts, 0.5/xsize,1-0.5/ysize,/normal, str2, color = black
    sgpsclose
    
    sgzopen, ofn+'png'
    sgtruecolor
    x0 = !d.x_ch_size*0.7
    y0 = !d.y_ch_size*0.6
    xyouts, 0,0,/device, str, color = black
    xyouts, x0,y0,/device, str, color = black
    sgzclose

end