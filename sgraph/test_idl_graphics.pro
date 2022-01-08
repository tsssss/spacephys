; compare to window render, raster outputs (png, jpg) are the same.
; pdf is similar, but file size is largest.
; ps and eps are different.
; conclusion: usable, not practical.

; set to use software rendering, it is more reliable than hardware rendering.
pref_set, strupcase('idl_gr_'+!d.name+'_renderer'), 1, /commit

; test data.
x = findgen(200)
y = sin(2*!dpi/25*x)*exp(-0.02*x)

w = window(dimensions = [400,300], font_name = 'Times', font_style = 'Normal', font_size = 12)
extra = {position:[0.2,0.15,0.95,0.9],xmajor:5,xminor:5,$
    xrange:[0,200], yrange:[-1,1], $
    xtitle:'Mass $\mu\Phi$', ytitle:'Time $\tau\overline{15}$'}
p = plot(x,y, /overplot, _extra = extra)
(p.axes)[0].ticklen = 1
(p.axes)[0].subticklen = 0.05
(p.axes)[0].gridstyle = 2
(p.axes)[0].tickfont_name = 'Times'
(p.axes)[1].ticklen = 1
(p.axes)[1].subticklen = 0.05
(p.axes)[1].gridstyle = 2
(p.axes)[1].tickfont_name = 'Times'

px2inch = 1d/101.6d
p.save, 'idl_graphics.pdf', height = 300*px2inch
p.save, 'idl_graphics.png', height = 300
p.save, 'idl_graphics.jpg', height = 300
p.save, 'idl_graphics.ps', height = 300
p.save, 'idl_graphics.eps', height = 300
end