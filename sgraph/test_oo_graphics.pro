; compare to window render, raster outputs (png, jpg) are the same.
; pdf and eps are almost same, font is slightly smaller.
; conclusion: very usable, but need to set everything yourself.

; set to use software rendering, it is more reliable than hardware rendering.
pref_set, strupcase('idl_gr_'+!d.name+'_renderer'), 1, /commit

; test data.
x = findgen(200)
y = sin(2*!dpi/25*x)*exp(-0.02*x)

; window settings.
xsz = 400 & ysz = 300
pos = [0.2,0.15,0.95,0.9]

w = obj_new('IDLgrWindow', retain = 2, dimensions = [xsz,ysz])
v = obj_new('IDLgrView')
m = obj_new('IDLgrModel')

; viewplane_rect for plot is in data unit, [x0,y0,dx,dy].
xr = [0d,200] & yr = [-1d,1] & zr = [0d,0]
rec = dblarr(4)
rec[2] = (xr[1]-xr[0])/(pos[2]-pos[0])
rec[0] = xr[0]-rec[2]*pos[0]
rec[3] = (yr[1]-yr[0])/(pos[3]-pos[1])
rec[1] = yr[0]-rec[3]*pos[1]
v.setproperty, viewplane_rect = rec

font = obj_new('IDLgrFont', 'Times', size = 10)
;font = obj_new('IDLgrFont', 'Symbol', size = 10)
xtitle = obj_new('IDLgrText', 'Mass '+tex2idl('$\mu\Phi$'), font = font, enable_formatting = 1)
ytitle = obj_new('IDLgrText', 'Time '+tex2idl('$\tau\overline{15}$'), font = font, enable_formatting = 1)

p = obj_new('IDLgrPlot', x,y, xrange = xr, yrange = yr)
xticklen = 0.03*(yr[1]-yr[0]) & yticklen = 0.03*(xr[1]-xr[0])
axb = obj_new('IDLgrAxis', 0, title = xtitle, $
    location = [xr[0],yr[0],zr[0]], gridstyle = 2, $
    range = [0,200], tickdir = 0, ticklen = yr[1]-yr[0], $
    subticklen = 0.6*0.03, /exact)
axt = obj_new('IDLgrAxis', 0, $
    location = [xr[0],yr[1],zr[0]], $
    range = [0,200], tickdir = 1, ticklen = xticklen, $
    subticklen = 0.6, tickformat='(A1)', /exact)
ayl = obj_new('IDLgrAxis', 1, title = ytitle, $
    location = [xr[0],yr[0],zr[0]], gridstyle = 2, $
    range = [-1,1], tickdir = 0, ticklen = xr[1]-xr[0], $
    subticklen = 0.6*0.03, /exact)
ayr = obj_new('IDLgrAxis', 1, $
    location = [xr[1],yr[0],zr[0]], $
    range = [-1,1], tickdir = 1, ticklen = yticklen, $
    tickformat = '(A1)', /exact)
v->add, m
m->add, p
m->add, axb & m->add, axt & m->add, ayl & m->add, ayr
axb.getproperty, ticktext = tmp
tmp.setproperty, font = font
axb.setproperty, ticktext = tmp
ayl.getproperty, ticktext = tmp
tmp.setproperty, font = font
ayl.setproperty, ticktext = tmp

; save to png.
w->draw, v
w->getproperty, image_data = im
write_png, 'oographics.png', im

; save to pdf.
px2inch = 1d/101.6d
pdf = obj_new('IDLgrPDF', dimensions = [xsz,ysz]*px2inch, units = 1)
pdf->addpage
pdf->draw, v
pdf->save, 'oographics.pdf'

; save to eps.
clip = obj_new('IDLgrClipboard', dimensions = [xsz,ysz], units = 0)
clip->draw, v, filename = 'oographics1.eps', /postscript
clip->draw, v, filename = 'oographics2.eps', /postscript, /vector
stop
obj_destroy, [v,w,pdf,clip]
end