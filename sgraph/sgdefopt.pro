
function sgdefopt, func

    if n_elements(func) eq 0 then func = 'plot'
    
    black = sgcolor('black') & white = sgcolor('white')
    
    if strlowcase(func) eq 'plot' then begin
        opt = {isotropic: 0, noerase:1, $
            xlog: 0, ylog: 0, ynozero: 1, $
            bakcground: white, color: black, $
            xstyle:1, ystyle:1, zstyle:1}
            
;            ; plot special.
;            max_value, minvalue, nsum, polar
;            ; char.
;            font, charthick, [xyz]charsize, charsize
;            ; tick.
;            [xyz]minor, [xyz]ticks, [xyz]tickv, [xyz]tickname, [xyz]tickformat, $
;            [xyz]tickinterval, [xyz]ticklayout, [xyz]ticklen, ticklen, [xyz]tickunits
;            ; axis.
;            [xyz]range, [xyz]thick, 
;            ; position.
;            position, data, device, normal, [xyz]margin, clip, noclip, $
;            ; annotation.
;            psym, symsize, linestyle, gridstyle, thick, $
;            [xyz]title, title, subtitle, $
;            ; other.
;            t3d, zvalue, nodata, background.
    endif
    
    return, opt
end


x = findgen(101)/100*2*!const.pi
y = 1.5*sin(x)
extra = {xstyle:1, ystyle:1, xminor:1}
erase, sgcolor('white')
plot, x, y, _extra = sgdefopt()
;   {background:sgcolor('white'), color:sgcolor('black')}
;    background = sgcolor('white'), color = sgcolor('black')
end
