
pro sgsetcoord, pos, type = type, _extra = extra

    if n_elements(pos) ne 4 then pos = [0d,0,1,1]
    if n_elements(type) eq 0 then type = 'linear'
    
    case strlowcase(type) of
        'linear':
        'log':
        'polar':
        'footpoint': begin  ; (lat,time).
            if n_elements(minlat) eq 0 then minlat = 50
            lim = [abs(minlat)-1d-14,0,90,360]
            map_set, name = 'AzimuthalEquidistant', 90, 0, 0, /isotropic, $
                /noerase, /noborder, position = pos, limit = lim, $
                latdel = 10, londel = 45, glinestyle = 1, $
                color = sgcolor('black'), $
                label = 1, latlab = 45
                
            end
        'aurora':
    endcase
end