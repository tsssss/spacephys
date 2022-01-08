;+
; Treat only latin letters, symbol will be dealt with latex.
; Cannot change font system because that should be determined outside here.
;-
function sgfont, str, fname, help = help

    if n_elements(fname) eq 0 then fname = 'helvetica'
    
    tmp = ['',' italic',' bold',' bold italic']
    ttnames = ['helvetica'+tmp, 'times'+tmp, 'courier'+tmp]
    ttcodes = [3,5,4,6, 7,8,15,16, 11,12,13,14]
    ttxpads = [0.3,0.3,0.3,0.3, 0.3,0.3,0.3,0.3, 0.3,0.3,0.3,0.3]
    ttypads = [0.4,0.4,0.4,0.4, 0.4,0.4,0.4,0.4, 0.4,0.4,0.4,0.4]
    
    psnames = ['helvetica'+['',' bold',' narrow',' narrow bold oblique'], $
        'times'+[' roman',' bold italic'], 'courier'+['',' oblique'], $
        'palatino'+['',' italic',' bold',' bold italic'], 'zapf dingbats', $
        'new centry schoolbook'+['',' bold'], 'avant garde book']
    pscodes = [3,4,5,6, 7,8, 11,12, 13,14,15,16, 10, 18,19, 17]
    
    henames = ['simplex roman','simplex greek','duplex roman','complex roman', $
        'complex greek','complex italic','gothic english','simplex script', $
        'complex script','gothic italian','gothic german','cyrillic', $
        'triplex roman','triplex italic']
    hecodes = [3,4,5,6,7,8,11,12,13,14,15,16,17,18]
    
    case !p.font of
        1: begin fnames = ttnames & fcodes = ttcodes & end
        0: begin fnames = psnames & fcodes = pscodes & end
        -1:begin fnames = henames & fcodes = hecodes & end
    endcase
    
    if keyword_set(help) then $
        return, transpose(string(fcodes,format='(I2)'))+'    '+transpose(fnames)
    
    idx = where(fname eq fnames, cnt)
    if cnt eq 0 then fcode = 3 else fcode = fcodes[idx] ; fall to default font.
    return, '!'+string(fcode,format='(I0)')+str+'!X'

end