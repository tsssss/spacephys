function xml_find_tag, tags, atts
    ntag = n_elements(tags)
    natt = n_elements(atts)
    idx = -1
    for i = 0, ntag-1 do begin
        flag = 0
        for j = 0, natt-1 do flag += stregex(tags[i],atts[j]) ne -1
        if flag eq natt then idx = [idx,i]  ; hit all attribute.
    endfor
    if n_elements(idx) eq 1 then return, -1 else return, idx[1:*]
end

; Use ps2pdf in GhostScript by default.
; If use epstopdf in Latex, then set epstopdf.
pro polar_read_ke_flux, fn, tr, eleyr0, ionyr0, $
    epstopdf = epstopdf, pstopdf = pstopdf, ps2pdf = ps2pdf
    
    ; default position, ion below ele, but device coord is 0 at top.
    if keyword_set(epstopdf) then begin ; for epstopdf.
        ionpos = [84.246d,358.724,449.292,680.400]
        elepos = [84.246d, 27.977,449.292,349.653]
    endif else begin                    ; for ps2pdf & pstopdf.
        ionpos = [98.246d,380.724,463.292,702.400]
        elepos = [98.246d, 49.977,463.292,371.653]
    endelse

    ; device coord is upside down, i.e., 0 is at top-left corner.
    eleyr = -eleyr0 & ionyr = -ionyr0
    
    ; read svg file.
    tags = sread_xml_parse_xml(fn)
    ntag = n_elements(tags)
    
    ; determine mode.
    is_ai = stregex(tags[1],'Illustrator') ge 0
    if is_ai then begin     ; svg by illustrator.
        ; find dash line.
        atts = ['stroke-dasharray="1.1339,2.8347"','stroke="#000000"']
        idx = xml_find_tag(tags,atts)
        ; find ion line plot, should find 2 tags that match.
        atts = ['stroke="#00FF00"','points=']   ; ion: green.
        idx = xml_find_tag(tags,atts)
        ; find the longest data series, because I manually combine separate
        ; segment into one whole series, or if originally there is only one
        ; series, there are two same lines.
        tmp = strlen(tags[idx]) & tmp = (where(tmp eq max(tmp)))
        idx = idx[tmp[0]]
        tmp = sread_xml_parse_tag(tags[idx], 'points', attribute = att)
        if n_elements(att) eq 0 then $
            tmp = sread_xml_parse_tag(tags[idx], 'd', attribute = att)
        tmp = double(strsplit(att.value,', "',/extract))
        x0 = tmp[0:*:2] & y0 = tmp[1:*:2]
        pos = ionpos & yr = ionyr & xr = time_double(tr)
        ionke = (y0-pos[1])/(pos[3]-pos[1])*(yr[1]-yr[0])+yr[0]
        t0 = (x0-pos[0])/(pos[2]-pos[0])*(xr[1]-xr[0])+xr[0]
        store_data, 'po_ion_keflux', t0, ionke, $
            limit = {yrange:ionyr0,ytitle:'(mW/m!U2!N)',ylabel:'KEi'}

        ; find electron line plot, should find 2 tags that match.
        atts = ['stroke="#FFFF00"','points=']   ; electron: yellow.
        idx = xml_find_tag(tags,atts)
        tmp = strlen(tags[idx]) & tmp = (where(tmp eq max(tmp)))
        idx = idx[tmp[0]]
        tmp = sread_xml_parse_tag(tags[idx], 'points', attribute = att)
        if n_elements(att) eq 0 then $
            tmp = sread_xml_parse_tag(tags[idx], 'd', attribute = att)
        tmp = double(strsplit(att.value,', "',/extract))
        x0 = tmp[0:*:2] & y0 = tmp[1:*:2]
        pos = elepos & yr = eleyr & xr = time_double(tr)
        eleke = (y0-pos[1])/(pos[3]-pos[1])*(yr[1]-yr[0])+yr[0]
        t0 = (x0-pos[0])/(pos[2]-pos[0])*(xr[1]-xr[0])+xr[0]
        store_data, 'po_ele_keflux', t0, eleke, $
            limit = {yrange:eleyr0,ytitle:'(mW/m!U2!N)',ylabel:'KEe'}
    endif else begin        ; svg by inkscape.
        ; to do.
    endelse
;    vars = ['po_ion_ke','po_ele_ke']
;    tplot, vars, trange = tr
end

tr = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
eleyr = 0.13*[-1,1] & ionyr = 0.25*[-1,1]
fn = '/Volumes/Works/works/ps_svg/1998_0925_ke_ai_from_pdf.svg'
fn = '/Volumes/Works/works/ps_svg/final/1998_0925_hydra_ke_special.svg'

tr = time_double(['1998-10-01/02:00','1998-10-01/04:00'])
eleyr = 0.12*[-1,1] & ionyr = 0.08*[-1,1]
fn = '/Volumes/Works/works/ps_svg/final/1998_1001_hydra_ke_special.svg'

tr = time_double(['1998-10-03/02:00','1998-10-03/04:30'])
eleyr = 0.12*[-1,1] & ionyr = 0.08*[-1,1]
fn = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun/1998_1001_03/1998_1001_03_ke_special.svg'

tr = time_double(['2000-09-02/07:00','2000-09-02/07:30'])
eleyr = 0.1*[-1,1] & ionyr = 0.03*[-1,1]
fn = googledir()+'/works/works/cusp/cusp list polar 2-4Re/2000_0902_07/2000_0902_07_ke_special.svg'

polar_read_ke_flux, fn, tr, eleyr, ionyr
end