;+
; Load project info, if no file exists, then initialize it.
; Keep updating as the project goes on.
;
; project. The basic project dictionary returned from init_project.
;-

pro field_aligned_conjunction_init_project, project

;---Information on searching events.
    ; mlon and mlat in deg.
    ; radius in km.
    radius = 1e3    ; km, roughly 10 deg.
    site = 'pokr'
    site_info = themis_read_mag_metadata(sites=site)
    search = dictionary($
        'site', site, $
        'mlon', site_info.mlon, $
        'mlat', site_info.mlat, $
        'radius', radius, $
        'probes', ['rbsp'+letters('b'),'arase','swarm'+letters('c'),'dmsp'+string([16,17,18],format='(I0)')]))
    project.search = search
    update_project, project

end

project = load_project('field_aligned_conjunction')
end