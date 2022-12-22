;+
; Check Chaston+2015 event on 2013-06-01 RBSP-A.
;-

time_range = time_double(['2013-06-01/03:00','2013-06-01/08:00'])
probe = 'a'

log = 1
unit = 'energy'
species = ['electron','proton','oxygen']
zrngs = [[4,8],[3.5,6],[3.5,6]]

hopel3 = sread_rbsp_hope_l3(time_range, probes=probe)
times = make_bins(time_range, 12, inner=1)

foreach species_name, species, type_id do begin
    foreach time, times do begin
        plot_hope_l3_pitch2d, time, species_name, unit=unit, $
            log=log, hopel3=hopel3, probe=probe, zrange=zrngs[*,type_id]
    endforeach
endforeach