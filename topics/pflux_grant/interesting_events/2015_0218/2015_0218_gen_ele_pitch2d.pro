

log = 1
unit = 'energy'
types = ['electron','proton','oxygen']
zrngs = [[4,8],[3.5,6],[3.5,6]]
types = ['electron']
zrngs = [[4,10]]


; 3 periods for RBSP-A.
probe = ['a']
time_range = time_double(['2015-02-18/01:30','2015-02-18/04:00'])
time_range = time_double(['2015-02-18/02:00','2015-02-18/02:15'])


hopel3 = sread_rbsp_hope_l3(time_range, probes=probe)
times = make_bins(time_range, 12, /inner)

foreach type, types, type_id do begin
    foreach time, times do begin
        plot_hope_l3_pitch2d, time, type, unit=unit, $
                log=log, hopel3=hopel3, probe=probe, zrange=zrngs[*,type_id]
    endforeach
endforeach

end
