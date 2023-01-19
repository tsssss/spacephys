utr = time_double(['2013-05-01/07:25','2013-05-01/07:50'])
uts = smkarthm(utr[0],utr[1],12,'dx')
;uts = time_double(['2013-05-01/07:36:43'])
log = 1
unit = 'energy'
types = ['electron','proton','oxygen','helium']
probes = ['b']

for k = 0, n_elements(probes)-1 do $
    for j = 0, n_elements(types)-1 do $
        for i = 0, n_elements(uts)-1 do $
            plot_hope_l3_pitch2d_polygon, uts[i], types[j], unit = unit, $
                log = log, hopel3 = hopel3, probe = probes[k]
end