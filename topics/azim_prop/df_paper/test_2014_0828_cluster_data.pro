
pres = ['tha','thd','the','c1','c2','c3','c4']+'_'

vars = pres+'b_sm'
deg = 180d/!dpi
foreach pre, pres do begin
    cotrans, pre+'b_gsm', pre+'b_sm', /gsm2sm
    get_data, pre+'b_sm', times, bsm
    bz = bsm[*,2]
    bxy = snorm(bsm[*,0:1])
    angle = atan(bz,bxy)*deg
    
    var = pre+'b_tilt'
    store_data, var, times, angle
    add_setting, var, /smart, {$
        display_type: 'scalar', $
        short_name: 'Tilt', $
        unit: 'deg'}
        
    var = pre+'bmag'
    store_data, var, times, snorm(bsm)
    add_setting, var, /smart, {$
        display_type: 'scalar', $
        short_name: '|B|', $
        unit: 'nT'}
endforeach

foreach var, ['c2','c3','c4']+'_b_tilt' do interp_time, var, to='c1_b_tilt'
stplot_merge, ['c1','c2','c3','c4']+'_b_tilt', newname='c_b_tilt', $
    colors=sgcolor(['red','green','blue','black']), labels='c'+['1','2','3','4']
foreach var, ['c2','c3','c4']+'_bmag' do interp_time, var, to='c1_bmag'
stplot_merge, ['c1','c2','c3','c4']+'_bmag', newname='c_bmag', $
    colors=sgcolor(['red','green','blue','black']), labels='c'+['1','2','3','4']

end