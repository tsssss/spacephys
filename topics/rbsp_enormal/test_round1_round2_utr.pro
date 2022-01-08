; round1: survey data, 11 Hz.
; round2: high res data, 32 Hz.
; check the time ranges.


id1s = psbl_de_id('spin_res')
id2s = psbl_de_id('detect_step')

utr0 = time_string(['2013-01-01','2016-01-01'])
store_data, 'test', utr0, [0,0]

ofn = shomedir()+'/test_round1_round2_overlap.pdf'
sgopen, ofn, xsize = 8, ysize = 3, /inch

tplot, 'test', trange = utr0, title = $
    'Event distribution, blue: round1 (11 Hz), red: round2 (32 Hz)'
timebar, id1s.utr, color = sgcolor('blue')
timebar, id2s.utr, color = sgcolor('red')

sgclose

end