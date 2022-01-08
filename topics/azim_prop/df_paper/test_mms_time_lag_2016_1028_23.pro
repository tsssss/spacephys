;+
; Test MMS time lag.
;-
;

event_id = '2016_1028_23'
time_range = time_double(['2016-10-28/22:30','2016-10-29/01:00'])
model = 't89'
deg = 180d/!dpi
re = 6378d
re1 = 1d/re
xyz = ['x','y','z']
rgb = sgcolor(['red','green','blue'])

dir = join_path([shomedir(),'Downloads','mms'])
files = file_search(join_path([dir,'*.cdf']))

data_file = join_path([shomedir(),'azim_prop','data',event_id+'_basic_data.tplot'])
tplot_restore, filename=data_file

foreach file, files do cdf2tplot, file
mms_probes = 'mms'+['1','2','3','4']
foreach probe, mms_probes do begin
    get_data, probe+'_fgm_r_gsm_srvy_l2', times, rgsm & rgsm = rgsm[*,0:2]*re1
    rsm = cotran(rgsm, times, 'gsm2sm')
    store_data, probe+'_r_sm', times, rsm
    add_setting, probe+'_r_sm', /smart, {$
        display_type: 'vector', $
        unit: 'Re', $
        short_name: 'R', $
        coord: 'SM', $
        coord_labels: xyz, $
        colors: rgb}
    store_data, probe+'_r_gsm', times, rgsm
    add_setting, probe+'_r_gsm', /smart, {$
        display_type: 'vector', $
        unit: 'Re', $
        short_name: 'R', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}

;---B model GSM.
    ntime = n_elements(times)
    bmodgsm = fltarr(ntime,3)
    for kk=0,ntime-1 do begin
        tilt = geopack_recalc(times[kk])
    
        ; in-situ position
        rx = rgsm[kk,0]
        ry = rgsm[kk,1]
        rz = rgsm[kk,2]
    
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, 2, rx,ry,rz, dbx,dby,dbz
        bmodgsm[kk,*] = [bx,by,bz]+[dbx,dby,dbz]
    endfor
    store_data, probe+'_bmod_gsm_'+model, times, bmodgsm
    add_setting, probe+'_bmod_gsm_'+model, /smart, {$
        display_type: 'vector', $
        unit: 'Re', $
        short_name: strupcase(model)+' B', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}
    
    
;---B GSM, SM.
    get_data, probe+'_fgm_b_gsm_srvy_l2_clean', times, bgsm & bgsm = bgsm[*,0:2]
    bsm = cotran(bgsm, times, 'gsm2sm')
    store_data, probe+'_b_sm', times, bsm
    add_setting, probe+'_b_sm', /smart, {$
        display_type: 'vector', $
        unit: 'Re', $
        short_name: 'B', $
        coord: 'SM', $
        coord_labels: xyz, $
        colors: rgb}
    store_data, probe+'_b_gsm', times, bgsm
    add_setting, probe+'_b_gsm', /smart, {$
        display_type: 'vector', $
        unit: 'Re', $
        short_name: 'B', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}
        
        
    b_tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
    ;sdespike, times, b_tilt
    bmodgsm = get_var_data(probe+'_bmod_gsm_'+model, at=times)
    bmodsm = cotran(bmodgsm, times, 'gsm2sm')
    bmod_tilt = atan(bmodsm[*,2],sqrt(bmodsm[*,0]^2+bmodsm[*,1]^2))*deg
    
    dis = snorm(get_var_data(probe+'_r_gsm', at=times))
    db_tilt = b_tilt-bmod_tilt
    index = where(dis ge 4 and db_tilt gt -100)
    yrange = minmax(db_tilt[index])
    store_data, probe+'_db_tilt', times, db_tilt, limits={$
        ytitle: '(deg)', labels:strupcase(probe), yrange: yrange}
endforeach


end