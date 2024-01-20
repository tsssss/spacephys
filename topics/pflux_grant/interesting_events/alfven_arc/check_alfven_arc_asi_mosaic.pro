
time_range = time_double(['2015-01-05/00:00','2015-01-05/02:00'])
sites = ['nrsq','gill','snkq','rank']
min_elevs = float([2.5,2.5,2.5,2.5])
merge_method = 'merge_elev'
calibration_method = 'moon'
zlog = 1
zrange = [5e2,1e4]


time_range = time_double(['2015-02-18/10:30','2015-02-18/11:05'])
sites = ['inuv']
min_elevs = float([0])
merge_method = 'max_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]

time_range = time_double(['2015-03-01/09:30','2015-03-01/10:50'])
sites = ['mcgr','fykn','atha']
min_elevs = float([1,1,1])
merge_method = 'merge_elev'
calibration_method = 'moon'
zlog = 1
zrange = [1e3,1e4]



time_range = time_double(['2015-03-02/10:30','2015-03-02/11:30'])
sites = ['mcgr','fykn','inuv']
min_elevs = float([1,1,1])
merge_method = 'merge_elev'
calibration_method = 'moon'
zlog = 1
zrange = [5e2,1e4]



time_range = time_double(['2018-10-09/10:00','2018-10-09/11:00'])
sites = ['kian','gako']
min_elevs = float([2.5,2.5])
resolutions = ['ast','asf']
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]


time_range = time_double(['2015-02-17/09:40','2015-02-17/10:40'])
sites = ['fykn','atha']
min_elevs = float([0,2.5])
resolutions = ['ast','asf']
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]

time_range = time_double(['2018-11-05/05:00','2018-11-05/07:30'])
sites = ['kian','mcgr','gako','fsim','pina','snkq','nrsq']
min_elevs = fltarr(n_elements(sites))+2.5
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]

time_range = time_double(['2018-09-13/11:00','2018-09-13/12:00'])
sites = ['gako','mcgr','inuv']
min_elevs = fltarr(n_elements(sites))+0
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]


time_range = time_double(['2018-09-22/11:00','2018-09-22/12:00'])
sites = ['kian','whit','fsim','inuv']
min_elevs = fltarr(n_elements(sites))+0
resolutions = ['ast','ast','asf','ast']
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]


time_range = time_double(['2018-10-08/23:00','2018-10-09/01:00'])
sites = ['nrsq']
min_elevs = float([1])
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]



time_range = time_double(['2015-03-17/08:40','2015-03-17/09:10'])
sites = ['whit','fykn','atha','fsmi','pina']
min_elevs = fltarr(n_elements(sites))+1
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
zrange = [1e2,1e4]



time_range = time_double(['2017-03-10/05:30','2017-03-10/06:00'])
time_range = time_double(['2017-03-10/09:30','2017-03-10/10:00'])
sites = ['fsim','fsmi','gill','kian','gako','whit']
min_elevs = [5,10,5,1]
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'moon_smooth'
zlog = 1
zrange = [1,1e4]

time_range = time_double(['2017-03-09/07:00','2017-03-09/08:00'])
sites = ['fsim','fsmi','gill','fykn']
min_elevs = [5,10,2.5,2.5]
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'moon'
zlog = 1
zrange = [5e2,1e4]


time_range = time_double(['2015-04-16/07:30','2015-04-16/09:00'])
sites = ['mcgr','whit']
min_elevs = [5,5]
resolutions = !null
merge_method = 'max_elev'
calibration_method = 'moon_smooth'
zlog = 0
zrange = [1e2,8e3]
zrange = [-1,1]*1e3

time_range = time_double(['2015-03-17/04:00','2015-03-17/09:00'])
time_range = time_double(['2015-03-17/06:00','2015-03-17/10:00'])
sites = ['inuv','whit','atha','fsim','fsmi','pina','kapu','snkq','gbay','nrsq']
min_elevs = [10,5,5,10,10,5,5,10,5,5]/2
resolutions = !null
merge_method = 'merge_elev'
calibration_method = 'simple'
zlog = 1
if zlog eq 0 then zrange = [-1,1]*4e3 else zrange = [5e1,1e4]
update = 0
;zrange = [1e1,1e4]


;time_range = time_double(['2015-02-17/09:00','2015-02-17/11:00'])
;sites = ['inuv','fsmi'];,'whit','fsmi','atha'];,'rank','snkq','kuuj']
;min_elev = !null
;merge_method = 'merge_elev'

mlt_image_var = themis_read_mlt_image(time_range, sites=sites, $
    min_elevs=min_elevs, resolutions=resolutions, update=update, $
    merge_method=merge_method, calibration_method=calibration_method)

get_data, mlt_image_var, times, mlt_images
index = where(finite(mlt_images,nan=1), count)
if count ne 0 then mlt_images[index] = 0
sgopen, 0, size=[1,1]*6

if zlog then begin
    zzs = alog10(mlt_images)
    zr = alog10(zrange)
endif else begin
    zzs = mlt_images
    zr = zrange
endelse

stop
angles = smkarthm(0,2*!dpi,40,'n')
foreach time, times, time_id do begin
    ct = 70
    if zlog eq 1 then ct = 49
    sgtv, bytscl(reform(zzs[time_id,*,*]),top=254, min=zr[0], max=zr[1]), position=[0,0,1,1], ct=ct
    plot, [-1,1],[-1,1], nodata=1, noerase=1, position=[0,0,1,1], xstyle=5, ystyle=5
    foreach rr, smkarthm(0,1,0.125,'dx') do plots, rr*cos(angles),rr*sin(angles), color=sgcolor('silver'), linestyle=1
    foreach rr, smkarthm(0,1,0.25,'dx') do plots, rr*cos(angles),rr*sin(angles), color=sgcolor('silver'), linestyle=0
    foreach tt, smkarthm(0,2*!dpi,25,'n') do plots, [0,1]*cos(tt),[0,1]*sin(tt), color=sgcolor('silver'), linestyle=1
    foreach tt, smkarthm(0,2*!dpi,13,'n') do plots, [0,1]*cos(tt),[0,1]*sin(tt), color=sgcolor('silver'), linestyle=0
    xyouts, 10,10, device=1, time_string(time)
    if time eq time_double('2015-04-16/08:03') then stop
    if time eq time_double('2015-04-16/08:06') then stop
    if time eq time_double('2017-03-09/07:45') then stop
endforeach


end