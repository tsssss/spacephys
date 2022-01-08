;+
; Check bad B0 data in some bins.
;-

probe = 'a'
bin_info_name = 'spherical_spatial'
prefix = 'rbsp'+probe+'_'

var_type = 'de_sm'
mlt_range = [8]+[-1,1]*0.5
mlat_range = -20+[0,5]
dis_range = 5.5+[0,0.5]

;; Check an empty bin.
;mlt_range = [-10]+[-1,1]*0.5
;mlat_range = -20+[0,5]
;dis_range = 4.0+[0,0.5]


given_ranges = list()
given_ranges.add, mlt_range
given_ranges.add, mlat_range
given_ranges.add, dis_range


bin_info = pflux_survey_load_bin_info(bin_info_name, probe=probe)
index_list = list()
foreach var, bin_info.vars, var_id do begin
    foreach bin, bin_info[var].bins do begin
        if product(bin.range-given_ranges[var_id]) ne 0 then continue
        index_list.add, bin.index
    endforeach
endforeach

pflux_survey_load_data, var_type, probe=probe
get_data, prefix+var_type, times, data

;bad_time_ranges = list()
;bad_time_ranges.add, time_double(['2014-08-20','2014-08-28'])
;bad_time_ranges.add, time_double(['2015-09-14','2015-09-17'])
;foreach bad_time_range, bad_time_ranges do begin
;    index = lazy_where(times, '[]', bad_time_range)
;    data[index,*] = !values.f_nan
;endforeach

time_index = lindgen(n_elements(times))
foreach index, index_list do time_index = set_inter(time_index, index)
times = times[time_index]
data = data[time_index,*]
print, mean(snorm(data))

end
