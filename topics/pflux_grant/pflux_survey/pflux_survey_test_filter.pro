
probe = 'b'
prefix = 'rbsp'+probe+'_'

vars = ['mlat','mlt','dis','pf_fac_norm']
foreach var, vars do pflux_survey_load_data, var, probe=probe

var = 'mlt'
range = 6+[-1,1]*0.5
mapping_index1 = pflux_survey_filter_data(var, range, probe=probe)

var = 'mlat'
range = [10,15]
range = [15,20]
mapping_index2 = pflux_survey_filter_data(var, range, probe=probe)

max_rec = max([max(mapping_index1),max(mapping_index2)])
mapping_index = bytarr(max_rec)+1
tmp_index = bytarr(max_rec)
tmp_index[mapping_index1] = 1
mapping_index = mapping_index and tmp_index
tmp_index = bytarr(max_rec)
tmp_index[mapping_index2] = 1
mapping_index = mapping_index and tmp_index

index = where(mapping_index eq 1, count)
var = prefix+'pf_fac_norm'
get_data, var, times, data
times = times[index]
data = data[index,*]

stop
end