

logfn = shomedir()+'/Google Drive/works/works/rbsp_de/de_at_gradb/list_de_at_gradb.log'
nheader = 3
headers = strarr(nheader)
nline = file_lines(logfn)-nheader
lines = strarr(nline)
openr, lun, logfn, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun

for i = 0, nline-1 do begin
    store_data, '*', /delete
    
    tfn = file_basename(lines[i])
    tprobe = strlowcase(strmid(tfn,19,1))
    t0 = strmid(tfn,0,14)
    ut0 = time_double(t0, tformat='YYYY_MMDD_hhmm')
    
    utr = ut0+[-0.5,1.5]*60
    de_at_gradb_survey_plot2, utr, tprobe, id = t0
    
    utr = ut0+[-1,1]*15*60
    de_at_gradb_survey_plot3, utr, tprobe, id = t0
endfor

end