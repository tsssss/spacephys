

logfn = shomedir()+'/Google Drive/works/works/rbsp_de/dipolarization/list_large_de_round3.log'
nheader = 3
headers = strarr(nheader)
nline = file_lines(logfn)-nheader
lines = strarr(nline)
openr, lun, logfn, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun

for i = 74, nline-1 do begin
    store_data, '*', /delete
    
    tfn = lines[i]
    id = strmid(tfn,0,16)
    tprobe = strmid(id,strlen(id)-1)
    
    t0 = time_double(strmid(id,0,9)+strmid(tfn,20,5),tformat='YYYY_MMDDhh:mm')
    t1 = time_double(strmid(id,0,9)+strmid(tfn,29,5),tformat='YYYY_MMDDhh:mm')
    if t1 lt t0 then t1+= 86400
    
    utr = [t0,t1]
    print, time_string(utr)
    de_at_gradb_survey_plot2, utr, tprobe, id = id

endfor

end