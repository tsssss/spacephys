;+
;-

function sreadtxt, ifile, start, nlines, ofile = ofile

  compile_opt idl2
  
  on_error, 2
  
  if ~file_test(ifile) then message, 'in file does not exist ...'
  
  totallines = file_lines(ifile)
  
  if n_elements(start) eq 0 then start = 0
  
  if start lt 0 then message, 'start line less than 0 ...'
 
  if n_elements(nlines) eq 0 then nlines = totallines - start

  file = strarr(nlines)
  
  openr, ilun, ifile, /get_lun
  
  if start gt 1 then begin
    header = strarr(start-1)
    readf, ilun, header
  endif
    
  if nlines le 0 then message, 'no text ...'
  
  tline = ''
  for ii = 0L, nlines - 1L do begin
    readf, ilun, tline
    file[ii] = tline
  endfor
  
  free_lun, ilun
  
  if n_elements(ofile) ne 0 then begin
    openw, olun, ofile, /get_lun
    for ii = 0L, nlines - 1L do printf, olun, file[ii]
    free_lun, olun
  endif
  
  return, file
  
end
