
gs_exe = 'gs'
testcmd = gs_exe+' -version'
spawn, testcmd, result

cmd = gs_exe + ' -sDEVICE=pdfwrite -q -dNOPAUSE -dBATCH'+ $
  ' -sPAPERSIZE=' + StrLowCase(pagetype) + $
  ' -sOutputFile="' + pdf_file + '" "' + ps_file + '"'
  
**** ps2pdf can be used in unix system.