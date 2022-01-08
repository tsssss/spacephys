;+
; Type:	procedure.
; Purpose: Print !path in lines to console.
; Parameters: none.
; Keywords: none.
; Notes: none.
; Dependence: none.
; History:
; 	2012-05-17, Sheng Tian, create.
;	2013-03-19, Sheng Tian, re-document.
;-

pro slprpath, path0

    if n_elements(path0) eq 0 then path0 = !path

	sep = path_sep(/search_path)

	paths = strsplit(path0, sep, /extract)
	npath = n_elements(paths)
	print, 'number of paths: ', npath

	for ii = 0, npath-1 do print, paths[ii]
end