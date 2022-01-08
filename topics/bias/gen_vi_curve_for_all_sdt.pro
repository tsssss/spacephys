
; all known sdt, newest first.
; Apr 02, 2014    ; this one didn't run.
dates = [$
;    'Nov 14, 2015', $
;    'Oct 14, 2015', $
;    'Jul 22, 2015', $
;    'Jan 10, 2015', $
;    'Dec 17, 2014', $
;    'Dec 03, 2014', $
;    'Jun 18, 2014', $
    'May 30, 2014', $
;    'Oct 30, 2014', $
;    'Sep 11, 2013', $
;    'Aug 07, 2013', $
;    'Jul 31, 2013', $
;    'Jul 17, 2013', $
;    'Jul 03, 2013', $
;    'May 08, 2013', $
;    'Apr 17, 2013', $
;    'Mar 27, 2013', $
;    'Feb 27, 2013', $
;    'Feb 06, 2013', $
    'Nov 16, 2013', $   ; 2012?
    'Dec 19, 2012']
dates = ['Aug 17, 2016']
ndate = n_elements(dates)
uts = dblarr(ndate)
for i = 0, ndate-1 do uts[i] = sfmdate(dates[i],'%b %d, %Y')
for i = 0, ndate-1 do dates[i] = stodate(uts[i],'%Y-%m-%d')
for i = 0, ndate-1 do begin
    plot_vi_curve, dates[i], probes = ['a','b'], /reload
endfor

end