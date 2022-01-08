
yr = 1990
mo = 05
dy = 19
hr = 1
mi = 3
sc = 40
ms = 555

t1 = systime(/sec)
for ii = 0, 10000 do $
  cdf_epoch, epoch, yr, mo, dy, hr, mi, sc, ms, /compute
t2 = systime(/sec)
print, t2-t1, format = '(F20)'

t1 = systime(/sec)
for ii = 0, 1000 do $
 jd = julday(mo, dy, yr, hr, mi, sc+ms/1000D)
t2 = systime(/sec)
print, t2-t1, format = '(F20)'

end