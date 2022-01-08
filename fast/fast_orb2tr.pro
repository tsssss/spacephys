;+
; Type: function.
; Purpose: Return time range of given fast orbit.
; Parameters:
;   orb, in, long, req. Given orbit. Valid value is [0,51323].
;       Valid meaningful orbit range is [330,51315].
; Keywords:
;   orbfile, in, string, opt. Specify the definitive orbit file.
;       '<root>/data/fast/orbit/definitive'.
; Return:
;   return, dblarr[2]. Time range in epoch.
; Notes: none.
; Dependence: slib.
; History:
;   2013-03-23, Sheng Tian, create.
;-
function line2time, line
    reads, line, yr, doy, hr, mi, sc, $
        format = '(19x,I4,1x,I3,2(1x,I2),1x,F)'
    tmp = sfmdoy(yr, doy)
    mo = tmp[0]
    dy = tmp[1]
    cdf_epoch, epoch, yr, mo, dy, hr, mi, sc, /compute_epoch
    return, sfmepoch(epoch, 'unix')
end

function fast_orb2tr, orb, orbfile = fn

    compile_opt idl2

    ; check orbit.
    if orb gt 51315 or orb lt 0 then message, 'wrong orbit ...'

    ; check file name.
    if n_elements(fn) eq 0 then $
        fn = sdiskdir('Research')+'/data/fast_cdf/orbit/definitive'

    ; open file.
    openr, lun, fn, /get_lun
    line = ''
    skip_lun, lun, 7*orb+2, /lines  ; 7 lines per orbit, and 2 lines in header.
    readf, lun, line
    t1 = line2time(line)
    skip_lun, lun, 6, /lines
    readf, lun, line
    t2 = line2time(line)
    free_lun, lun
    return, [t1,t2]

end
