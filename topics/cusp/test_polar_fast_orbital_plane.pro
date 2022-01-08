;+
; This code studies the precession speed of Polar and FAST orbital plane, 
; Polar is 360 deg/year, FAST is 1.5*360 deg/year.
;-
pro test_polar_fast_orbital_plane

    ; [date, polar plane, fast plane]
    dat = [$
        time_double('1997-01-01'), 08, 01, $
        time_double('1997-02-01'), 05, 10, $
        time_double('1997-03-01'), 03, 07, $
        time_double('1997-04-01'), 01, 04, $
        time_double('1997-05-01'), 11, 01, $
        time_double('1997-06-01'), 09, 11, $
        time_double('1997-07-01'), 06, 07, $
        time_double('1997-08-01'), 05, 04, $
        time_double('1997-09-01'), 03, 01, $
        time_double('1997-10-01'), 01, 11, $
        time_double('1997-11-01'), 11, 08, $
        time_double('1997-12-01'), 09, 05, $
        
        time_double('1998-01-01'), 19.0, 13.5, $
        time_double('1999-01-01'), 17.5, 14.5, $
        time_double('2000-01-01'), 18.5, 15.0, $
        time_double('2001-01-01'), 18.0, 15.0, $
        time_double('2002-01-01'), 18.0, 15.0, $
        time_double('2003-01-01'), 17.5, 16.0, $
        time_double('2004-01-01'), 17.0, 16.0, $
        time_double('2005-01-01'), 15.5, 16.0]
        
        
    ; 1997 to 2008.
    nyr = 12 & yrs = 1997+findgen(nyr)
    nmo = 12 & mos = 1+findgen(nmo)
    
    mltp = 8D & dtp = 10D   ;-2D
    mltf = 1D & dtf = 09D   ;-3D
    
    for i = 0, nyr-1 do begin
        for j = 0, nmo-1 do begin
            tp = (mltp+(i*12+j)*dtp) mod 12
            tf = (mltf+(i*12+j)*dtf) mod 12
            print, string(yrs[i],format='(I4)')+'-'+string(mos[j],format='(I02)'), $
                '    '+string(tp,format='(I2)')+'    '+string(tf,format='(I2)')
        endfor
    endfor

end
