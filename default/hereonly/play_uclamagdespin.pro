
; I test whether the returned magnetic field is different 
; if I feed ucla_mag_despin with {Mag1dc_S, Mag2dc_S, Mag3dc_S}
; or {MagXDC, MagYDC, MagZDC} or {MagXDC, MagYDC, MagZDC}.
;
; test data is orbit 3020, interested time range is 
; trange = time_double(['1997-05-27/18:10','1997-05-27/18:25']))


pro play_uclamagdespin, inm

  print, 'Load 3 types of B field and compare output of ucla_mag_despin.'
  print, 'Use FAST orbit 03020.'
  print, 'Time range: 1997-05-27 18:10 UT to 1997-05-27 18:25UT.'

  nmmag = '_Mag'
  nmmagdc = '_MagDC'
  nmmagdcs = '_MagDC_S'
  trange = time_double(['1997-05-27/18:10','1997-05-27/18:25'])
  if n_elements(inm) eq 0 then $
    inm = 'dB_fac_v'

  ; in sdt, open a window and add plot:
  ; Fields-Survey -> {Mag1dc_S, Mag2dc_S, Mag3dc_S},
  ; AttitudeCtrl -> SUN.
  print, 'Load Mag1dc_S, Mag2dc_S, Mag3dc_S, SUN in SDT'
  stop
  ucla_mag_despin
  get_data, inm, data = dbmagdcs, limits = lmmagdcs
  store_data, inm+nmmagdcs, data = dbmagdcs, limit = lmmagdcs
  options, inm+nmmagdcs, 'ytitle', 'Magdc_S!CdB_fac_v!C(nT)'

  ; in sdt, open a window and add plot:
  ; Fields-Survey -> {MagX, MagY, MagZ},
  ; AttitudeCtrl -> SUN.
  print, 'Load MagX, MagY, MagZ, SUN in SDT'
  stop
  ucla_mag_despin
  get_data, inm, data = dbmag, limits = lmmag
  store_data, inm+nmmag, data = dbmag, limits = lmmag
  options, inm+nmmag, 'ytitle', 'Mag!CdB_fac_v!C(nT)'

  ; in sdt, open a window and add plot:
  ; Fields-Survey -> {MagXDC, MagYDC, MagZDC},
  ; AttitudeCtrl -> SUN.
  print, 'Load MagXDC, MagYDC, MagZDC, SUN in SDT'
  stop
  ucla_mag_despin
  get_data, inm, data = dbmagdc, limits = lmmagdc
  store_data, inm+nmmagdc, data = dbmagdc, limits = lmmagdc
  options, inm+nmmagdc, 'ytitle', 'Magdc!CdB_fac_v!C(nT)'

  loadct3, 45
  device, decomposed = 0
  title = 'UCLA_Mag_Despin 3 input B field comparison'
  tplot, [inm+'_*'], trange = trange, title = title
  pstplot, [inm+'_*'], trange = trange, title = title, $
    filename = '~/tmp/backup_ucla_mag_despin_input_Bfield_compare.ps'
  tplot_save, [inm+'_*'], filename = '~/tmp/backup_ucla_mag_despin_input_Bfield_compare'

end
