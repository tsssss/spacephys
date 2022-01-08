
; before use, start sdt. 
;
; for fa_fields_despin
; mode	name			DQD
; D/U	V5-V8_S			V5-V8_S
; D	V1-V2_S			V1-V2_S
; U	V1-V4_S			V1-V4_S
; U	V4_S			V4_S
; U	V8_S			V8_S
; D/U	1032-spinPhase		SMPhase_FieldsSurvey0	
; 
; For ucla_mag_despin
; mode	name		DQD
; D	MagX		MagXYZ
; D	MagY		MagXYZ
; D	MagZ		MagXYZ
; D	Sun		AttitudeCtrl


pro play_fafieldsdespin
 
  vnms = ['_V12','_V158']
  ebnm = 'E_NEAR_B'
  evnm = 'E_ALONG_V'
 
  fa_fields_despin, /use_v158
  get_data, ebnm, data = eb_v158
  get_data, evnm, data = ev_v158

  fa_fields_despin

  store_data, ebnm+vnms[1], data = eb_v158
  store_data, evnm+vnms[1], data = ev_v158

end
