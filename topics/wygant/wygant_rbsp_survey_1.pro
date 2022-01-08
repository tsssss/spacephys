;+
; Hope moments, E/B fields.
; Density (electron and all ions).
; Electron differential energy spectrogram.
; Magnetic field (|B|, Bxyz GSM).
; Electric field (E0xyz GSM).
; Bulk velocity (GSM).
; Electron energy flux (GSM).
; All Ion energy flux (GSM).
;-
pro wygant_rbsp_survey_1, utr0, tprobe

    hope_mom = sread_rbsp_hope_moments(utr0, probe=tprobe)
    
    

end