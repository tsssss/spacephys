;+
; Type: procedure.
; Purpose: Plot time-frequency spectrogram.
; Parameters: <+++>.
; Keywords: <+++>.
; Notes: <+++>.
; Dependence: slib.
; Author: Sheng Tian.
; History: 2014-01-11, Sheng Tian, create.
;-

pro splot_spect, vname
    compile_opt idl2

    get_data, vname, tmp = tmp
    xx = tmp.x      ; x-axis coord.
    yy = tmp.v      ; scales or Fourier frequency.
    zz = tmp.y      ; spectrogram.
    cc = tmp.c      ; cone of influence.
end
