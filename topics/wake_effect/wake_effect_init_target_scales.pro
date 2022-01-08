;+
; Return a dictionary of target scales.
;-
function wake_effect_init_target_scales, spin_rate

    spin_freq = 1d/spin_rate
    base = 2d^(1d/4)
    freq_range = [0.5,8]/spin_rate
    selected_freqs = smkgmtrc(min(freq_range),max(freq_range),base,'dx')
    nselected_freq = n_elements(selected_freqs)
    selected_periods = 1d/selected_freqs
    wavelet_info = wavelet_info()
    t2s = wavelet_info[7]
    selected_scales = selected_periods*t2s

    return, dictionary($
        'base', base, $                 ; separation between terms.
        'spin_rate', spin_rate, $       ; spin period in sec.
        'spin_freq', spin_freq, $       ; spin freq in Hz.
        'nscale', nselected_freq, $     ; # of target scales.
        'scales', selected_scales, $    ; scales in sec.
        'freqs', selected_freqs, $      ; freqs in Hz.
        'periods', selected_periods, $  ; periods in sec.
        't2s', t2s, $                   ; convert period to scale.
        's2t', 1d/t2s)                  ; convert scale to period.
end
