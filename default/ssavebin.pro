
pro ssavebin, fn, var

    if n_elements(fn) eq 0 then message, 'no file name ...'
    if file_test(fn) eq 0 then stouch, fn
    openw, lun, fn, /append, /get_lun
    writeu, lun, var
    free_lun, lun

end

nrec = 100
fn = shomedir()+'/test_sine_wave.bin'
txs = findgen(nrec+1)/nrec*2*!dpi
tys = sin(txs)

ssavebin, fn, txs
ssavebin, fn, tys

end
