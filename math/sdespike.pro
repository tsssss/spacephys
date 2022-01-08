;+
; Type: procedure.
; Purpose: Remove spikes in given m-dimension line plots. 
; Parameters:
;   f0, in/out, dblarr[n,m]/dblarr[n]/string, req. Data, string for tplot var.
;   t0, in/out, dblarr[n], opt. Time or x data. If f0 is tplot var, then t0 
;       will be loaded, otherwise t0 can be set or generated with indgen.
;       t0 must be in uniform time.
; Keywords:
;   width, in, int, opt. Width to enclose spike for sure, default is 20.
;   nsigma, in, float, opt. Set threshold to select large deriv due to spikes.
; Notes: The width is usually good enough so you do not need to set. Do not set
;   it too short. The nsigma is set to 5 by default, lowing it will fix more
;   subtle spikes at the risk of including large amplitude and high freq wave.
; Dependence: slib.
; History:
;   2014-10-01, Sheng Tian, create.
;   2015-02-25, Sheng Tian, update algorithm.
;-

pro sdespike, t0, f0, width = width, nsigma = nsgma, ratio = ratio

    ; prepare t0 and f0.
    if n_params() eq 1 then begin
        if size(t0,/type) eq 7 then begin   ; tplot var.
            vname = t0 & get_data, vname, t0, f0
        endif else f0 = t0                  ; set f0 not t0.
    endif

    tmp = size(f0) & if tmp[0] gt 2 then message, 'more than 2 dimensions ...'
    ndim = (tmp[0] eq 1)? 1: tmp[2]
    nrec = n_elements(f0)/ndim
    if n_elements(vname) eq 0 and n_params() eq 1 then t0 = findgen(nrec)

    ; emperical settings.
    if n_elements(width) eq 0 then width = nrec/20  ; changed from 1000, Sheng, 2019-07-14.
    if n_elements(nsgma) eq 0 then nsgma = 10
    if n_elements(ratio) eq 0 then ratio = 0.6  ; use 60% mediocre points.
    ratio1 = (1-ratio)*0.5
    nsgma1 = nsgma
    nsgma2 = 10                 ; tell large deriv, larger to loose condition.
    wd = fix(6d/(t0[1]-t0[0]))  ; a padding width, 1 spin period of data.
    weight = 0.5*sin(smkarthm(0,!dpi,nrec,'n'))+1d  ; less picky in the middle.

    ; work on each dimension.
    for i = 0, ndim-1 do begin        
        ; original data.
        tf0 = f0[*,i]
        
        ; prepare a corrected deriv.
        df0 = [0d,tf0[1:nrec-1]-tf0[0:nrec-2]]      ; derivative or data.
        ; moving stddev for deriv.
        dfstddev = smvstddev(df0, t0, width = width, ratio = ratio1, /quadratic)
        ; remove spikes in deriv, which deviates a lot from the moving stddev.
        ; it works b/c deriv fluctuates around 0.
        idx = where(abs(df0) le nsgma1*dfstddev)
        df1 = interpol(df0[idx], t0[idx], t0)       ; the corrected deriv.
        deldf = df0-df1
        idx = where(deldf ne 0, cnt)    ; points of large derivative.
        if cnt eq 0 then break          ; no large derivative, so no spike.
        
        ; deviation from moving mean.
        del = tf0-smvmean(tf0, t0, width = width, ratio = ratio1, /quadratic)
        ; deviation significantly and span a large range.
        delstddev = smvstddev(del, t0, width = width, ratio = ratio1)
        
        ; idx saves potential bad points that deviate a lot from moving mean.
        idx = where(abs(del) gt nsgma2*delstddev)
        sp1 = 0ull     ; sp1/sp2 is the start/end index of idx.
        while sp1 lt n_elements(idx) do begin
            sp2 = sp1   ; find continous chunk of potential bad points.
            while sp2 lt n_elements(idx)-1 do $
                if idx[sp2+1]-idx[sp2] eq 1 then sp2 = sp2+1 else break
            i1 = idx[sp1] & i2 = idx[sp2]   ; start/end index of the chunk.
            if i2-i1 lt 2 then begin
                sp1 = sp2+1
                continue     ; the chunk is too short.
            endif
            ; only those points coincides with large derivatives are bad.
            ; use moving stddev of deriv, b/c wave has large deriv.
            tmp = stddev(df1[i1-1>0:i2+1<nrec-1])   ; local stddev of deriv.
            tmp*= weight[i1-1>0:i2+1<nrec-1]    ; omit large deriv in middle.
            tidx = where(abs(deldf[i1-1>0:i2+1<nrec-1]) gt nsgma2*tmp, cnt)

;red = sgcolor('red')
;!p.multi = [0,1,4] & !x.style = 1
;tpad = 10*wd
;ti1 = (i1-tpad)>0 & ti2 = (i2+tpad)<nrec-1 ; start/end index for bad points.
;plot, tf0, ytitle = 'orig data (full)'
;oplot, i1*[1,1], !y.crange, color = red
;oplot, i2*[1,1], !y.crange, color = red
;plot, (deldf)[ti1:ti2], ytitle = 'bad deriv'
;plot, tf0[ti1:ti2], ytitle = 'orig data'
;oplot, (i1-ti1)*[1,1], !y.crange, color = red
;oplot, (i2-ti1)*[1,1], !y.crange, color = red
;plot, del[ti1:ti2], ytitle = 'deviation', $
;    yrange = [-1,1]*3*nsgma2*max(delstddev[ti1:ti2])
;oplot, (nsgma2*weight*delstddev)[ti1:ti2], color = red
;oplot,-(nsgma2*weight*delstddev)[ti1:ti2], color = red
;!p.multi = -1
;print, 'chunck index:', i1-10*wd, i2+10*wd, nrec, cnt
;stop

            sp1 = sp2+1            
            if cnt gt 0 then begin
                ti1 = i1-wd>0
                if ti1 eq i1 then continue  ; can't recover if spike at start.
                tmp = df1[ti1:i2]  ; the corresponding chunk of deriv.
                tmp[i1-ti1:i2-ti1] -= mean(tmp[i1-ti1:i2-ti1])  ; int to 0.
                ; recover by int using the modified deriv and init value.
                tmp[0] = tf0[ti1]
                for j = 1, n_elements(tmp)-1 do tmp[j] = tmp[j-1]+tmp[j]
                tf0[i1:i2] = tmp[i1-ti1:i2-ti1]
            endif
        endwhile
        
        ; remove point spike.
        sdespike_median, t0, tf0, width = 5
        
        f0[*,i] = tf0
    endfor
    
    if n_elements(vname) eq 1 then store_data, vname, t0, f0 $
    else if n_params() eq 1 then t0 = f0
end
