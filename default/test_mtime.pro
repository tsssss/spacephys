
pro test_mtime

    case !version.os_family of
        'unix': begin       ; from ut, can get ut, then convert to lt, mtime input is lt.
            fn = shomedir()+'/test.txt'
            spawn, 'touch '+fn
            finfo1 = file_info(fn)
            print, '**** file created on'
            print, 'local time: ', systime()
            print, 'univs time: ', time_string(double(finfo1.mtime))
            
            mtime0 = 'Thu, 02 May 2013 01:24:35 GMT'
            print, '**** change mtime to: ', mtime0
            
            spawn, 'date -ujf "%a, %d %b %Y %H:%M:%S GMT" "'+mtime0+$
                '" +"%s"', mtime1
            mtime1 = double(mtime1[0])
            print, 'univs time: ', time_string(mtime1)
            
            spawn, 'date +"%z"', dt
            dt = double(dt)*36d
            mtime1 += dt
            print, 'local time: ', time_string(mtime1)
            
            mtime2 = time_string(mtime1,tformat='YYYYMMDDhhmm.ss')
            spawn, 'touch -cmt '+mtime2+' '+fn
            finfo2 = file_info(fn)
            print, '**** current mtime: ', time_string(double(finfo2.mtime))
            end
        'Windows': begin    ; from ut, can get lt directly, mtime input is lt.
            fn = shomedir()+'\test.txt'
            spawn, 'powershell New-Item -name '+file_basename(fn)+' -path '+file_dirname(fn)+' -type "file"', tmp, /noshell
            finfo1 = file_info(fn)
            print, '**** file created on'
            print, 'local time: ', systime()
            print, 'univs time: ', time_string(double(finfo1.mtime))
            
            mtime0 = 'Thu, 02 May 2013 01:24:35 GMT'    ; 5 hr for summer.
            mtime0 = 'Sun, 01 Feb 2015 08:21:33 GMT'    ; 6 hr for winter.
            print, '**** change mtime to: ', mtime0
            
            spawn, 'powershell $a = Get-Date -Date \"'+mtime0+'\" -UFormat %s;$a', mtime1
            mtime1 = double(mtime1[0])
            print, 'local time: ', time_string(mtime1)
            
            spawn, 'powershell $a = Get-Date -Date \"'+mtime0+'\" -UFormat %Z;$a', dt
            dt = double(dt)*3600d
            mtime2 = mtime1-dt
            print, 'univs time: ', time_string(mtime2)
            
            mtime3 = time_string(mtime2,tformat='YYYY-MM-DD hh:mm:ss')
            spawn, 'powershell $a = Get-Date -Date \"'+mtime3+'\" ;(ls '+fn+').LastWriteTimeUtc = $a', /noshell
            finfo2 = file_info(fn)
;            print, '**** current mtime: ', time_string(double(finfo2.mtime))
            spawn, 'powershell (ls '+fn+').LastWriteTimeUtc', tmp, /noshell
            print, '**** current mtime: ', tmp
            end
    endcase

end