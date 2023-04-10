pro binary_read,array,filename=filename,_extra=extra

if not(keyword_set(filename)) then begin
   filename=''
   filename = dialog_pickfile(filter=filename)
endif
if filename eq '' then return

if file_search(filename) ne '' then begin
   openr,unit,filename,/get_lun,_extra=extra
   readu,unit,array
   close,unit
   free_lun,unit
endif else message,'No file: '+filename,/cont

end

