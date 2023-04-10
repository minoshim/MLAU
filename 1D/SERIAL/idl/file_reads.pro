function file_reads,filenames

files=findfile(filenames)
nn=n_elements(files)

for n=0,nn-1 do begin
   tmp=read_ascii(files[n])
   if (n eq 0) then str=tmp else str=[str,tmp]
endfor

return,str.(0)
end
