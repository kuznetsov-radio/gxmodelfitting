function InSav, o, name
 s=o->Names()
 res=0
 for i=0, n_elements(s)-1 do if strcmp(s[i], name, /fold_case) then res=1
 return, res
end 