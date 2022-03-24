function interpol_localized, y, x, xnew
 if n_elements(x) eq 2 then return, interpol(y, x, xnew) $
 else if n_elements(x) eq 3 then return, interpol(y, x, xnew, /quadratic) $
 else return, interpol(y, x, xnew, /spline)
end