AGNinterp_thinned = function(lum = 1e+44, ct = 60, al = 4, be = -0.5, ta = 1, rm = 60,
                     an = 30, em_factor = 1, Fritz = NULL, Temple_emline_reduced = NULL){
  
if(is.null(Fritz)){
  data('Fritz', envir = environment())
}

if(is.null(Temple_emline_reduced)){
  data('Temple_emline_reduced', envir = environment())
}

ctmix = interp_quick(ct, Fritz$ct)
almix = interp_quick(al, Fritz$al)
bemix = interp_quick(be, Fritz$be)
tamix = interp_quick(ta, Fritz$ta)
rmmix = interp_quick(rm, Fritz$rm)
anmix = interp_quick(an, Fritz$an)


slice = Fritz$Aspec[c(ctmix[1:2]), c(almix[1:2]), c(bemix[1:2]),
                           c(tamix[1:2]), c(rmmix[1:2]), c(anmix[1:2]), ]

weights = rep(1, 64)
weights = weights * rep(ctmix[3:4], each = 1, times = 32)
weights = weights * rep(almix[3:4], each = 2, times = 16)
weights = weights * rep(bemix[3:4], each = 4, times = 8)
weights = weights * rep(tamix[3:4], each = 8, times = 4)
weights = weights * rep(rmmix[3:4], each = 16, times = 2)
weights = weights * rep(anmix[3:4], each = 32, times = 1)

tempmat = matrix(as.numeric(slice), 64, 178)
agn_spectrum = lum * colSums(tempmat * weights)

Temple_emline_reduced$flux = Temple_emline_reduced$flux * agn_spectrum[44]
emission_spectrum = em_factor * Temple_emline_reduced$flux

combined_spectrum = agn_spectrum
combined_spectrum[41:93] = agn_spectrum[41:93] + emission_spectrum

return(data.frame(wave = Fritz$Wave, lum = combined_spectrum))
}
