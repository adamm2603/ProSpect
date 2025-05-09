AGNinterp_emlines_ = function(lum = 1e+44, ct = 60, al = 4, be = -0.5, ta = 1, rm = 60, an = 30, em_factor = 1, Fritz_interp = NULL, EmSpec = NULL)
{
    
if ((is.null(Fritz_interp) & em_factor != 0) | inherits(Fritz_interp, 'Fritz_interp')) {
Fritz_interp = NULL
data("Fritz_interp", envir = environment())
}

  
if (is.null(EmSpec)) {
# Load or define the emission spectrum if not provided
data("Temple_emline", envir = environment())
}


ctmix = interp_quick(ct, Fritz_interp$ct)
almix = interp_quick(al, Fritz_interp$al)
bemix = interp_quick(be, Fritz_interp$be)
tamix = interp_quick(ta, Fritz_interp$ta)
rmmix = interp_quick(rm, Fritz_interp$rm)
anmix = interp_quick(an, Fritz_interp$an)


slice = Fritz_interp$Aspec[c(ctmix[1:2]), c(almix[1:2]), c(bemix[1:2]),
            c(tamix[1:2]), c(rmmix[1:2]), c(anmix[1:2]), ]

weights = rep(1, 64)
weights = weights * rep(ctmix[3:4], each = 1, times = 32)
weights = weights * rep(almix[3:4], each = 2, times = 16)
weights = weights * rep(bemix[3:4], each = 4, times = 8)
weights = weights * rep(tamix[3:4], each = 8, times = 4)
weights = weights * rep(rmmix[3:4], each = 16, times = 2)
weights = weights * rep(anmix[3:4], each = 32, times = 1)


tempmat = matrix(as.numeric(slice), 64, 3236)

# Compute the AGN spectrum
agn_spectrum = lum * colSums(tempmat * weights)
#Scale up the emission lines to match the continuum (the emission line spectrum is
#normalised to be equal to 1 at the peak of the CIV emission line.
Temple_emline$flux = Temple_emline$flux * agn_spectrum[400]
# Apply the emission factor to the emission spectrum
emission_spectrum = em_factor * Temple_emline$flux

combined_spectrum = agn_spectrum
# Combine the AGN spectrum with the emission spectrum
combined_spectrum[41:3156] = agn_spectrum[41:3156] + emission_spectrum

return(data.frame(wave = Fritz_interp$Wave, lum = combined_spectrum))
}
