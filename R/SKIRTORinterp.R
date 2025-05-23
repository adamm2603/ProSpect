SKIRTOR_interp = function(lum = 1e+44, t = 7, p = 1, q = 1, oa = 40, R = 20,
                             i=60, em_factor = 1, SKIRTOR = NULL, Temple_emline = NULL){
  
  if(is.null(SKIRTOR)){
    data('SKIRTOR', envir = environment()) 
  }
  
  if(is.null(Temple_emline)){
    data('Temple_emline', envir = environment())
  }
  
  tmix = interp_quick(t, SKIRTOR$t)
  pmix = interp_quick(p, SKIRTOR$p)
  qmix = interp_quick(q, SKIRTOR$q)
  oamix = interp_quick(oa, SKIRTOR$oa)
  Rmix = interp_quick(R, SKIRTOR$R)
  imix = interp_quick(i, SKIRTOR$i)
 
  slice = SKIRTOR$Aspec[c(tmix[1:2]), c(pmix[1:2]), c(qmix[1:2]),
                      c(oamix[1:2]), c(Rmix[1:2]), 1, c(imix[1:2]), ]
   
  weights = rep(1, 64)
  weights = weights * rep(tmix[3:4], each = 1, times = 32)
  weights = weights * rep(pmix[3:4], each = 2, times = 16)
  weights = weights * rep(qmix[3:4], each = 4, times = 8)
  weights = weights * rep(oamix[3:4], each = 8, times = 4)
  weights = weights * rep(Rmix[3:4], each = 16, times = 2)
  weights = weights * rep(imix[3:4], each = 32, times = 1)
  
  tempmat = matrix(as.numeric(slice), 64, 3224)
  agn_spectrum = (colSums(tempmat * weights))/SKIRTOR$Wave
  
  Temple_emline$flux = Temple_emline$flux * agn_spectrum[433]
  # Apply the emission factor to the emission spectrum
  emission_spectrum = em_factor * Temple_emline$flux
  combined_spectrum = agn_spectrum 
  combined_spectrum[35:3150] = agn_spectrum[35:3150] + emission_spectrum
  
  
  agn_spectrum = (combined_spectrum/(3.839e33 * 1E11))*lum
  return(data.frame(wave=SKIRTOR$Wave, lum = agn_spectrum)) 

}