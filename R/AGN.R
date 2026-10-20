Fritz_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, al = 4, be = -0.5,
                        Fritz = NULL){
  if(is.null(Fritz)){
    Fritz=NULL
    data('Fritz', envir = environment())
  }

  ctmix = interp_quick(ct, Fritz$ct)
  almix = interp_quick(al, Fritz$al)
  bemix = interp_quick(be, Fritz$be)
  tamix = interp_quick(ta, Fritz$ta)
  rmmix = interp_quick(rm, Fritz$rm)
  anmix = interp_quick(an, Fritz$an)

  slice = Fritz$Aspec[
    c(ctmix[1:2]),
    c(almix[1:2]),
    c(bemix[1:2]),
    c(tamix[1:2]),
    c(rmmix[1:2]),
    c(anmix[1:2]),
  ]

  weights = rep(1,64)
  weights = weights * rep(ctmix[3:4], each=1, times=32)
  weights = weights * rep(almix[3:4], each=2, times=16)
  weights = weights * rep(bemix[3:4], each=4, times=8)
  weights = weights * rep(tamix[3:4], each=8, times=4)
  weights = weights * rep(rmmix[3:4], each=16, times=2)
  weights = weights * rep(anmix[3:4], each=32, times=1)
  tempmat = matrix(as.numeric(slice),64,178)
  return(data.frame(wave=Fritz$Wave, lum=lum * colSums(tempmat * weights)))
}

SKIRTOR_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
                          SKIRTOR = NULL){

  if(is.null(SKIRTOR)){
    data('SKIRTOR', envir = environment())
  }


  tmix = interp_quick(ta, SKIRTOR$ta)
  pmix = interp_quick(p, SKIRTOR$p)
  qmix = interp_quick(q, SKIRTOR$q)
  oamix = interp_quick(ct, SKIRTOR$ct)
  rmmix = interp_quick(rm, SKIRTOR$rm)
  imix = interp_quick(an, SKIRTOR$an)

  slice = SKIRTOR$Aspec[c(tmix[1:2]), c(pmix[1:2]), c(qmix[1:2]),
                        c(oamix[1:2]), c(rmmix[1:2]), 1, c(imix[1:2]), ]

  weights = rep(1, 64)
  weights = weights * rep(tmix[3:4], each = 1, times = 32)
  weights = weights * rep(pmix[3:4], each = 2, times = 16)
  weights = weights * rep(qmix[3:4], each = 4, times = 8)
  weights = weights * rep(oamix[3:4], each = 8, times = 4)
  weights = weights * rep(rmmix[3:4], each = 16, times = 2)
  weights = weights * rep(imix[3:4], each = 32, times = 1)

  tempmat = matrix(as.numeric(slice), 64, 132)
  agn_spectrum = (colSums(tempmat * weights))/SKIRTOR$Wave


  agn_spectrum = (agn_spectrum/(3.839e33 * 1E11))*lum
  return(data.frame(wave=SKIRTOR$Wave, lum = agn_spectrum))

}

AGN_angle_averaged_SED <- function(model, AGNlum, AGNct, AGNrm, AGNta,
                                   AGNal, AGNbe, AGNp, AGNq) {
    angles <- model$an
    if (length(angles) < 2L || anyNA(angles) ||
        is.unsorted(angles, strictly = TRUE) ||
        angles[1L] != 0 || tail(angles, 1L) != 90) {
        stop("AGN viewing angles must increase from 0 to 90 degrees")
    }
    mu <- cos(angles * pi / 180)
    weights <- c((mu[1L] - mu[2L])/2,
                 (mu[seq_len(length(mu) - 2L)] - mu[3:length(mu)])/2,
                 (mu[length(mu) - 1L] - mu[length(mu)])/2)
    weights <- weights/sum(weights)


    spectrum_at <- function(angle) {
        if (inherits(model, "Fritz")) {
            Fritz_interp(lum = AGNlum, ct = AGNct, al = AGNal,
                         be = AGNbe, ta = AGNta, rm = AGNrm,
                         an = angle, Fritz = model)
        } else {
            SKIRTOR_interp(lum = AGNlum, ct = AGNct, ta = AGNta,
                           p = AGNp, q = AGNq, rm = AGNrm,
                           an = angle, SKIRTOR = model)
        }
    }
    first <- spectrum_at(angles[1L])
    mean_lum <- weights[1L] * first$lum
    for (i in seq.int(2L, length(angles))) {
        spectrum <- spectrum_at(angles[i])
        if (!identical(spectrum$wave, first$wave)) {
            stop("AGN viewing angles must share a wavelength grid")
        }
        mean_lum <- mean_lum + weights[i] * spectrum$lum
    }
    data.frame(wave = first$wave, lum = mean_lum)
} 
