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


ProSpect_component_on_grid <- function(component, wave) {
    if (is.null(component)) {
        return(rep(0, length(wave)))
    }
    if (identical(component$wave, wave)) {
        return(component$lum)
    }
    10^approx(log10(component$wave), log10(component$lum),
              xout = log10(wave), yleft = -Inf, yright = -Inf)$y
}


AGN_angle_averaged_SED = function (model, AGNlum, AGNct, AGNrm, AGNta,
                                   AGNal, AGNbe, AGNp, AGNq) {
    angles = model$an
    if (length(angles) != 10 || anyNA(angles) || any(angles !=
        seq(0, 90, by = 10))) {
        stop("AGN viewing angles must be 0, 10, ..., 90 degrees")
    }
    if (!inherits(model, "Fritz") && !inherits(model, "SKIRTOR")) {
        stop("AGN model must be Fritz or SKIRTOR")
    }
    weights = c(0.0868240888334652, 0.171010071662834, 0.163175911166535,
        0.150383733180435, 0.133022221559489, 0.11161889704895,
        0.0868240888334652, 0.0593911746138847, 0.0301536896070458,
        0.00759612349389599)
    for (i in 1:10) {
        if (inherits(model, "Fritz")) {
            spectrum = Fritz_interp(lum = AGNlum, ct = AGNct,
                al = AGNal, be = AGNbe, ta = AGNta, rm = AGNrm,
                an = angles[i], Fritz = model)
        }
        else {
            spectrum = SKIRTOR_interp(lum = AGNlum, ct = AGNct,
                ta = AGNta, p = AGNp, q = AGNq, rm = AGNrm, an = angles[i],
                SKIRTOR = model)
        }
        if (i == 1) {
            wave = spectrum$wave
            mean_lum = weights[i] * spectrum$lum
        }
        else {
            if (!identical(spectrum$wave, wave)) {
                stop("AGN viewing angles must share a wavelength grid")
            }
            mean_lum = mean_lum + weights[i] * spectrum$lum
        }
    }
    return(data.frame(wave = wave, lum = mean_lum))
}
