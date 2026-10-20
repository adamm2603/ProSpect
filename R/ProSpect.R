ProSpectSED = function(SFH = SFHfunc,
                          z = 0.1,
                          tau_birth = 1,
                          tau_screen = 0.3,
                          tau_AGN = 1,
                          pow_birth = -0.7,
                          pow_screen = -0.7,
                          pow_AGN = -0.7,
                          alpha_SF_birth = 1,
                          alpha_SF_screen = 3,
                          alpha_SF_AGN = 0,
                          AGNlum = 0,
                          sparse = 5,
                          speclib = NULL,
                          Dale = NULL,
                          AGN = NULL,
                          filtout = NULL,
                          filters = 'all',
                          Dale_M2L_func = NULL,
                          returnall = TRUE,
                          H0 = 67.8,
                          OmegaM = 0.308,
                          OmegaL = 1 - OmegaM,
                          waveout = seq(2, 9.35, by = 0.01),
                          ref,
                          unimax = 13.8e9,
                          agemax = NULL,
                          LumDist_Mpc = NULL,
                          addradio_SF = FALSE,
                          addradio_AGN = FALSE,
                          Te_SF = 1e4,
                          ff_frac_SF = 0.1,
                          ff_power_SF = -0.1,
                          sy_power_SF = -0.8,
                          Te_AGN = 1e4,
                          ff_frac_AGN = 0.1,
                          ff_power_AGN = -0.1,
                          sy_power_AGN = -0.8,
                          AGNct = 40,
                          AGNrm = 60,
                          AGNan = 30,
                          AGNta = 1,
                          AGNal = 4,
                          AGNbe = -0.5,
                          AGNp = 1,
                          AGNq = 1,
                          Eb = 0,
                          L0 = 2175.8,
                          LFWHM = 470,
                          IGMabsorb = 0,
                          Inoue14_LAFcoef = NULL,
                          Inoue14_DLAcoef = NULL,
                          AGN_gal_screen = TRUE,
                          AGN_emission = c("viewing_angle", "angle_averaged"),
                          ...) {
  AGN_emission = match.arg(AGN_emission)

  #call = match.call()

  if(!is.null(waveout)){
    waveout_max = max(waveout)
  }else{
    waveout_max = 9.35
  }

  if ('emission' %in% names(list(...))) {
    if (list(...)$emission & missing(waveout)) {
      waveout = NULL
    }
  }

  if(!is.null(SFH)){
    tau_birth = .interval(tau_birth, 0, 10, reflect = FALSE)
    tau_screen = .interval(tau_screen, 0, 10, reflect = FALSE)
    pow_birth = .interval(pow_birth, -2, 0, reflect = FALSE)
    pow_screen = .interval(pow_screen, -2, 0, reflect = FALSE)
    alpha_SF_birth = .interval(alpha_SF_birth, 0.0625, 4, reflect = FALSE)
    alpha_SF_screen = .interval(alpha_SF_screen, 0.0625, 4, reflect = FALSE)
  }

  if(!is.null(AGN)){
    tau_AGN = .interval(tau_AGN, 0, 10, reflect = FALSE)
    pow_AGN = .interval(pow_AGN, -2, 0, reflect = FALSE)
    alpha_SF_AGN = .interval(alpha_SF_AGN, 0.0625, 4, reflect = FALSE)
  }

  if(!is.null(SFH)){

    Stars = SFH(
      z = z,
      tau_birth = tau_birth,
      tau_screen = tau_screen,
      pow_birth = pow_birth,
      pow_screen = pow_screen,
      sparse = sparse,
      speclib = speclib,
      filters = NULL,
      unimax = unimax,
      agemax = agemax,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM,
      ...
    )

    if(!isFALSE(Dale)){
      Dust_Birth = Dale_interp(alpha_SF = alpha_SF_birth, Dale = Dale)
      Dust_Screen = Dale_interp(alpha_SF = alpha_SF_screen, Dale = Dale)

      SED_Bdust_Sdust = Dust_Birth$Aspec * Stars$lumtot_birth + Dust_Screen$Aspec *
        Stars$lumtot_screen
      SED_Stars_Bdust_Sdust = addspec(
        wave1 = Stars$wave_lum,
        flux1 = Stars$lum_atten,
        wave2 = Dust_Screen$Wave,
        flux2 = SED_Bdust_Sdust,
        extrap = 0,
        waveout = waveout
      )
    }else{
      Final = data.frame(wave = Stars$wave_lum, lum = Stars$lum_atten)
      Dust_Birth = NULL
      Dust_Screen = NULL
      SED_Bdust_Sdust = NULL
      SED_Stars_Bdust_Sdust = NULL
    }

    if (!is.null(Dale_M2L_func) & returnall) {
      dustlum_birth = Stars$lumtot_birth
      dustlum_screen = Stars$lumtot_screen
      dustmass_birth = Stars$lumtot_birth / Dale_M2L_func(alpha_SF_birth)
      dustmass_screen = Stars$lumtot_screen / Dale_M2L_func(alpha_SF_screen)
    } else{
      dustlum_birth = 0
      dustlum_screen = 0
      dustmass_birth = 0
      dustmass_screen = 0
    }

    ## adding radio contribution using just FIR generated from star-formation
    if (addradio_SF) {
      Final = radiocont(
        wave = SED_Stars_Bdust_Sdust$wave,
        flux = SED_Stars_Bdust_Sdust$flux,
        z = 0,
        Te = Te_SF,
        ff_frac = ff_frac_SF,
        ff_power = ff_power_SF,
        sy_power = sy_power_SF,
        wavesamp = seq(6, waveout_max, by=0.1),
        flux_in = 'wave',
        flux_out = 'wave'
      )
    }else if(!isFALSE(Dale)) {
      Final = SED_Stars_Bdust_Sdust
    }
  } else {
    Dust_Birth = NULL
    Dust_Screen = NULL
    SED_Bdust_Sdust = NULL
    SED_Stars_Bdust_Sdust = NULL
    dustlum_birth = NULL
    dustlum_screen = NULL
    dustmass_birth = NULL
    dustmass_screen = NULL
    Stars = NULL
    Final = NULL
  }

  if (is.null(AGN) | AGNlum == 0) {
    #Final = Final
    AGN = NULL
    dustlum_AGN = 0
    dustmass_AGN = 0
  } else {
    if(isFALSE(Dale)){
      stop('Dale cannot be FALSE when using an AGN model!')
    }

    if (AGN_emission == "angle_averaged" &&
      !inherits(AGN, "Fritz") && !inherits(AGN, "SKIRTOR")) {
      stop("Angle-averaged AGN emission requires a Fritz or SKIRTOR model")
        }

    if(inherits(AGN, 'Fritz') | inherits(AGN, 'SKIRTOR')) {
      # Save the original template before interpolation when angle-averaged
      # emission is requested.
      if (AGN_emission == "angle_averaged") {
        AGN_model = AGN
      }

      # Interpolate the selected line-of-sight AGN spectrum.
      if (inherits(AGN, 'Fritz')) {
        selected = Fritz_interp(
          lum = AGNlum,
          ct = AGNct,
          al = AGNal,
          be = AGNbe,
          ta = AGNta,
          rm = AGNrm,
          an = AGNan,
          Fritz = AGN
        )
      } else {
        selected = SKIRTOR_interp(
          lum = AGNlum,
          ta = AGNta,
          p = AGNp,
          q = AGNq,
          ct = AGNct,
          rm = AGNrm,
          an = AGNan,
          SKIRTOR = AGN
        )
      }

      dustlum_AGN = NA
      dustmass_AGN = NA
      host_atten = NULL

      if (AGN_emission == "angle_averaged") {
        # The line-of-sight spectrum is the direct observed component, while
        # the angle-averaged spectrum is used for the host-absorbed energy
        # budget and the associated dust re-emission.
        mean_SED = AGN_angle_averaged_SED(
          AGN_model, AGNlum, AGNct, AGNrm, AGNta,
          AGNal, AGNbe, AGNp, AGNq
        )

        if (isTRUE(AGN_gal_screen)) {
          direct = CF_screen_atten(
            wave = selected$wave,
            flux = selected$lum * .erg_to_lsol,
            tau = tau_screen,
            pow = pow_screen,
            Eb = Eb,
            L0 = L0,
            LFWHM = LFWHM
          )

          host_atten = atten_emit(
            wave = mean_SED$wave,
            flux = mean_SED$lum * .erg_to_lsol,
            tau = tau_screen,
            pow = pow_screen,
            alpha_SF = alpha_SF_screen,
            Dale = Dale,
            Dale_M2L_func = Dale_M2L_func,
            waveout = waveout,
            Eb = Eb,
            L0 = L0,
            LFWHM = LFWHM
          )

          AGN_final = addspec(
            wave1 = selected$wave,
            flux1 = direct$flux,
            wave2 = host_atten$emit$wave,
            flux2 = host_atten$emit$flux,
            extrap = 0,
            waveout = waveout
          )
        } else {
          # No host-galaxy screen: retain only the selected direct spectrum.
          AGN_final = data.frame(
            wave = selected$wave,
            flux = selected$lum * .erg_to_lsol
          )
        }
      } else {
        # Standard viewing-angle mode: apply the host screen to the selected
        # line-of-sight SED only when requested.
        if (isTRUE(AGN_gal_screen)) {
          host_atten = atten_emit(
            wave = selected$wave,
            flux = selected$lum * .erg_to_lsol,
            tau = tau_screen,
            pow = pow_screen,
            alpha_SF = alpha_SF_screen,
            Dale = Dale,
            Dale_M2L_func = Dale_M2L_func,
            waveout = waveout,
            Eb = Eb,
            L0 = L0,
            LFWHM = LFWHM
          )
          AGN_final = host_atten$final
        } else {
          AGN_final = data.frame(
            wave = selected$wave,
            flux = selected$lum * .erg_to_lsol
          )
        }
      }

      # Record dust absorption by the host screen only when that screen is on.
      if (isTRUE(AGN_gal_screen) &&
          !is.null(Dale_M2L_func) && returnall &&
          !is.null(host_atten)) {
        dustlum_screen = if (is.null(dustlum_screen)) 0 else dustlum_screen
        dustmass_screen = if (is.null(dustmass_screen)) 0 else dustmass_screen
        dustlum_screen = dustlum_screen + host_atten$total_atten
        dustmass_screen = dustmass_screen + host_atten$dustmass
      }

      # Subtract the AGN-associated radio continuum unless explicitly requested
      # to add it back, matching the existing ProSpectSED behaviour.
      AGN_final = radiocont(
        wave = AGN_final$wave,
        flux = AGN_final$flux,
        z = 0,
        Te = Te_AGN,
        ff_frac = ff_frac_AGN,
        ff_power = ff_power_AGN,
        sy_power = sy_power_AGN,
        wavesamp = seq(6, waveout_max, by = 0.1),
        flux_in = 'wave',
        flux_out = 'wave',
        subtractonly = !addradio_AGN
      )

      AGN = AGN_final
      if (is.null(Final)) {
        Final = AGN
      } else if (length(Final$flux) == length(AGN$flux)) {
        Final = data.frame(wave = Final$wave, flux = Final$flux + AGN$flux)
      } else {
        Final = addspec(
          wave1 = Final$wave,
          flux1 = Final$flux,
          wave2 = AGN$wave,
          flux2 = AGN$flux
        )
      }
      colnames(AGN)[2] = 'lum'
    } else{
      #Use old model
      #First we attenuate by the hot torus
      AGN = atten_emit(
        wave = AGN$Wave,
        flux = AGN$Aspec * AGNlum * .erg_to_lsol,
        tau = tau_AGN,
        pow = pow_AGN,
        alpha_SF = alpha_SF_AGN,
        Dale = Dale,
        Dale_M2L_func = Dale_M2L_func,
        waveout = waveout
      ) #no bump for the hot torus part

      if (!is.null(Dale_M2L_func) & returnall) {
        dustlum_AGN = AGN$total_atten
        dustmass_AGN = AGN$dustmass
      }
      # Apply the host-galaxy screen only when requested.
      if (isTRUE(AGN_gal_screen)) {
        AGN = atten_emit(
          wave = AGN$final$wave,
          flux = AGN$final$flux,
          tau = tau_screen,
          pow = pow_screen,
          alpha_SF = alpha_SF_screen,
          Dale = Dale,
          Dale_M2L_func = Dale_M2L_func,
          waveout = waveout,
          Eb = Eb,
          L0 = L0,
          LFWHM = LFWHM
        )

        if (!is.null(Dale_M2L_func) & returnall) {
          dustlum_screen = dustlum_screen + AGN$total_atten
          dustmass_screen = dustmass_screen + AGN$dustmass
        } else {
          dustlum_AGN = 0
          dustmass_AGN = 0
        }
      } else if (is.null(Dale_M2L_func) | !returnall) {
        dustlum_AGN = 0
        dustmass_AGN = 0
      }

      ## subtracts off AGN contribution to the radio continuum unless you specifically request to add it back
      AGN$final = radiocont(
        wave = AGN$final$wave,
        flux = AGN$final$flux,
        z = 0,
        Te = Te_AGN,
        ff_frac = ff_frac_AGN,
        ff_power = ff_power_AGN,
        sy_power = sy_power_AGN,
        wavesamp = seq(6, waveout_max, by = 0.1),
        flux_in = 'wave',
        flux_out = 'wave',
        subtractonly = !addradio_AGN # whether to add AGN radio or just subtract Dale radio
      )

      AGN = AGN$final
      if (is.null(Final)) {
        Final = AGN
      } else if (length(Final$flux) == length(AGN$flux)) {
        Final = data.frame(wave = Final$wave, flux = Final$flux +
                             AGN$flux)
      } else {
        Final = addspec(
          wave1 = Final$wave,
          flux1 = Final$flux,
          wave2 = AGN$wave,
          flux2 = AGN$flux
        )
      }
      colnames(AGN)[2] = 'lum'
    }
  }
  colnames(Final)[2] = 'lum'

  ## Previous ProSpect had IGM absorb done on the rest frame
  ## Now do to the observer frame later on
  # if (IGMabsorb > 0) {
  #   sel = which(Final$wave < 1215.67)
  #   Final$lum[sel] = Final$lum[sel] * (1 - IGMabsorb)
  #   sel = which(Final$wave < 911.75)
  #   Final$lum[sel] = 0
  # }

  if (is.null(filtout) & !is.null(filters)) {
    if (filters[1] == 'all') {
      cenwave = NULL
      data('cenwave', envir = environment())
      filters = cenwave$filter
    }else if(filters[1] == 'GAMA'){
      filters = c(
        'FUV_GALEX',
        'NUV_GALEX',
        'u_VST',
        'g_VST',
        'r_VST',
        'i_VST',
        'Z_VISTA',
        'Y_VISTA',
        'J_VISTA',
        'H_VISTA',
        'K_VISTA',
        'W1_WISE' ,
        'W2_WISE',
        'W3_WISE',
        'W4_WISE',
        'P100_Herschel',
        'P160_Herschel',
        'S250_Herschel' ,
        'S350_Herschel',
        'S500_Herschel'
      )
    }else if(filters[1] == 'WAVES'){
      filters = c(
        'u_VST',
        'g_VST',
        'r_VST',
        'i_VST',
        'Z_VISTA',
        'Y_VISTA',
        'J_VISTA',
        'H_VISTA',
        'K_VISTA',
        'W1_WISE' ,
        'W2_WISE'
      )
    }

    filtout = list()
    for (i in filters) {
      filtout = c(filtout, list(approxfun(getfilt(i))))
    }
    names(filtout) = filters
  }

  if (z > 0 & !is.null(filtout)) {
    Flux = Lum2Flux(
      wave = Final$wave,
      lum = Final$lum,
      z = z,
      H0 = H0,
      OmegaM = OmegaM,
      OmegaL = OmegaL,
      ref = ref,
      LumDist_Mpc = LumDist_Mpc
    )
    Flux$flux = convert_wave2freq(flux_wave = Flux$flux * .cgs_to_jansky,
                                  wave = Flux$wave)

    if (IGMabsorb > 0 & is.numeric(IGMabsorb)){
      sel = which(Flux$wave/(1+z) < 1215.67)
      Flux$flux[sel] = Flux$flux[sel] * (1 - IGMabsorb)
      sel = which(Flux$wave/(1+z) < 911.75)
      Flux$flux[sel] = 0
    }else if (IGMabsorb == "Inoue14"){
      Flux$flux = Flux$flux * Inoue14_IGM(Flux$wave, z, Inoue14_LAFcoef, Inoue14_DLAcoef)
    }

    photom_out = {}
    for (i in 1:length(filtout)) {
      photom_out = c(photom_out,
                     bandpass(
                       flux = Flux$flux,
                       wave = Flux$wave,
                       filter = filtout[[i]]
                     ))
    }
  } else if (z > 0 & is.null(filtout)) {
    Flux = Lum2Flux(
      wave = Final$wave,
      lum = Final$lum,
      z = z,
      H0 = H0,
      OmegaM = OmegaM,
      OmegaL = OmegaL,
      ref = ref,
      LumDist_Mpc = LumDist_Mpc
    )
    Flux$flux = convert_wave2freq(flux_wave = Flux$flux * .cgs_to_jansky,
                                  wave = Flux$wave)
    if (IGMabsorb > 0 & is.numeric(IGMabsorb)){
      sel = which(Flux$wave/(1+z) < 1215.67)
      Flux$flux[sel] = Flux$flux[sel] * (1 - IGMabsorb)
      sel = which(Flux$wave/(1+z) < 911.75)
      Flux$flux[sel] = 0
    }else if (IGMabsorb == "Inoue14"){
      Flux$flux = Flux$flux * Inoue14_IGM(Flux$wave, z, Inoue14_LAFcoef, Inoue14_DLAcoef)
    }

    photom_out = Flux
  } else if (z <= 0 & !is.null(filtout)) {
    Flux = cbind(wave = Final$wave,
                 flux = Final$lum * .lsol_to_absolute)
    photom_out = photom_flux(Flux, outtype = 'magAB', filters = filtout)
    #photom_out = photom_lum(Final, z=0, outtype = 'magAB', filters = filtout) same as above, but redundant
  } else{
    Flux = NULL
    photom_out = NULL
  }

  if (returnall) {
    StarsAtten = data.frame(wave = Stars$wave_lum, lum = Stars$lum_atten)
    StarsUnAtten = data.frame(wave = Stars$wave, lum = Stars$lum_unatten)

    if(!isFALSE(Dale)){
      DustEmit = data.frame(wave = Dust_Screen$Wave, lum = SED_Bdust_Sdust)
    }else{
      DustEmit = NULL
    }

    output = list(
      Photom = photom_out,
      FinalFlux = Flux,
      FinalLum = Final,
      StarsAtten = StarsAtten,
      StarsUnAtten = StarsUnAtten,
      DustEmit = DustEmit,
      AGN = AGN,
      Stars = Stars,
      dustmass = c(
        birth = dustmass_birth,
        screen = dustmass_screen,
        AGN = dustmass_AGN,
        total = sum(c(
          dustmass_birth, dustmass_screen, dustmass_AGN
        ), na.rm = TRUE)
      ),
      dustlum = c(
        birth = dustlum_birth,
        screen = dustlum_screen,
        AGN = dustlum_AGN,
        total = sum(c(
          dustlum_birth, dustlum_screen, dustlum_AGN
        ), na.rm = TRUE)
      ),
      #call = call,
      z = z,
      filters = filters,
      filtout = filtout,
      cosmo = list(H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    )
    class(output) = 'ProSpectSED'
    return(output)
  } else{
    return(photom_out)
  }
}

