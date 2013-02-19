pro plot_rms_spec
;; Plots the RMS along the time series for each wavelength in the
;; spectrum

  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  restore,'data/timedata.sav'
