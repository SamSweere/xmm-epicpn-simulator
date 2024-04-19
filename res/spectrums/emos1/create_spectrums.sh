#!/bin/bash

wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/m1tffg_events.fits
wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/m1mffg_events.fits
wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/m1kffg_events.fits

for MODE in m1tffg m1mffg m1kffg
do
  evselect table="$MODE"_events.fits energycolumn=PI spectrumset="$MODE"_spectrum.fits specchannelmin=0 specchannelmax=20479 withspecranges=yes withspectrumset=yes
  rm "$MODE"_events.fits
done