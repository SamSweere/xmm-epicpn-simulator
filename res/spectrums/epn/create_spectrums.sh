#!/bin/bash

wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/pntffg_events.fits
wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/pnmffg_events.fits
wget https://xmm-tools.cosmos.esa.int/external/xmm_calibration/background/bs_repository/pnkffg_events.fits

for MODE in pntffg pnmffg pnkffg
do
  evselect table="$MODE"_events.fits energycolumn=PI spectrumset="$MODE"_spectrum.fits specchannelmin=0 specchannelmax=20479 withspecranges=yes withspectrumset=yes
  rm "$MODE"_events.fits
done