#!/bin/csh
# Numeric fix
setenv LC_NUMERIC "en_US.UTF-8"

set nArgs=2
set ScriptVersion=1.3
set ScriptName=cmp_psf_SIXTE.csh
if ($#argv <  $nArgs) then 
	echo
	echo "  Script  :  ${ScriptName} --  version ${ScriptVersion}"
	echo 
	echo "  Purpose :  Generates SIXTE-compatible calibration files for the XMM psf"
	echo
	echo "  Usage   :  ${ScriptName}  CIFfil OutFile [instrument=M1/M2/PN] [model=...] [psfrad=...] [stretch=...]"
	echo "                                [energy='... ... ...'] [offaxis='... ... ...'] [azimuth='... ... ...']"
	echo
	exit
endif

# Main parameters
set ccf="$1"
set outPSF="$2"
if !(-e ${ccf}) then
	echo "Input calibration index '${ccf}' does not exist"
	exit
endif

# Defaults (sampling from SIXTE instruemnt file v1.1.0)
set binning=-1
set ImWidth=2.5
set lst=None
set model=ELLBETA
set InstName=Unknown
set StretchFactor=1.0
set energy_array=( 300 600 1000 3000 6000 9000 12000 )
set theta_array=( 0 210 420 600 720 900 1200 )
set phi_array=( 0 )
@ i=0 
while ($i < 89 )
	@ i++
	@ newVal=$phi_array[$i]
	@ newVal+=4
	set phi_array=( ${phi_array} $newVal )
end

# Read optional arguments
if ($#argv > $nArgs) then
  while ($nArgs < $#argv )
      @ nArgs++
      set optstr="${argv[$nArgs]}" 
      # Check option format
      set check=`echo $optstr | grep =`
      set opttype=`echo $optstr | cut -d= -f1`
      set optval=`echo $optstr | cut -d= -f2`
      if ( ("$check" != "$optstr") || \
           ("a${opttype}a" == 'aa') || ("a${optval}a" == 'aa') ) then
     	  	echo "Invalid format for pipeline option: '$optstr'."
		exit
      endif  
      # Set option
      set check=0
      set opttype=`echo $optstr | cut -d= -f1 | tr '[a-z]' '[A-Z]'`
      if ($opttype == 'INSTRUMENT') then
      		set InstOK=0
		set optval=`echo $optval | tr '[a-z]' '[A-Z]'`
		if (($optval == M1) || ($optval == MOS1) || ($optval == EMOS1)) then 
			set InstOK=1
			set XRT=1
			set InstName=MOS1
      		endif
		if (($optval == M2) || ($optval == MOS2) || ($optval == EMOS2)) then 
			set InstOK=1
			set XRT=2
			set InstName=MOS2
      		endif
		if (($optval == PN) || ($optval == EPN)) then 
			set InstOK=1
			set InstName=PN
			set XRT=3
      		endif
		if ($InstOK != 1) then
			echo "Unkown instrument '${optval}'"
			exit
		endif
		set check=1
      endif  	
      #if ($opttype == 'EVENTLIST') then
	#	set lst=$optval
	#	if !(-e ${lst}) then
	#		echo "Input event list '${lst}' does not exist"
	#		exit
	#	endif	
      	#	set check=1
      #endif  	
      #if ($opttype == 'PIXSIZE') then
      #		set binning=$optval
      #		set check=1
      #endif  	
      if ($opttype == 'PSFRAD') then
		set ImWidth=$optval
		set check=1
      endif  	
      if ($opttype == 'STRETCH') then
		set StretchFactor=$optval
		set check=1
      endif  	
      if ($opttype == 'ENERGY') then
		set energy_array=( $optval )
		set check=1
      endif  	
      if ($opttype == 'OFFAXIS') then
		set theta_array=( $optval )
		set check=1
      endif  	
      if ($opttype == 'AZIMUTH') then
		set phi_array=( $optval )
		set check=1
      endif  	
      if ($opttype == 'MODEL') then
		set model=`echo $optval | tr '[a-z]' '[A-Z]' `
		set check=1
      endif  	
      # Final check
      if ($check == 0) then
      	echo "Unknown optional argument: '$opttype'"
      	exit
      endif
  end
endif

# Check PSF model
if ( ($model != LOW) && ($model != MEDIUM) && ($model != EXTENDED) && ($model != HIGH) && ($model != ELLBETA)) then
	echo "Unknown psf model '${model}'"
	exit
endif
if (($model == MEDIUM) || ($model == HIGH)) then
	if ($binning != -1) then
		echo 'PSFGEN imposes pixels psf 1.1" for psf images in '${model}' - ignoring the BINNING option"
	endif
	set binning=1.1
endif
if (($model == EXTENDED) || ($model == ELLBETA)) then
	if ($binning != -1) then
		echo 'PSFGEN imposes pixels psf 1.0" for psf images in '${model}' - ignoring the BINNING option"
	endif
	set binning=1.0
endif
#if (($model == HIGH) && ($binning == -1)) then
#	echo 'Pixel size was not specified - using the default of 1.1" for model '"'HIGH'"
#	set binning=1.1
#endif

# Check CCF file
if (`ftlist ${ccf} H | grep CALINDEX | wc -l` != 1) then
	echo "Could not find Calibration index extension in CIF file"
	exit
endif
setenv SAS_CCF $ccf
# Get MISCDATA file
set MISCDnum=`fdump infile=${ccf}'[CALINDEX]' outfile=STDOUT columns=SCOPE,TYPEID,ISSUE rows=- prhead=no more=yes | grep XMM | grep MISCDATA | awk '{printf("%04i\n",$4)}'`
if ("a${MISCDnum}a" == 'aa') then
	echo Could not find XMM MISCDATA entry in CIF
	exit  
endif
set expandCCFPATH=( `echo ${SAS_CCFPATH} | tr ':' ' '` )
set MISCDATAfile=`find ${expandCCFPATH}  -name XMM_MISCDATA_${MISCDnum}.CCF`
if ("a${MISCDATAfile}a" == 'aa') then
	echo Could not find XMM MISCDATA file in CCF path
	exit  
endif
if (`ftlist ${MISCDATAfile} H | grep MISCDATA | wc -l` != 1) then
	echo "Format of CCF XMM MISCDATA file is not recognized"
	exit
endif
# Get observation date
#fkeypar fitsfile=${SAS_CCF}'[CALINDEX]' keyword=OBSVDATE
#if (`pget fkeypar exist` != 'yes') then 
#	echo "Could not get observation date from input CIF"
#	exit
#endif
#set ObsDate=`pget fkeypar value | cut -d"'" -f2 | cut -dT -f1`
#set ObsTime=`pget fkeypar value | cut -d"'" -f2 | cut -dT -f2`

# Check instrument and set focal length
if ($InstName == Unknown) then
	if ($lst == None) then
		echo 'Instrument must be specified if no event-list is given as input'
		exit
	else 
		fkeypar fitsfile=${lst}'[EVENTS]' keyword=INSTRUME
		if (`pget fkeypar exist` != 'yes') then 
			echo "Could not find instrument name in file '${lst}'"
			exit
		endif
		switch (`pget fkeypar value`)
			case "'EMOS1'":  
				set XRT=1
				set InstName=MOS1
				breaksw
			case "'EMOS2'":  
				set XRT=2
				set InstName=MOS2
				breaksw
			case "'EPN'":  
				set XRT=3
				set InstName=PN
				breaksw
			default:  
				set outname=`pget fkeypar value`
				echo "Instrument `pget fkeypar value` is not a valid EPIC name"
				exit
				breaksw
		
		endsw
		
	endif
else 

	# Check that list and instrument are compatible
	if ($lst != None) then
		fkeypar fitsfile=${lst}'[EVENTS]' keyword=INSTRUME
		if (`pget fkeypar exist` != 'yes') then 
			echo "Could not find instrument name in file '${lst}'"
			exit
		endif
		if ( `pget fkeypar value | cut -d"'" -f2` != "E${InstName}" ) then
			echo Input list is not compatible with selected instrument
			exit
		endif
	endif

endif
set FOCAL_LENGTH=`fdump infile=${MISCDATAfile}'[MISCDATA]' outfile=STDOUT prhead=no more=yes columns=INSTRUMENT_ID,PARM_ID,PARM_VAL rows=- fldsep='|' pagewidth=200 | grep XRT${XRT} | grep FOCAL_LENGTH | awk -F'|' '{printf("%.5f\n",$4/1000.0)}'`

set SAS_binning=`echo $binning | awk '{printf("%s",$1/0.05)}'`
if ( ${SAS_binning}  != `echo $binning | awk '{print $1/0.05}'`) then 
	echo 'Please choose a binning that is a multiple of 0.05"'
	exit
endif
@ half=`echo $binning $ImWidth | awk 'function ceil(x,y){y=int(x); y=(x>y?y+1:y); return(y)} {print ceil($2*60/$1)}'`
@ CRPIX=${half} + 1
@ width=${half}
@ width*=2
set CDELT=`echo $binning $FOCAL_LENGTH $StretchFactor | awk '{print $1*$3/3600.0/180.0*3.14159265*$2}'`
echo
echo "Computing XMM PSF ($model) for instrument EPIC-${InstName} in ${binning}"'" pixels:'
echo
echo "  For a radius of ${ImWidth}', images will be ${width} pixels wide, with a scale of ${CDELT} m/pixel"
echo

# Create temporry folder
set datestr=`date '+%s%N'`
set tmpDir=/tmp/XMMpsfScript${datestr}
if !(-d ${tmpDir}) mkdir ${tmpDir}

# Get reference image for psfgen 
if ( ( ${model} == HIGH ) && ( $binning != 1.1 ) ) then
	echo "  Extracting reference image ..."
	set RefImg=${tmpDir}/Refimg${datestr}.fits
	evselect table=$lst withimageset=yes withimagedatatype=yes \
   		expression="((PI in [500:2000]) && (PATTERN in [0:4]) && (FLAG & 0x2fb002c)==0 )" \
   		filtertype=expression imageset=${RefImg} imagebinning=binSize \
   		xcolumn='X' ximagebinsize=${SAS_binning} ycolumn='Y' yimagebinsize=${SAS_binning} \
   		imagedatatype='Real64' updateexposure=yes writedss=yes filteredset=${tmpDir}/tmp_list.fits -V 0 -w 0
		fparkey value=0.0 fitsfile=${RefImg}'[0]' keyword=PA_PNT add=no
		/bin/rm -f ${tmpDir}/tmp_list.fits
	echo
	set ImString=" image=${RefImg}"
else
	set ImString=''
endif

# Extract PSF
@ i=0
@ nHDU=-1
while (${i} < $#energy_array)
	
	# Set energy
	@ i++
	@ j=0
	set energy=${energy_array[$i]}
	
	while (${j} < $#theta_array )
		
		# Set  offaxis
		@ j++
		@ k=0
		set theta=${theta_array[$j]}
	
		while (${k} < $#phi_array )

			# Set  offaxis
			@ k++
			set phi=${phi_array[$k]}
			@ nHDU +=1
			
			# Parameters
			echo "  Obtaining PSF image for E=${energy}eV, Theta=${theta}"'"'", Phi=${phi}deg ..." 
			set psf_img=${tmpDir}/PsfImg_${energy}eVthet${theta}arcminphi${phi}deg.fits
	
			# Get image from PSF Gen
			set phi_SAS=`echo $phi | awk '{ phi = 90.0 - $1 ; if(phi<0){phi+=360}; print phi/180.0*3.14159265}'`
			
			psfgen instrument=E${InstName} coordtype=TEL x=${theta} y=${phi_SAS} output=${psf_img} level=${model} xsize=${width} ysize=${width} energy=${energy} ${ImString} -V 0 -w 0
			# Normalize and convert to double
			ftstat ${psf_img}+0 > /dev/null
			fcarith infile=${psf_img} const=`pget ftstat sum`  outfil=${psf_img} ops=DIV copyprime=no clobber=yes datatype=d
			# Add OGIP keywords
			#fparkey fitsfile=${psf_img}+0 keyword=HISTORY  value="created by cmp_SIXTE_psf.csh based on the ELLBETA model in psfgen" add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CDES0001 value='ELLBETA model PSF' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CVST0001 value='00:00:00' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CVSD0001 value='2000-01-01' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CBD30001 value="PHI( ${phi}.000000)deg" add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CBD20001 value="THETA( ${theta}.000000)arcsec" add=yes insert=COMMENT > /dev/null
			set EnergyFormatted=`echo $energy | awk '{printf("%.3f",$1)}'`
			fparkey fitsfile=${psf_img}+0 keyword=CBD10001 value="ENERGY( ${EnergyFormatted})keV" add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CCNM0001 value='2D_PSF  ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CDTP0001 value='DATA    ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CCLS0001 value='BCF     ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=SUMRCTS  value=1. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CHANTYPE value='PI	 ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CHANMAX  value=-99. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CHANMIN  value=-99. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=ENERG_HI value=1. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=ENERG_LO value=1. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CDELT2   value=${CDELT}  add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CDELT1   value=${CDELT}  add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CRVAL2   value=0. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CRVAL1   value=0. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CRPIX2   value=${CRPIX}. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CRPIX1   value=${CRPIX}. add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CUNIT2   value='m       ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CUNIT1   value='m       ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CTYPE2   value='DETY    ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=CTYPE1   value='DETX    ' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=PHI      value=${phi}. add=yes insert=COMMENT comm='Azimuth in degree' > /dev/null
			set Theta_amin=`echo $theta | awk '{printf("%.3f", $1/60.0)}'`
			fparkey fitsfile=${psf_img}+0 keyword=THETA    value=${Theta_amin} add=yes insert=COMMENT comm='Off-axis angle in arcmin' > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=ENERGY   value=${energy}. add=yes insert=COMMENT comm='Energy in eV' > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=BACKGRND value=0.0 add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=FILTER   value='NONE	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=INSTRUME value="${InstName}" add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=TELESCOP value='XMM-Newton' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUCLAS4 value='NET	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUCLAS3 value='PREDICTED' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUCLAS2 value='PSF	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUCLAS1 value='IMAGE	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUVERS  value='1.0.0	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUDOC   value='CAL/GEN/92-027' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=HDUCLASS value='OGIP	' add=yes insert=COMMENT > /dev/null
			fparkey fitsfile=${psf_img}+0 keyword=EXTNAME  value="${energy}eVthet${theta}arcsecphi${phi}deg" add=yes insert=COMMENT > /dev/null
			#fparkey fitsfile=${psf_img}+0 keyword=GCOUNT   value=1 comm="required keyword; must = 1" add=yes insert=COMMENT > /dev/null
			#fparkey fitsfile=${psf_img}+0 keyword=PCOUNT   value=0 comm="required keyword; must = 0" add=yes insert=COMMENT > /dev/null

			fthedit infile=${psf_img}+0 keyword=COMMENT operation=deleteall
			
			if ($nHDU == 0) then 
				/bin/mv -f $psf_img $outPSF
			else
				fappend infile=${psf_img}+0 outfile=$outPSF pkeywds=yes history=no
				/bin/rm -f $psf_img
			endif
			
		end
	end
end

/bin/rm -rf ${tmpDir}

echo
echo Done.
echo
