#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "Sarfia"


Function newWiener(waveName)

String waveName
string/g globalwave = waveName
	Variable /G threshold
	Variable /G tau//=0.06
	Variable /G taurise//=0.001
	Variable/C noiseVariance
	Variable sr = 1/(dimdelta($waveName, 0))
	Variable deltat=dimdelta($waveName, 0)
	Variable Nsmoothpoints = 21
	
	Variable wind = .2
	//////////////////////////////////////////////////////
	Variable/g unitaryEvent = 0.145	//Estimate of the amplitude of an individual event (usually 0.2 for bipolar cells)
	/////////////////////////////////////////////////////
	
	Variable scalingfactor=unitaryEvent/tau  //This scaling factor gets the amplitudes right
	
	duplicate/o $waveName, xsignal, h, xfiltered, wave_smth
	Smooth/S=2 Nsmoothpoints, wave_smth				//Savitsky-Golay

	//STEP 1: Make the impulse response, h 	
	h=unitaryEvent*(1-exp(-x/taurise))*exp(-x/tau)
	//h=unitaryEvent*exp(-x/tau)
	Variable scaleunitary=(unitaryevent/wavemax(h))
	Duplicate/O h kernel_smth	
	
	string firstFilterName = waveName+ "_k1"
	duplicate/o h, $firstFilterName
	
	//STEP 2: Make the FT of the impulse response (H) and the power spectrum (HMAGSQR)
	FFT/OUT=1/DEST=HF  h
	FFT/OUT=4/DEST=HMAGSQR  h	
	
	//STEP 3: Make FT of the original signal to get mean power density of original signal, S
	FFT/OUT=1/DEST=SF  wave_smth
	
	//STEP 4: Make FT of the noise to get mean power density, N
	String RawHistName=waveName+"_H"					//Histogram of values in the DF/F wave
	String fitRawHistName="fit_"+RawHistName		//For fitting a Gaussian to the first part _H wave so as to estimate baseline SD 
	Make/N=0/O $RawHistName
	Histogram/B=5 $waveName,$RawHistName;DelayUpdate 
	FindPeak/Q/P $RawHistName
	Make/O/D/N=0 levelsWave
	Variable level = 0.8*waveMax($RawHistName)		//Set points for fit at either side of 0.5 the peak (avoids tail in distribution, reflecting responses)
	FindLevels/Q/N=2/D=levelsWave $RawHistName, level
	K0 = 0													//Constrain fit to baseline = 0
	Make/O/N=4 W_coef
	CurveFit/Q/X=1/H="1000"/TBOX=768 gauss $RawHistName[0,x2pnt($RawHistName,levelsWave[1])]
	noiseVariance = (W_coef[3]/2.355)^2   //FWHM = 2.355 * sd
	print W_coef[3]
	
	//STEP 5: Calculate SNR(f)
	Duplicate/O/C SF SNR
	Duplicate/O/C SF N
	N=noiseVariance
	SNR = (unitaryEvent*tau*HMAGSQR)/N   //I don't know quite why I need to scale in this way! SNR is frequency-dependent.  Gaussian power spectrum is flat

	//STEP 6: Calculate the FT of the Wiener deconvolution filter, G
	Duplicate/O/C SF G
	G = (HF^-1)*( HMAGSQR / (HMAGSQR + SNR^-1) )
		
	//STEP7: Convolution filtering of the original, xsignal, and the Wiener filter
	FFT/OUT=1/DEST=XSIG  xsignal
	Duplicate/O/C XSIG XCONVG
	XCONVG=XSIG*G
	
	IFFT/DEST=xfiltered  XCONVG
	xfiltered*=scalingFactor
	
	String smoothwaveName = waveName +"_S"
	String DeconWaveName = waveName +"_D"
	Duplicate/O xfiltered $DeconWaveName
	Duplicate/O wave_smth $smoothWaveName
	
	killwaves SNR, SF, N, levelsWave, kernel_smth, HMAGSQR, HF, h, G, xCONVG,xfiltered,xSig,xsignal
	killwaves w_coef, wave_smth
end





////////////////////// aveTempWeiner(waveName,filterName) /////////////////////////

Function aveTempWiener(waveName,filterName)

String waveName
String filterName
	
	Variable /G threshold
	Variable /G tau//=0.06
	Variable /G taurise//=0.001
	Variable/C noiseVariance
	Variable sr = 1/(dimdelta($waveName, 0))
	Variable deltat=dimdelta($waveName, 0)
	Variable Nsmoothpoints = 21
	
	Variable wind = .2
	//////////////////////////////////////////////////////
	Variable unitaryEvent = 0.145	//Estimate of the amplitude of an individual event (usually 0.2 for bipolar cells)
	/////////////////////////////////////////////////////
	
	Variable scalingfactor=unitaryEvent/tau  //This scaling factor gets the amplitudes right
	
	duplicate/o $waveName, xsignal, h, xfiltered, wave_smth
	Smooth/S=2 Nsmoothpoints, wave_smth				//Savitsky-Golay

	//STEP 1: Make the impulse response, h 	
	duplicate/o $filterName, h
	//h=unitaryEvent*exp(-x/tau)
	Variable scaleunitary=(unitaryevent/wavemax(h))
	Duplicate/O h kernel_smth	
	
	string firstFilterName = waveName+ "_k1"
	//duplicate/o h, $firstFilterName
	
	//STEP 2: Make the FT of the impulse response (H) and the power spectrum (HMAGSQR)
	FFT/OUT=1/DEST=HF  h
	FFT/OUT=4/DEST=HMAGSQR  h	
	
	//STEP 3: Make FT of the original signal to get mean power density of original signal, S
	FFT/OUT=1/DEST=SF  wave_smth
	
	//STEP 4: Make FT of the noise to get mean power density, N
	String RawHistName=waveName+"_H"					//Histogram of values in the DF/F wave
	String fitRawHistName="fit_"+RawHistName		//For fitting a Gaussian to the first part _H wave so as to estimate baseline SD 
	Make/N=0/O $RawHistName
	Histogram/B=5 $waveName,$RawHistName;DelayUpdate 
	FindPeak/Q/P $RawHistName
	Make/O/D/N=0 levelsWave
	Variable level = 0.8*waveMax($RawHistName)		//Set points for fit at either side of 0.5 the peak (avoids tail in distribution, reflecting responses)
	FindLevels/Q/N=2/D=levelsWave $RawHistName, level
	K0 = 0													//Constrain fit to baseline = 0
	Make/O/N=4 W_coef
	CurveFit/Q/X=1/H="1000"/TBOX=768 gauss $RawHistName[0,x2pnt($RawHistName,levelsWave[1])]
	noiseVariance = (W_coef[3]/2.355)^2   //FWHM = 2.355 * sd
	print W_coef[3]
	
	//STEP 5: Calculate SNR(f)
	Duplicate/O/C SF SNR
	Duplicate/O/C SF N
	N=noiseVariance
	wave curveFits
	SNR = (unitaryEvent*tau*HMAGSQR)/N   //I don't know quite why I need to scale in this way! SNR is frequency-dependent.  Gaussian power spectrum is flat

	//STEP 6: Calculate the FT of the Wiener deconvolution filter, G
	Duplicate/O/C SF G
	G = (HF^-1)*( HMAGSQR / (HMAGSQR + SNR^-1) )
		
	//STEP7: Convolution filtering of the original, xsignal, and the Wiener filter
	FFT/OUT=1/DEST=XSIG  xsignal
	Duplicate/O/C XSIG XCONVG
	XCONVG=XSIG*G
	
	IFFT/DEST=xfiltered  XCONVG
	xfiltered*=scalingFactor
	
	String smoothwaveName = waveName +"_SNew"
	String DeconWaveName = waveName +"_DNew"
	Duplicate/O xfiltered $DeconWaveName
	Duplicate/O wave_smth $smoothWaveName
	
	//	STEP 8. Find amplitudes of events found by deconvolution that cross threshold
	//////////////////////////////////////////////////////////
	//threshold = .3//*unitaryEvent
	//////////////////////////////////////////////////////////
	variable s = dimsize(xfiltered,0)
	make/o/n=(s) ampEvents, locationEvents
	variable i 
	variable j = 0
	for (i=1;i<s-1;i+=1)
		if (xfiltered[i] >xfiltered[i-1] && xfiltered[i+1] < xfiltered[i] && xfiltered[i] > threshold)
			ampEvents[j] = xfiltered[i]
			locationEvents[j] = i * deltax(xfiltered)
			j+=1
		endif
	endfor
	
	deletepoints j,s-j,ampEvents
	deletepoints j,s-j,locationEvents

	variable nEvents = dimsize(ampEvents,0)
	

	String eventswaveName=wavename + "_ENew"			//_E contains time of each event in tickswave _T
	String tickswaveName=wavename + "_T"			//Series of 1's!!!!! Simply for indicating where events have been detected when plotted vs _E
	String amplitudewaveName=wavename + "_ANew"	//_A contains amplitude of each event in tickswave _T
	String ampHistwaveName=wavename + "_AHNew"		//_AH contains amplitude histogram 

	Duplicate/O ampEvents $amplitudewaveName
	Duplicate/O locationEvents $eventswaveName
	
	display/K=1/W=(35,276,772,481) $smoothwaveName
	AppendToGraph/R $amplitudewaveName vs $eventswaveName;DelayUpdate
	SetAxis/A/E=1 right
	SetAxis/A/E=1 left
	ModifyGraph mode($amplitudewaveName)=1,rgb($amplitudewaveName)=(0,0,0);DelayUpdate
	ModifyGraph msize($amplitudewaveName)=2
	ModifyGraph rgb($amplitudewaveName)=(0,0,0)
	ModifyGraph noLabel(right)=2,axThick(right)=0	
	ModifyGraph zero(left)=4,zeroThick(left)=2
	Label left "\\F'Symbol'\\F'SymbolPi'D\\F'Helvetica'F/F";DelayUpdate
	Label bottom "Time (s)"					

	String nEventsStr= "Found "+ num2str(nEvents) +" events"
	TextBox/C/N=text0/F=0/A=LT nEventsStr

	wavestats/q $amplitudewaveName					// histogram of amplitudes and display
	variable bins=V_max/0.02
	make/o/n=(bins) $ampHistwaveName
	Histogram/B={0,0.02,bins} $amplitudewaveName,$ampHistwaveName
	display/k=1/W=(432,45,827,253)/K=1 $ampHistwaveName
	ModifyGraph mode=5,rgb=(0,0,0), hbFill=2
	Label bottom "Amplitude (ÆF/F)"
	Label left "Number of Events"
	SetAxis/Z bottom 0,V_max
	

	wave W_StatsQuantiles, W_sigma
	killwaves/Z coefs, W_findlevels,W_sigma,W_StatsQuantiles
	KillWaves/Z amps, ones, tempfit 
	KillWaves/Z ws, W_coef, levelsWave, wd, wd_smth



end
