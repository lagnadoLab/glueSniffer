#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.




Function threshGUI(deconMatName)
string deconMatName
duplicate/o $deconMatName, data
variable/g nROI = dimsize(data,1)
variable/g onROI = 0
variable maxDat = wavemax(data)
variable/g threshVal = 0
string/g wn = "dataWave"
duplicate/o data, $wn
string/g threshValName = (deconMatName[0,strlen(deconMatName)-10]) + "_thresh"
make/o/n=(nROI) threshVals
duplicate/o threshVals, $threshValName
MakeWindow()
end

////////////////////
Function makeWindow()

string/g wn
duplicate/o $wn, data
variable/g onROI
variable/g nEvents = 0



/// Panel Window
duplicate/o/RMD=[][onROI] data, waveDat
variable maxDat = wavemax(waveDat)

Display/K=1/W=(100,0,800,400)/N=threshyThingy waveDat as "Thresholding GUI"
ModifyGraph margin(bottom)=100
/// Slider
Slider threshSlide, pos+={25,350}, size={450,250},vert=0,value=0, proc=threshSliderProc, limits = {0,maxDat,maxDat/1000}
//Value Display
ValDisplay threshDisp value = #"threshVal"
ValDisplay nEvents value = #"nEvents"

//Next 
Button nextButt size = {100,25}, proc = nextButtProc
end




/////////// SLIDER 
Function threshSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable/g threshval = sa.curval
				doWindowShit()
							
			endif
			break
	endswitch
	return 0
End

Function doWindowShit()
string/g wn
variable/g threshVal
variable/g onROI
duplicate/o $wn, data
duplicate/o/RMD=[][onROI] data, waveDat
duplicate/o waveDat, tempWave, tempAmps, tempEvents
variable dW = dimsize(tempWave,0)
variable/g nEvents
variable i
nEvents = 0
print onROI
for (i=1;i<dW-1;i+=1)
	if (tempWave[i] > threshVal && tempWave[i] > tempWave[i-1] && tempWave[i]>tempWave[i+1])
		 nEvents +=1
		 tempAmps[nEvents] = tempWave[i]
		 tempEvents[nEvents]=pnt2x(tempWave,i)
	endif
endfor
deletePoints nEvents,dW, tempAmps
deletePoints nEvents,dW, tempEvents
removeFromGraph/z tempAmps
appendToGraph  tempAmps vs tempEvents
modifygraph mode(tempAmps)=1, rgb(tempAmps)=(0,0,0)

killwaves/z tempAmps, tempEvents
end

Function nextButtProc(ba) : ButtonControl
struct WMButtonAction&ba
	switch(ba.eventCode)
		case 2:
		string/g threshValName
		variable/g threshVal
			dowindow/k threshyThingy
			variable/g OnROI
			variable/g nROI
			duplicate/o $threshValName, threshVals
			threshVals[onROI] = threshVal
			duplicate/o threshVals, $threshValName
			killwaves/z threshVals

			if (onROI<nROI-1)
				onROI +=1
				makeWindow()
			else
				killwaves/z data, dataWave, deconMat, decWave, tempA, tempAmps, tempE,tempEvents,tempWave,thr,waveDat
			endif
			
			if(V_flag==1)
				Abort
			endif  
			

	endswitch
end
