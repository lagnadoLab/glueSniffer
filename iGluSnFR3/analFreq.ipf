#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



////////////////////////// Frequency GUI ////////////////////////////
/////////panel //////////

window analFreqWin(): Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(600,0,950,171) as "Frequency Analysis"
	Button button1,pos={17.00,14.00},size={120.00,20.00},proc=freqAnalButt,title="Analyze"
	Button button1,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button1a,pos={159.00,14.00},size={120.00,20.00},proc=freqGatherButt,title="Gather"
	Button button1a,fSize=14,fStyle=1,fColor=(65015,39905,23740)

end

//////////////////////////////////////////////////////////
//////////////////// gather butt ////////////////////////////
////////////////////////////////////////////////////////////////

function freqGatherButt(ba): ButtonControl
	struct wmbuttonaction & ba
	switch(ba.eventcode)
		case 2:

			string baseName
			prompt baseName, "Choose experiment Name"
			doprompt "Choose Experiment Name", baseName
			if (V_flag==1)
				abort
			endif
			string baseFolder = getDatafolder(1)
			gatherFreq(baseName, baseFolder)
			break
		case -1:
			break
		endswitch
	return 0
end
	
	
	
	
////////////////////////////////////////////////////////////////
//////////////////////////// gather freq //////////////////
///////////////////////////////////////////////////////

function gatherFreq(baseName,basefolder)
string basename,basefolder
setdatafolder baseName
killwaves/z allQperC, allEperC

//gather QperC
string qPerCS = "qPerC_R" + "*"
string qPerCList = wavelist(qperCs,";","")
//gather EperC
string epercS = "ePerC_R"+"*"
string ePerCList = wavelist(eperCs,";","")


//gather Temporal stuff
string gVSstring = "*"+ "_gVS"
string gTJstring = "*" + "_gTJ"
string qVSstring = "*" + "_qVS"
string qTJstring= "*" + "_qTJ"

string gvsList = wavelist(gvsstring,";","")
string gtjlist = wavelist(gtjstring,";","")
string qvslist = wavelist(qvsstring,";","")
string qtjlist = wavelist(qtjstring,";","")


variable nI = itemsinlist(eperCList)
make/o/n =(11,15) qTJNan = 0
make/o/n=(11,15) allQTJ,allQVS
allqtj = 0
allqvs = 0
variable onRun = 0
make/o/n=(nI) gTJNan=0


killwaves/z allGvs,allGtj
do

	string qpercN = stringfromlist(onRun,qpercList)
	string epercn = stringfromlist(onRun,eperClist)
	string gvsN= stringfromList(onRun,gvslist)
	string gtjN = stringfromList(onRun, gtjlist)
	string qvsN = stringfromlist(onRun,qvslist)
	string qtjN = stringFromList(onRun,qtjlist)
	
	duplicate/o $qperCn, qperCTemp
	duplicate/o $ePerCn, ePerCTEmp
	duplicate/o $gvsN, gvsTemp
	duplicate/o $gtjn, gtjtemp
	duplicate/o $qtjN, qtjtemp
	duplicate/o $qvsN, qvsTemp
	concatenate/np=0 {qperCTemp}, allQpERC
	concatenate/np=0 {eperCTemp}, allEpERC
	concatenate/np=1 {gvsTemp}, allGvs
	concatenate/np=1 {gtjTemp}, allGtj
	
	variable nStim = dimsize(gtjTemp,0)
	variable i,j
	for (i=0;i<nStim;i+=1)
		if (numType(gtjTemp[i])==1)
			gtjnan[i]+=1
		endif
	endfor
	
	for (i=0;i<nStim;i+=1)
		for (j=0;j<15;j+=1)
			if (numtype(qvsTemp[i][j])==0)
				qtjNan[i][j]+=1
				allQTJ[i][j]+=qtjTemp[i][j]
				allQVS[i][j]+=qvsTemp[i][j]
			endif
		endfor
	endfor
	

	
	onRun+=1
while (onRUn < nI)


matrixop/o qMean = meancols(allQperC)
matrixop/o eMean = meanCOls(allEperC)
matrixop/o qVar = varCols(allQPerC)
matrixop/o eVar = varCols(allEPerC)
duplicate/o evar, ese
duplicate/o qvar, qse
wave stimO

for (i=0;i<dimsize(qVar,1);i+=1)
	qSE[i] = sqrt(qVar[i])/sqrt(nI * 2 / stimO[i])
	eSE[i] = sqrt(eVar[i])/sqrt(nI *2 / stimO[i])
endfor
matrixop/o gTJ = allGtj^t
matrixop/o gVS = allGVS^t
killwaves/z allGtj, allGVs
matrixop/o gtjMean = meancols(gtj)
matrixop/o gVSmean = meancols(gvs)
matrixop/o gvsVar = varCols(gvs)
matrixop/o gtjVar = varCols(gtj)
make/o/n=(dimsize(gtjVar,0)) gtjse, gvsSe
for (i=0;i<dimsize(gtjVar,0);i+=1)
	if (gtjVar[i]>0)
		gtjse[i]= sqrt(gtjvar[i])/sqrt(nI)
		gvsse[i] = sqrt(gvsvar[i])/sqrt(nI)
	endif
endfor

for (i=0;i<11;i+=1)
	for (j=0;j<15;j+=1)
		if (qtjNan[i][j]>0)
			allQTJ[i][j]/=qtjNan[i][j]
		endif
	endfor
endfor

killwaves/z allgvs,allgtj, gvstemp, gtjtemp, qperCtemp, epercTemp,qvstemp,qtjtemp
end




///////////////////////////// analyzeButt//////////////////////////
function freqAnalButt(ba) : buttonControl
	struct WMButtonAction & ba
	switch(ba.eventcode)
		case 2:
			// Select Events Wave
			string eList = wavelist("*E",";","")
			string eventName
			prompt eventName, "Select Event Wave", popup, eList
			doprompt "Pick Event Wave", eventName
			if (V_flag==1)
				Abort
			Endif
			
			// Use events string to grab amps string
			string ampName = eventName[0,strlen(eventName)-2] + "AQ"
			
				// Choose stim Number and run number
			
			variable  freqNum =0
			prompt freqNum, "Which frequency protocol?"
			doprompt "Choose frequency protocol", freqNum
			if (v_flag==1)
				abort
			endif
			
			variable  runnum =0
			prompt runNUm, "Which run??"
			doprompt "Choose run?", runNum
			if (v_flag==1)
				abort
			endif
			
			variable  stimStart =10
			prompt stimstart, "Stim Start??"
			doprompt "Stim Start?", stimstart
			if (v_flag==1)
				abort
			endif
			
			
			
			// CHoose name to save all analysis under
			string baseName = "freqAnalysis"
			prompt baseName, "Choose experiment Name"
			doprompt "Choose Experiment Name", baseName
			if (V_flag==1)
				abort
			endif
			if (dataFolderExists(baseName)==0)
				print "Folder doesn't exist - making it"
				newdatafolder $baseName
			else
				print "Folder exists - adding data to it"
			endif
			string baseFolder = getDatafolder(1)
			analFreq(baseName,eventName,ampname,freqNum,runNum,stimStart,basefolder)
			vstjGlobalFreq(eventname, ampname,basename,freqNum,runNum,stimstart,basefolder)
			break
		case -1:
			break
		endswitch
	
	return 0
end
			

			



////////////////////////////////////////////////////////////////////////////////


function analFreq(basename,eventName, ampName, freqNum,runNum,stimstart,basefolder)

	//declare inputs
	string baseName,eventName, ampName
	variable freqNum, runNum,stimstart
	string basefolder

	// duplicate events
	duplicate/o $eventName, events
	duplicate/o $ampName, amps
	
	///// set folder
	string newPathA = getDataFolder(1,$baseName)
	movewave events, $newPatha
	movewave amps, $newPatha
	setdataFolder baseName
	
	//pull protocol Stuff
	variable stimDur = frequencyProtocols(freqNum,stimstart)
	wave stimO
	wave startTimes
	variable nStims = 12

	duplicate/o stimO, stimInd
	stimInd = p
	
	//sort stimuli by frequency
	sort stimO, stimO,startTimes, stimInd



	// intialize some things
	make/o/n=(nStims) nCyc
	variable i,j


	
for (i=0;i<nStims;i+=1)
		nCyc[i] = stimDur * stimO[i] // get nCyc vector
	endfor

	variable maxCyc = wavemax(nCyc)

	make/o/n=(maxCyc,nStims) qPerC = NaN 
	duplicate/o qPerC, ePerC

	//count quanta/cycle and events/cycle for each stimulus frequency
	for (i=0;i<nStims;i+=1)
		if (stimO[i] > 0)
		for (j=0;j<maxCyc;j+=1)
			qPerC[j][i] = countFreqAmps(startTimes[i]+((j)/stimO[i]) ,startTimes[i]+((j+1)/stimO[i]),events,amps)
			ePerC[j][i] = countFreqEvs(startTimes[i]+((j)/stimO[i]),startTimes[i]+((j+1)/stimO[i]),events,amps)
				if (j>=nCyc[i])
					qPerC[j][i] = NaN
					ePerC[j][i] = NaN
				endif
		endfor
		endif
	endfor

	//define how to svae shit
	string qPerCS = "qPerC_R" + num2str(runNum) + "_stimNum" + num2str(freqNum)
	string ePerCs = "ePerC_R" + num2str(runNUm) + "_stimNum" + num2str(freqNum)

	//save shit
	duplicate/o qPerC, $qperCs
	duplicate/o eperC, $EpercS
	
	//clean stuff up
	killwaves/z qperC,eperC,ncyc,protocolMat,events,amps


	/////////// Move back to base folder
	setdatafolder baseFolder
	

end


//////////////////////////////////////////////////////
/////////////////////////// Get Protocol Things ///////////
////////////////////////////////////////////////////////////


function frequencyProtocols(protocolNum,startTime)

variable protocolNum,startTime

make/o/n=12 stimO


if (protocolNum==1)
	stimO[0] = {0,0.5,1,2,5,8,10,15,20,25,30,0}
elseif (protocolNum==2)
	stimO[0] = {0,8,15,10,0.5,1,20,2,30,5,25,0}
elseif (protocolNum==3)
	stimO[0]={0,15,0.5,8,2,10,1,20,30,5,25,0}
elseif (protocolNum==4)
	stimO[0]={0,5,2,1,0.5,30,25,20,15,10,8,0}
elseif (protocolNum==5)
	stimO={0,30,25,20,15,10,8,5,2,1,0.5,0}
endif


make/o/n=1 zero = 0
make/o/n=(11) tempStart
tempStart = {0,6,12,18,24,30,36,42,48,54,60}
killwaves/z startTimes
tempStart+=startTIme
concatenate/np=0 {zero,tempStart}, startTimes




variable stimDur = 6
killwaves/z zero, tempStart
return stimDur

end




/////////////////////////////////////////////////////////////////////
////////////////////// Pull Freq EVENTS ////////////////////////////
////////////////////////////////////////////////////////
function countFreqEvs(stimStart,stimStop,events,amps)

variable stimStart,stimStop
wave events,amps

variable nEv = dimsize(events,0)

variable i
variable nStimEvs = 0

for (i=0;i<nEv;i+=1)
	if (events[i]>stimStart && events[i] < stimStop)
		nStimEvs+=1
	endif
endfor
return nStimEvs
end



///////////////////////////////////////////////////
//////////////////// Pull Freq AMPS /////////////
///////////////////////////////////////////////////


function countFreqAmps(stimStart,stimStop,events,amps)

variable stimStart,stimStop
wave events,amps

variable nEv = dimsize(events,0)

variable i
variable nStimAMps = 0

for (i=0;i<nEv;i+=1)
	if (events[i]>stimStart && events[i] < stimStop)
		nStimAmps+=amps[i]
	endif
endfor
return nStimAmps
end

/////////////////////////////////////////////////////////////////////
//////////////////////// Vector Strength/Temp for Freq/////////////
/////////////////////////////////////////////////////////////////////
function vstjGlobalFreq(eventname, ampname,basename, freqNum,runNum,startTIme,basefolder)

	//declare inputs
	string eventname, ampname,basename
 	variable freqNum,runNum,startTIme
 	string basefolder
 	
 	// make events from events wave
	duplicate/o $eventName, events
 	duplicate/o $ampName, amps
 	
	///// set folder
	string newPathA = getDataFolder(1,$baseName)
	movewave events, $newPatha
	movewave amps, $newPatha
	setdataFolder baseName

 	
//pull protocol Stuff
variable stimDur = frequencyProtocols(freqNum,startTime)
wave stimO
wave startTimes
duplicate/o stimO, stimInd

sort stimO, stimO, startTimes, stimInd
	
	// declare some variables
	variable nStim = dimsize(stimO,0)
	variable maxQ = 15
	
	// initialize vectors for storage
	make/o/n=(nStim) nEvs, vecStr,tempJitt
	
	variable i,j,k
	for (i=0;i<nStim;i+=1)
		nEvs[i] = countFreqEvs(startTimes[i],startTimes[i]+stimDur,events,amps) //count events in range
		make/o/n=(nEvs[i]) tempPhase,sinPart,cosPart //initialize a wave for all temp data
		pullEvs(startTimes[i],startTimes[i]+stimDur,events,amps) // pull events and amps from range
		wave stimEvents, stimAmps //declare the events and amps from range
		if (dimsize(stimEvents,0)!=nEvs[i]) //make sure shit didn't fuck up
			print "stimEvents erro"
		endif
		if (nEvs[i] > 0) // must have at least 2 points for this
			for (j=0;j<nEvs[i];j+=1) // cycle over events in range
				tempPhase[j] = mod(stimEvents[j]-startTimes[i],1/stimO[i]) // find their phase in s
				sinPart[j] = sin(tempPhase[j]* 2 * pi * stimO[i]) // find sin part
				cosPart[j] = cos(tempPhase[j] * 2 * pi * stimO[i]) // and cos part
			endfor
			vecStr[i] = 1/nEvs[i] * sqrt(sum(sinPart)^2+ sum(cosPart)^2) //compute vs
			tempJitt[i] = sqrt(2 *(1-vecStr[i])/(2*pi*stimO[i])) //use formula to compute tj
		else //if less than 2, can't do analysis
			vecStr[i] = NaN
			tempJitt[i] = NaN
		endif
		
		make/o/n=(nStim,maxQ) qTJ, qVS, nQEv
		
		
		///quantal part
		if (nEvs[i]>0)
		for (j=0;j<maxQ;j+=1)
			variable onWhichQ = 0
			make/o/n=(nEvs[i]) qPhase, qCos, qSin
			for (k=0;k<nEvs[i];k+=1)
				if (stimAmps[k]==j+1)
					qPhase[onWhichQ] = mod(stimEvents[k]-startTimes[i],1/stimO[i]) //phase
					qSin[onWhichQ] = sin(qPhase[onWhichQ]* 2 * pi * stimO[i])//cose part\
					qCos[onWhichQ] = cos(qPhase[onWhichQ] * 2 * pi * stimO[i])//sin prt
					onWhichQ+=1
				endif
			endfor
			deletepoints onWhichQ, nEvs[i]-onWhichQ, qPhase,qSin,qCos
			qVS[i][j] = 1/onWhichQ * sqrt(sum(qSin)^2 + sum(qCos)^2)
			qTJ[i][j] = sqrt(2 * (1-qVS[i][j]))/(2 * pi *stimO[i])
			nQEv[i][j] = onWhichQ
			endfor
		else
			qVS[i][j] = NaN
			qTJ[i][j] = NaN
		endif

	endfor
	
	
	for (i=0;i<nStim;i+=1)
		for (j=0;j<maxQ;j+=1)
			if (nQev[i][j] < 2)
				qVS[i][j] = NaN
				qTJ[i][j] = NaN
			endif
		endfor
	endfor
	
	
	//define how to save stuff
	string gVSstring = "Run_" + num2str(runNum) + "_gVS"
	string gTJstring = "Run_" + num2str(runNum) + "_gTJ"
	string qVSstring = "Run_" + num2str(runNum) + "_qVS"
	string qTJstring= "Run_" + num2str(runNum) + "_qTJ"
	string nqevS = "Run_" + num2str(runNum) + "_qnEv"
	string nEvsS = "Run_" + num2str(runNum) + "_nEv"
	
	
	//save stuff
	duplicate/o NEvs, $nEvss
	duplicate/o qVS, $qvsstring
	duplicate/o nqev, $nqevS
	duplicate/o qtj, $qtjstring
	duplicate/o vecStr, $gVSstring
	duplicate/o tempjitt, $gTJstring
	
	
	//clean up
	killwaves/z vecStr, tempJitt, sinPart, cosPart, qvs, qtj,nqev,qphase,qsin,qcos,events, amps,stimevents,stimamps,tempPhase
	
	
	
	
	
	
	
	
	/////////// Move back to base folder
	setdatafolder baseFolder
	

end





//////////////////////////////////////////////////////

function pullEvs(start,stop,events,amps)
	variable start,stop
	wave events, amps
	variable nEvents = dimsize(events,0) // count all events
	
	variable i
	variable nStimEvs=0
	
	
	///count events during stimulus
	for (i=0;i<nEvents;i+=1)
		if (events[i] > start && events[i]< stop)
			nStimEvs+=1
		endif
	endfor
	
	if (nStimEvs>0)
		for (i=0;i<nEvents;i+=1)
			if (events[i]>start && events[i-1] < start) 
				variable startInd = i
				break
			endif
		endfor
		
		if (events[nEvents-1] < stop)
			// last event occured within stim
			variable stopInd = nEvents
		else
			for (i=0;i<nEvents;i+=1)
				if (events[i]> stop && events[i-1] < stop)
					 stopInd = i-1
				endif
			endfor
		endif
		
		duplicate/o/r=[startInd,stopInd] events, stimEvents
		duplicate/o/r=[startind,stopInd] amps, stimAmps
	else
		make/o/n=0 stimEvents,stimAmps	
	endif
		
		
		
	end
			
	
	
	
////////////////////////////////////////////////////////////
////////////////////// section evenst by freqyency


function sectionEvents(eventname, ampname,basename, freqNum,runNum,startTIme,basefolder)

	//declare inputs
	string eventname, ampname,basename
 	variable freqNum,runNum,startTIme
 	string basefolder
 	
 	// make events from events wave
	duplicate/o $eventName, events
 	duplicate/o $ampName, amps
 	
	///// set folder
	string newPathA = getDataFolder(1,$baseName)
	movewave events, $newPatha
	movewave amps, $newPatha
	setdataFolder baseName

 	
//pull protocol Stuff
variable stimDur = frequencyProtocols(freqNum,startTime)
wave stimO
wave startTimes
duplicate/o stimO, stimInd

sort stimO, stimO, startTimes, stimInd
	
	// declare some variables
	variable nStim = dimsize(stimO,0)
	variable maxQ = 15
	variable i,j,k
	
	make/o/n=(nStim,maxQ) nQEvs

	
	///Extract phase times for each quantal event	
	for (i=0;i<nStim;i+=1)
		for (j=0;j<maxQ;j+=1)
			nQevs[i][i] = pullSpecificEvents(startTimes[i],startTimes[i]+stimDur,j+1,events,amps)
			wave specEvents
			duplicate/o specEvents, phase
			
			for (k=0;k<nQEvs[i][j];K=1)
				phase[k] = mod(specEvents[k]-startTimes[i],1/stimO[i])
			endfor
			
		endfor
	endfor
end
	
	
	
	
	
end

	
		//////////////////////////////////////////////////////////
	//////////////////////////specific events ///////////////
	///////////////////////////////////////////////////////////
	
function pullSpecificEvents(startTime,stopTime,specAmp,events,amps)
	variable startTIme,stopTIme,specAmp
	wave events, amps
	killwaves/z specEvents
	make/o/n=1 tempAmpTime
	
	variable nEvs = dimsize(events,0) // number of events
	variable nSpec=0
	variable i,f
	for (i=0;i<nEvs;i+=1)
		if (amps[i]==specAmp)
			nSpec+=1
			tempAmpTime=events[i]
		endif
		concatenate/np {tempAmpTime}, specEvents
	endfor
end
			