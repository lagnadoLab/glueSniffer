#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//////////////////////////////////////////////////// Quanta Stuff /////////////////////////////////////////
Function quantaStuff(ampName,timeName,stimName,protocolNum)
// Inputs are :
	// ampName is the name of the wave of glutamatergic amps  (ends in _A in glueSniffer Code)
	// timeName is the name of the wave of glutamatergic event times (_E)
	// stimName is the name of the 12 element vector wave of the order of the frequencies presented - always starts and ends with 0 (for spontaneous)
		//For example, if you presented the protocol with order [.5,1,2,5,8,10,15,20,25,30], then stimOrder would be a vector = [0,.5,1,2,5,8,10,15,20,25,30,0]
			//To make this : use this format of code
					// make/o/n=12 stimOrder
					// stimOrder[0] = {0, (insert your ten frequencies), 0}
	// Protocol number is the run number (integer) of the experiment (for example the second run of the ten stimuli)

	// Initialize some stuff
	string ampName, timeName,stimName
	variable protocolNum
	duplicate/o $ampName, amps
	duplicate/o $timeName, events
	duplicate/o $stimName, stimOrder, cycleTime
	cycleTime = 1/cycleTime// length of cycles for each stim frequency
	variable stimBlocks = 10 // number of stim blocks
	variable stimBlockTime = 6 //each block lasts this long
	variable startSA = 5 // no stim at beginning lasts this long
	variable stopSA = 5// no stim at end lasts this long
	variable nEv = dimsize(events,0) // how many events detected
	variable maxQuanta = wavemax(amps) //highes amplitude event

	// Initialize some other stuff
	make/o/n=(stimBlocks + 2) freqs, startTimes, nCyc// Initialize vectors for storing frequencies, startTimes, and number of cycles
	startTimes[0] = {0,5,11,17,23,29,35,41,47,53,59,65}// This is the start times of each stim block
	variable maxFreq = wavemax(stimOrder) // Highest frequency
	variable maxSamps = maxFreq * stimBlockTime // maximum number of cycles for any stim
	make/o/n=(maxSamps,stimBlocks+2,maxQuanta) quantaMat // Make a 3D Matrix (samples x stimulus x quanta) to store stuff
	make/o/n=(stimBlocks+2,maxQuanta) meanMat,varMat,pqMat // 2D matrix (stims x quanta) for mean and vars of nums
	make/o/n=(maxSamps,stimBlocks+2) sumMat, sumEv // 
	pqMat = 0
	
	// Set up matrices to record data - any indices of the matrix that won't be filled 
	// (since the number of cycles is different for each) set to NaN to prevent being included
	// in mean and var calculations
	quantaMat =NaN // Initiliaze both
	sumMat = NaN //
	sumEv = NaN
	make/o/n=(stimBlocks+2) numPerStim, quantPerStim
	quantPerStim = 0
	numPerStim = 0
	variable i,j,k,l
	for (i=1;i<stimBlocks+1;i+=1) //Loop over stim frequencies
		nCyc[i] = stimOrder[i] * stimBlockTime
		quantaMat[0,nCyc[i]-1][i][] = 0
		sumMat[0,nCyc[i]-1][i] = 0
		sumEv[0,nCyc[i]-1][i]=0
	endfor
	
	
	
	
	for (i=0;i<nEv;i+=1) //Loop over events
		for (j=1;j<stimBlocks+1;j+=1) //Over stim Blocks
			make/o/n=(nCyc[j]) cycStarts // find the start times of each cyce
			if (events[i] > startTimes[j] && events[i] < startTimes[j] + stimBlockTime) // Find the number of events in the entire stim block
				numPerStim[j] += 1
				quantPerStim[j]+= amps[i] // number of quanta in each block
				for (k=0;k<maxQuanta;k+=1)
					if (amps[i]==k)
						pqMat[j][k] +=1
					endif
				endfor
			endif
			for (k=0;k<nCyc[j];k+=1) // Loop over all cycles within each stim block
				cycStarts[k] = startTimes[j] + (k)*cycleTime[j]
				if (events[i] >cycStarts[k] && events[i] < cycStarts[k] +cycleTime[j]) // an event occured in the cycle
					sumMat[k][j] +=amps[i]
					sumEv[k][j]+=1
					for (l=0;l<maxQuanta;l+=1) //an event with quanta l 
						if (amps[i]==l+1)
							quantaMat[k][j][l]+=1
						endif
					endfor
				endif
			endfor
		endfor
	endfor
	make/o/n=(stimBlocks+2) meanSum,varSum, meanEv, varEv
	
	// Compute mean and vars from data
	for (i=1;i<stimBlocks+1;i+=1)
		for (k=0;k<maxQuanta;k+=1)
			duplicate/o/RMD=[0,max(0,nCyc[i]-1)][i][k] quantaMat, tempVEc
			meanMat[i][k] = mean(tempVec)
			varMat[i][k] = variance(tempVEc)
			duplicate/o/RMD=[0,max(0,nCyc[i]-1)][i] sumMat, tempSum
			meanSum[i]=mean(tempSum)
			varSum[i] = variance(tempSum)
			duplicate/o/RMD=[0,max(0,nCyc[i]-1)][i] sumEv, tempEv
			meanEv[i] = mean(tempEv)
			varEv[i] = variance(tempEv)
			
			
		endfor
	endfor
	
	duplicate/o pqMat, blahMat, sumQMat, blahMat2
	blahMat = 0
	blahMat2 = 0
	
	sumdimension/d=1/dest=sumpqMat pqMat
	for (i=0;i<stimBlocks+2;i+=1)
		for (j=0;j<maxQuanta;j+=1)
			sumQMat[i][j] *= j
			blahMat[i][j] = pqMat[i][j]/sumPqMat[i]
		endfor
	endfor
			
	sumdimension/d=1/dest = sumQMatNor sumQMat
	for (i=0;i<stimBlocks+2;i+=1)
		for (j=0;j<maxQuanta;j+=1)
			blahMat2[i][j] = sumQMat[i][j]/sumQMatNor[i]
		endfor
	endfor
	
			
			
			

	

	// Record shit in a logical manner
	string meanEventName = ampName[0,strlen(ampName)-3] +"R" + num2str(protocolnum) + "_Emat_Mean"
	string varEventName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_EMat_Var"
	string quantaMatName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qmat"
	string meanMatName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qmat_Mean"
	string varMatName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qmat_Var"
	string sumMatName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qsum"
	string sumMeanName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qsum_mean"
	string sumVarName = ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qsum_var"
	string sumPerStimName =ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_stimSum"
	string quantPerStimName = ampName[0,strlen(ampName)-3] +"R" + num2str(protocolnum) + "_stimQsum"
	string meanMatSumName =ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qmat_MeanSum"
	string varMatSumName =ampName[0,strlen(ampName)-3] + "R" + num2str(protocolnum) + "_Qmat_VarSum"
	duplicate/o numPerStim, $sumPerStimname
	duplicate/o quantPerStim, $quantPerStimName
	duplicate/o quantaMat, $quantaMatName
	duplicate/o meanMat, $meanMatName
	duplicate/o varMat, $varMatName
	duplicate/o sumMat, $sumMatName
	duplicate/o meanSum, $sumMeanName
	duplicate/o varSum, $sumVarName
	duplicate/o meanEv, $meanEventName
	duplicate/o varEv, $varEventName
	
	
	// Transform the mean and variance estimates to units of entire block
	for (i=0;i<stimBlocks+2;i+=1)
		meanMat[i] *= stimOrder[i]
		varMat[i]*=stimOrder[i]^2
	endfor
	duplicate/o meanMat, $meanMatSumName
	duplicate/o varMat, $varMatSumName
	
	
	// Clean all the shit up
	killwaves/z quantaMat, meanMat, varMat, sumMat, meanSum, varSum, numPerStim,quantPerStim
	killwaves/z cycleTime, cycStarts, amps, events, times, startTimes, nCyc, tempVec, tempEv, meanEv, varEv
	
	// Graph some stuff
	Display/k=1 $meanMatName[][0] vs stimOrder as "Mean Events Per Cycle"
	Label left "nQuantal Events Per Cycle"
	Label bottom "Frequency (Hz)"
 	i=1
 	do
		AppendToGraph $meanMatName[][i] vs StimOrder
		i+=1
	while (i<maxQuanta)
	
	Display/k=1 $varMatName[][0] vs stimOrder as "Variance Events PerCycle"
	Label left "Varianc of nQuantal Events per Cycle"
	Label bottom "Frequency (Hz)"
	i = 1
	do 
		AppendToGraph $varMatName[][i] vs StimOrder
		i+=1
	while (i<maxQuanta)
	
	
	
	Display/k=1 $meanMatSumName[][0] vs StimOrder as "Mean Events Per Stimulus Block"
	i=1
	Label Left "Mean nQuantal Events per Stim Block"
	Label bottom "Frequency (Hz)"
	do
		AppendToGraph $meanMatSumName[][i]
		i+=1
	while (i<maxQuanta)
	
	Display/k=1 $sumMeanName vs stimOrder as "Mean Quanta Per Cycle"
	Label Left "Mean Sum of Quanta Per Cycle"
	Label Bottom "Frequency"
	Display/k=1 $sumVarName vs stimOrder as "Variance Quanta Per Cycle"
	Label Left "Variance Sum of Quanta Per Cycle"
	Label Bottom "Frequency (Hz)"
end




Function tempPrec(dataName,stimOrderName,protocolNum,azNum)
string dataName, stimOrderName
variable protocolNum,azNum
duplicate/o $dataName, data
duplicate/o $stimOrderName, stimOrder
// tempPrec looks at the delay between the first spike and the cycle start
// Inputs are :
	// data is a wave of glutamatergic event times (ends in _E in glueSniffer Code)
	// stim Order is a 12 element vector wave of the order of the frequencies presented - always starts and ends with 0 (for spontaneous)
		//For example, if you presented the protocol with order [.5,1,2,5,8,10,15,20,25,30], then stimOrder would be a vector = [0,.5,1,2,5,8,10,15,20,25,30,0]
			//To make this : use this format of code
					// make/o/n=12 stimOrder
					// stimOrder[0] = {0, (insert your ten frequencies), 0}
	// Protocol number is the run number (integer) of the experiment (for example the second run of the ten stimuli)
duplicate/o stimOrder, cycleTime
cycleTime= 1/cycleTime
killwaves/z allDelays
variable nEvents= dimsize(data,0)
// Jose Code
// First 5 seconds are without stimuli, last 5 also
// Each stim block is 6 seconds
// 10 Frequencies
variable stimBlocks = 10
variable stimBlockTime = 6
variable startSA = 5
variable stopSA = 5
string delayVec = "DelayVector"

make/o/n=(stimBlocks + 2) freqs, startTimes, nCyc
startTimes[0] = {0,5,11,17,23,29,35,41,47,53,59,65}
variable maxFreq = wavemax(stimOrder)
variable maxDelays = maxFreq * stimBlockTime
make/o/n=(maxDelays,stimBlocks+2) delayMat
delayMat =NaN
variable k = 0
variable i,j
variable onWhichCyc = 0
for (i=1;i<11;i+=1) //Loop over stim frequencies
	nCyc[i] = stimOrder[i] * stimBlockTime
	make/o/n=(nCyc[i]) cycStarts
	onWhichCyc = 0
	for (j=0;j<nCyc[i];j+=1)
	k=0
		cycStarts[j] = startTimes[i] + (j)*cycleTime[i]
		
		do
			if (data[k]>cycStarts[j]&& data[k-1] < cycStarts[j])
				delayMat[onWhichCyc][i] = data[k]-cycStarts[j]
				if (delayMat[onWhichCyc][i] > cycleTime[i])
					delayMat[onWhichCyc][i] = NaN
				else
					onWhichCyc+=1
				endif
				break
			endif
			k+=1
		while (k<nEvents)
		
endfor
string delayVecType = delayVec+"_run"+num2str(protocolNum)+"Az"+num2str(azNUm)+"_stim"+ num2str(stimOrder[i])
duplicate/o/RMD=[][i] delayMat, tempVec

deletePoints onWhichCyc,(maxDelays-onWhichCyc),tempVec
duplicate/o tempVec, $delayVecType

concatenate/np {tempVec}, allDelays
string allDName = "allDelays_az"+num2str(azNum)


string delayVecHist = delayVecType + "_H"
make/o/n=100 $delayVecHist; delayUpdate
Histogram/B={0,.001,100} $delayVecType, $delayVecHist
Display/k=1 $delayVecHist
string res = num2str(dimsize($delayVecType,0)) + " of "+ num2str(nCyc[i]) +" responded"
textbox res
endfor
duplicate/o allDelays, $allDName
killwaves/z tempVec, delayMat, data, cycStarts,startTimes,cycleTime,nCyc

string allDelayHistName = "allDelays"+num2str(azNUm)+"H"
make/o/n=100 $allDelayHistName; delayupdate
histogram/b={0,.001,100} $allDName, $allDelayHistName
display/k=1 $alldelayHistName

end










////////////////// Contrast Temp Proc
Function contrastTemporal(wavename,frequency,protocolNum)
string wavename
variable frequency, protocolNum

duplicate/o $wavename, evts

killwaves/z jitter, startTimes

variable tStart = 10
variable tStop = 40
variable stimFreq = 5
variable cycLen = 1/stimFreq
variable nCyc = (tStop-tStart)/ cycLen
variable nEvents = dimsize(evts,0)

make/o/n=(nCyc) startTimes, jitter
variable i,j
variable onWhich = 0
for (i=0;i<nCyc;i+=1)
	startTimes[i] = tStart + (i*cycLen)
	for (j=1;j<nEvents;j+=1)
		if (evts[j] > startTimes[i]  && evts[j-1] <startTimes[i])
			if (evts[j] -startTimes[i] < cycLen)
				jitter[onWhich] = evts[j]-startTimes[i]
				onWhich+=1
			endif
		endif
	endfor
endfor
deletepoints onWhich,nCyc-onWhich, jitter

string jitName = (wavename[0,strlen(wavename)-3]) + "_Jitt_R" + num2str(protocolNum)
duplicate/o jitter, $jitName
killwaves/z startTimes, jitter
end




///////////////// Quant Stuff Contrast //////////////
Function contrastQuantStuff(wavename,ampName,stimFreq,protocolNum)
string ampName,waveName
variable stimFreq, protocolNum
duplicate/o $ampName, amps
duplicate/o $waveName, events

variable maxQuant = wavemax(amps) // Maximum Quanta in trace
variable nEvents = dimsize(amps,0)



// Stimulus cycle stuff
variable cycleTime = 1/stimFreq
make/o/n=(2) startStop
startStop[0] = {10,40} // Start and stop time of stim
variable nCyc = (startStop[1]-startStop[0])*stimFreq

make/o/n=(nCyc) cycStarts
make/o/n=(nCyc,maxQuant) quantaMat
quantaMat = 0
make/o/n=(nCyc) sumQuanta, eventsMat
variable i,j,k
for (i=0;i<nEvents;i+=1)
	for (j=0;j<nCyc;j+=1)
		cycStarts[j] = startStop[0] + (j*cycleTime)
		if (events[i] >cycStarts[j] && events[i] < cycStarts[j]+cycleTime)
			eventsMat[j] +=1
			sumQuanta[j]+=amps[i]
			for (k=1;k<maxQuant+1;k+=1)
				if (amps[i]==k)
					quantaMat[j][k-1] +=1
				endif
			endfor
		endif
	endfor
endfor

///Mean an var stuff for cycle stuff
make/o/n=(maxQuant) meanQ,varQ
for (i=1;i<maxQuant+1;i+=1)
	duplicate/o/rmd=[][i] quantaMat, tempVec
	meanQ[i-1] = mean(tempVec)
	varQ[i-1] = variance(tempVec)
endfor

make/o/n=1 qSumMean, qSumVar,eSumMean,eSumVar
qSumMean[0] = mean(sumQuanta)
qSumVar[0] = variance(sumQuanta)
eSumMean[0] = mean(eventsMat)
eSumVar[0] = variance(eventsMat)

string quantaMatN = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Qmat"
string quantaMean = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_QMean"
string quantaVar = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Qvar"
string sumQuantaN = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Qsum"
string sumEvents = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Esum"
string quantaSumVar = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Qsum_var"
string quantaSumMean = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Qsum_Mean"
string sumEventsMean = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Esum_mean"
string sumEventsVar = ampName[0,strlen(ampName)-3] + "R" +num2str(protocolNum) + "_Emat_Var"


duplicate/o qSumMean, $quantaSumMean
duplicate/o qSumVar, $quantaSumVar
duplicate/o eSumMean, $sumEventsMean
duplicate/o eSumVar, $sumEventsVar
duplicate/o meanQ, $quantaMean
duplicate/o varQ, $quantaVar
duplicate/o quantaMat, $quantaMatN
duplicate/o sumQuanta, $sumQuantaN
duplicate/o eventsMat, $sumEvents



killwaves/z meanQ, varQ, quantaMat, sumQuanta, eventsMat, qSumMean, qSumVar, eSumMean, eSumVar
killwaves/z amps, cycStarts, startStop, tempVec, events
end


function vectorStrengthContrast(eventsName,protocolNum)
string eventsName
variable protocolNum



variable stimStart = 10
variable stimStop = 50
variable stimPeriod = .2
string ampName = (eventsName[0,strlen(eventsName)-1])

duplicate/o $eventsName, events, stimEvents
duplicate/o $ampName, amps, stimAmps



//Extract events during stim period
variable nEvents = dimsize(events,0)
variable i 
for (i=0;i<nEvents;i+=1)
	if (events[i] < stimStart | events[i] > stimStop)
		deletepoints i,1,stimAmps,stimEvents
	endif
endfor



//Compute time since peak by modulo
variable nStimEvents = dimsize(stimEvents,0)
make/o/n=(nStimEvents) phaseTimes,cosPart, sinPart
for (i=0;i<nStimEvents;i+=1)
	phaseTimes[i] = mod(stimEvents[i]-stimStart,stimPeriod)
	cosPart[i] = cos(phaseTimes[i]*2*pi/stimPeriod)
	sinPart[i] = sin(phaseTimes[i]*2*pi/stimPeriod)
endfor
variable vectorStrength = 1/nStimEvents * sqrt(sum(sinPart)^2+sum(cosPart)^2)
variable temporalJitter = sqrt(variance(phaseTimes))


string tempString = (eventsName[0,strlen(eventsName)-3]) + "_temp_"+"R" + num2str(protocolNum)
make/o/n=2 tempValues
tempValues[0] = {vectorStrength,temporalJitter}


variable tempGuess = sqrt(2 * (1-vectorStrength)) / (2 * pi * (1/stimPeriod))

print "Estimated Temporal Jitter from Vector Strength is " + num2str(tempGuess)
print "Calculated Temporal Jitter is " + num2str(temporalJitter)
duplicate/o tempValues, $tempString
//killwaves/z events, amps, phaseTimes,cosPart,sinPart,stimAmps,stimEvents

end