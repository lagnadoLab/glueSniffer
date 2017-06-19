#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


///////////////////// Temporal Jitter/ Vector Strength for Frequency Tuning Protocols///////////
function timingFrequency(eventName,ampName,stimName,protocolNum)
string eventName, ampName, stimName
variable protocolNum


duplicate/o $eventName, events
duplicate/o $ampName, amps
duplicate/o $stimName, freqs


make/o/n=(12) startTimes, freqTempJitt, freqVecStren
startTimes[0] = {0,5,11,17,23,29,35,41,47,53,59,65}
variable maxGrandQ = waveMax(amps)
make/o/n=(12,maxGrandQ+1) qFreqTempJitt,qFreqVecStren, nQev

// Pull events within window
variable i,j,k,nStimEvents, tempJitt, vecStren, maxQuanta,onWhichQuanta

for (i=0;i<11;i+=1)
	pullEvents(startTimes[i],startTimes[i+1],events,amps) //Pull events within time window
	wave stimEvents,stimAmps
	nStimEvents = dimsize(stimEvents,0)
	make/o/n=(nStimEvents) tempPhaseTimes, sinPart, cosPart
	for (j=0;j<nStimEvents;j+=1)
		tempPhaseTimes[j] = mod(stimEvents[j]-startTimes[i],1/freqs[i])
		sinPart[j] = sin(tempPhaseTimes[j] * 2 * pi * freqs[i])
		cosPart[j] = cos(tempPhaseTimes[j] * 2 * pi * freqs[i])
	endfor
	tempJitt = sqrt(variance(tempPhaseTimes))
	vecStren = 1/nStimEvents * sqrt(sum(sinPart)^2 + sum(cosPart)^2)
	freqTempJitt[i] = tempJitt
	freqVecStren[i] = vecStren
	
	
	
	
	///// Quanta Timing Stuff
	
	maxQuanta = wavemax(stimAmps)
	for (k=1;k<maxQuanta+1;k+=1)
	make/o/n=(nStimEvents) qPhaseTimes, qSinPart, qCosPart
	onWhichQuanta = 0
		for (j=0;j<nStimEvents;j+=1)
			if (stimAmps[j]==k)
				qPhaseTimes[onWhichQuanta] = mod(stimEvents[j]-startTimes[i],1/freqs[i])
				qSinPart[onWhichQuanta] = sin(qPhaseTimes[onWhichQuanta]*2*pi*freqs[i])
				qCosPart[onWhichQuanta] = cos(qPhaseTimes[onWhichQuanta]*2*pi*freqs[i])
				onWhichQuanta+=1
			endif
		endfor
		deletepoints onWhichQuanta, nStimEvents-onWhichQuanta, qPhaseTimes,qSinPart,qCosPart
		qFreqTempJitt[i][k] = sqrt(variance(qPhaseTimes))
		qFreqVecStren[i][k] = 1/onWhichQuanta * sqrt(sum(qSinPart)^2+sum(qCosPart)^2)
		nQev[i][k] = onWhichQuanta
	endfor
	
				
	
			
endfor
string vecString = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+"_VS"
	string tempJittString = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+"_TJ"
	string globalTJFreq = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+"_GTJ"
	string globalVSFreq = eventName[0,strlen(eventName)-2] +"R" + num2str(protocolNum)+"_GVS"
	duplicate/o qFreqTempJitt, $tempJittString
	duplicate/o qFreqVecStren, $vecString
	duplicate/o freqTempJitt, $globalTJFreq
	duplicate/o freqVEcStren, $globalVSFreq
	killwaves/z events, amps, stimEvents,stimAmps, qFreqTempJitt,qFreqVecStren,freqTempJitt,freqVEcStren
	killwaves/z cosPart,freqs, sinPart, tempSinPart, tempPhaseTimes,tempPhaseEvents, tempCosPart
	killwaves/z startTimes,qPhaseTimes,qCosPart,qSinPart,nQevents
end


////////////////////Temporal Jitter/ Vector Strength for Contrast Protocols///////////
function timingContrast(eventName,ampName,freq,protocolNum)
string eventName, ampName
variable freq,protocolNum

duplicate/o $eventName, events
duplicate/o $ampName, amps

variable stimStart = 10
variable stimStop = 40

pullEvents(stimStart,stimStop,events,amps)
wave stimAmps, stimEvents
variable nStimevents = dimsize(stimEvents,0)
make/o/n=(nStimEvents) phaseTimes, sinPart, cosPart
variable i
for (i=0;i<nStiMEvents;i+=1)
	phaseTimes[i] = mod(stimEvents[i]-stimStart,1/freq)
	sinPart[i] = sin(phaseTimes[i] * 2 *pi * freq)
	cosPart[i] = cos(phaseTimes[i] * 2 * pi * freq)
endfor

variable tempJitt = sqrt(variance(phaseTimes))
variable vecStren = 1/nStimEvents  * sqrt(sum(sinPart)^2 + sum(cosPart)^2)


variable maxQ = wavemax(stimAmps)
make/o/n=(maxQ+1) qTJ, qVS
variable onWhichInd
variable j
for (i=0;i<maxQ+1;i+=1)
	make/o/n=(nStimEvents) qPhase,qSin,qcOS
	onWhichInd = 0
	for (j=0;j<nStimEvents;j+=1)	
		if (stimAmps[j]==i)
			qPhase[onWhichInd] = mod(stimEvents[j],1/freq)
			qSin[onWhichInd] = sin(2 * pi * qPhase[onWhichInd] * freq)
			qCos[onWhichInd] = cos(2 * pi * qPhase[onWhichInd] * freq)
			onWhichInd+=1
		endif
	endfor
	deletepoints onWhichInd, nStimEvents-onWhichInd, qPhase,qSin,qCos
	qTJ[i] = sqrt(variance(qPhase))
	qVS[i] = 1/onWhichInd * sqrt(sum(qSin)^2 + sum(qCos)^2)
	
endfor

	string vecString = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+"_VS"
	string tempJittString = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+"_TJ"
	string tempString = eventName[0,strlen(eventName)-2]  +"R" + num2str(protocolNum)+ "_temp"
	duplicate/o qTJ, $tempJittString
	duplicate/o qVS, $vecString
	make/o/n=2 tempStuff
	tempStuff = {vecStren,tempJitt}
	duplicate/o tempStuff, $tempString
	killwaves/z qPhase,qCos,qSin,qTJ,qVS,tempStuff,phaseTimes,cosPart,sinPart, stimAmps,stimEvents,events,amps
end


////////////////////PULL EVENTS CODE//////////////////////////////
//Extracts amplitudes and events within the range of [stimStart,stimStop]
function pullEvents(stimStart,stimStop,events,amps)
variable stiMStart,stimStop
wave events, amps
killwaves/z stimEvents,stimAmps
duplicate/o events, stimEvents
duplicate/o amps, stimAmps

variable nEvents = dimsize(events,0)
variable i 


for (i=1;i<nEvents;i+=1)
	if (events[i] > stimStart && events[i-1]<stimStart)
		variable startInd = i
		break
	endif
endfor


if (events[nEvents-1]<stimStop)
	variable stopInd = nEvents
else
	findlevel/q/p/edge=1 events, stimStop
	 stopInd = V_levelx -1
endif

duplicate/o/r=[startInd,stopInd] events, stimEvents
duplicate/o/r=[startInd,stopInd] amps, stimAmps
end

