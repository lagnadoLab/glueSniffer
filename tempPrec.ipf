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
	
			
			
			
	
	display blahMat[1][]
	for (i=2;i<stimBlocks+1;i+=1)
		appendTograph blahMat[i][]
	endfor
	
	display blahMat2[1][]
	for (i=2;i<stimBlocks+1;i+=1)
		appendTograph blahMat2[i][]
	endfor
	
	

	// Record shit in a logical manner
	string meanEventName = ampName[0,strlen(ampName)-3] + "Emat_Mean"
	string varEventName = ampName[0,strlen(ampName)-3] + "EMat_Var"
	string quantaMatName = ampName[0,strlen(ampName)-3] + "Qmat"
	string meanMatName = ampName[0,strlen(ampName)-3] + "Qmat_Mean"
	string varMatName = ampName[0,strlen(ampName)-3] + "Qmat_Var"
	string sumMatName = ampName[0,strlen(ampName)-3] + "Qsum"
	string sumMeanName = ampName[0,strlen(ampName)-3] + "Qsum_mean"
	string sumVarName = ampName[0,strlen(ampName)-3] + "Qsum_var"
	string sumPerStimName =ampName[0,strlen(ampName)-3] + "stimSum"
	string quantPerStimName = ampName[0,strlen(ampName)-3] +"stimQsum"
	string meanMatSumName =ampName[0,strlen(ampName)-3] + "Qmat_MeanSum"
	string varMatSumName =ampName[0,strlen(ampName)-3] + "Qmat_VarSum"
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
 	i=1
 	do
		AppendToGraph $meanMatName[][i] vs StimOrder
		i+=1
	while (i<maxQuanta)
	
	Display/k=1 $varMatName[][0] vs stimOrder as "Variance Events PerCycle"
	i = 1
	do 
		AppendToGraph $varMatName[][i] vs StimOrder
		i+=1
	while (i<maxQuanta)
	
	
	
	Display/k=1 $meanMatSumName[][0] vs StimOrder as "Mean Events Per Stimulus Block"
	i=1
	do
		AppendToGraph $meanMatSumName[][i]
		i+=1
	while (i<maxQuanta)
	
	Display/k=1 $sumMeanName vs stimOrder as "Mean Quanta Per Cycle"
	Display/k=1 $sumVarName vs stimOrder as "Variance Quanta Per Cycle"
end




Function tempPrec(dataName,stimOrderName,protocolNum)
string dataName, stimOrderName
variable protocolNum
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
		cycStarts[j] = startTimes[i] + (j)*cycleTime[i]
		
		do
			if (data[k]>cycStarts[j])
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
string delayVecType = delayVec+"_stimRun"+num2str(protocolNum)+"_stim"+ num2str(stimOrder[i])
duplicate/o/RMD=[][i] delayMat, tempVec

deletePoints onWhichCyc,(maxDelays-onWhichCyc),tempVec
duplicate/o tempVec, $delayVecType
endfor
killwaves/z tempVec, delayMat, data, cycStarts,startTimes,cycleTime,nCyc



end
