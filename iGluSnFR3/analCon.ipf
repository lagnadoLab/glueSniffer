#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



//////////////////////////////////////////////////////////////
//////////////////// GUI for contrast experiment analysis///////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////Panel ////////////////

Window analCon() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(600,0,950,171) as "Contrast Analysis"
	Button button1,pos={17.00,14.00},size={120.00,20.00},proc=simpleContrastButt,title="Simple Contrast"
	Button button1,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button1a,pos={159.00,14.00},size={120.00,20.00},proc=gatherSimpleButt,title="Simple Gather"
	Button button1a,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button2, pos={17.00,140.00},size={120.00,20.00}, proc = complexContrastButt, title = "Complex Contrast"
	Button button2,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button3, pos={159.00,140.00},size={120.00,20.00}, proc = gatherComplexButt, title = "Complex Gather"
	Button button3,fSize=14,fStyle=1,fColor=(65015,39905,23740)
end








/////////////////////////////////////////////////////////////////////
///////////////////////////// Gather complex Butt/////////////////
////////////////////////////////////////////////////////////////////
function gatherComplexButt(ba) : ButtonControl
	Struct WMButtonAction & ba
	switch(ba.eventcode)
		case 2:

			string baseName
			prompt baseName, "Choose experiment Name"
			doprompt "Choose Experiment Name", baseName
			if (V_flag==1)
				abort
			endif
			string baseFolder = getDatafolder(1)
			gatherComplex(baseName,baseFolder)
			break
		case -1:
			break
		endswitch
	return 0
end




/////////////////////////////////////////////////////////////////////
///////////////////////////// Gather Simple Butt/////////////////
////////////////////////////////////////////////////////////////////
function gatherSimpleButt(ba) : ButtonControl
	Struct WMButtonAction & ba
	switch(ba.eventcode)
		case 2:

			string baseName
			prompt baseName, "Choose experiment Name"
			doprompt "Choose Experiment Name", baseName
			if (V_flag==1)
				abort
			endif
			string baseFolder = getDatafolder(1)
			gatherSimple(baseName, baseFolder)
			break
		case -1:
			break
		endswitch
	return 0
end

//////////////////////////////////////////////////////////////////
//////////////////////// SImple Contrast /////////////////////////
//////////////////////////////////////////////////////////////////

Function simpleContrastButt(ba): ButtonControl
	Struct WMButtonAction & ba
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
			
			variable  contrastNum =0
			prompt contrastNum, "Which Contrast?"
			doprompt "Choose Contrast", contrastNum
			if (v_flag==1)
				abort
			endif
			
			
			// CHoose name to save all analysis under
			string baseName = "complexC"
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
			analSimple(baseName,eventName,ampName,contrastNum,baseFolder)
			break
		case -1:
			break
		endswitch
	
	return 0
end
			
			
			
//////////////////////////////////////////////////////////
////////////////////// Anal Simple ///////////////////////
//////////////////////////////////////////////////////////
function analSimple(baseName,eventName,ampName,contrastNum,baseFolder)

	string basename, eventname, ampname
	variable contrastNum
	string baseFolder
	duplicate/o $eventname, events
	duplicate/o $ampName, amps

	///// set folder
	string newPathA = getDataFolder(1,$baseName)
	movewave events, $newPatha
	movewave amps, $newPatha
	setdataFolder baseName
	
	// Set some hard parameters
	variable maxQ = 15
	variable nEv = dimsize(events,0)



	///Protocol params

	variable stimStart = 10
	variable stimStop =40
	variable freq = 5
	variable nCyc = (stimStop-stimStart)*freq

	make/o/n=(nCyc) qPerC, ePerC
	qPerC = 0
	ePerC = 0

	variable i,j, cycStart

	for (i=0;i<nCyc;i+=1)
		cycStart = stimStart + (i/freq)
		qPerC[i] = qCounts(cycStart,cycStart+(1/freq), events,amps)
		ePerC[i] = eCounts(cycStart,cycStart+(1/freq), events,amps)
	endfor
			
	

	make/o/n=(maxQ) nnQ

	for (i=0;i<nEv;i+=1)
		for (j=0;j<maxQ;j+=1)
			if (events[i] > stimStart && events[i] < stimStop && amps[i]==j)
				nnQ[j-1]+=1
			endif
		endfor
	endfor


	/// Find minimum iei
	make/o/n=1 minIeiVal
	minIeiVal[0] = findMinIEI(events)
	
	
	spontaneousStuff(events,amps,contrastNum)
	
	
	/////////////////////////Define how to save shit ///////////////////////////////
	string ePerCName = "C" + num2str(contrastNum) + "_cycE"
	string qPerCname = "C" + num2str(contrastNum) + "_cycQ"
	string saveOutNnq = "C" + num2str(contrastNum) + "_nnQ"

	string minIEIout = "C" + num2str(contrastNum) + "_minIEI"
	duplicate/o ePerC, $ePerCname
	duplicate/o qPerC, $qPerCname
	duplicate/o nnq, $saveOutNnq
	duplicate/o minIEIVal, $minIEIout
	
		/////////////// Kill the leftovers /////////////
	killwaves/z qPerC, ePerC,nnQ,stimO, sortInd, tempE, tempQ,tempNnq,sortedEperC,sortedQperC,sortedNnq
	killwaves/z events, amps, minIEIVal


/////////// Move back to base folder
	setdatafolder baseFolder
	
end









//////////////////////////////////////////////////////////////
/////////////////////////////////Gather Simple////////////////
/////////////////////////////////////////////////////////////////


function gatherSimple(baseName, baseFolder)
	string baseName, baseFolder
	setdatafolder baseName
	killwaves/z allCycQ, allCycE,allNnq, minIEIs
	
	string ecBase = "*" + "cycE"
	string ecList = wavelist(ecBase,";","")
	
	string qcBase = "*" + "cycQ"
	string qcList = wavelist(qcBase, ";","")
	
	string nnqBase = "*" + "nnQ"
	string nnQList = wavelist(nnqBase,";","")
	
	string minIEIbase = "*" + "_minIEI"
	string minIEIlist = wavelist(minIEIBase,";","")
	
	string saBase = "*" + "_SA"
	string saList = wavelist(saBase,";","")
	
	///loop over all runs
	variable nI = itemsinlist(ecList)
	variable onRun = 0
	
	
	variable contrastN
	string cN
			
	make/o/n=(nI) contrasts

	
	
	onRun = 0
	make/o/n=2 meanSA
	meanSA = 0
	killwaves/z allCycQ, allCycE, allNnq allIEI,sortContInd
	do
	
	
	
		string ecName = stringfromlist(onRun, ecList)
		string qcName = stringfromlist(onRun, qcList)
		string nnqName = stringfromlist(onRun, nnqList)
		string minIEIName = stringfromlist(onRun, minIEIbase)
		string saName = stringfromlist(onRun, saList)

		
		cN = ecName[1,strlen(ecName)-6]
		contrastN = str2num(cN)
		contrasts[onRun] = contrastN
		duplicate/o $ecName, ec
		duplicate/o $qcName, qc
		duplicate/o $nnqName, nnQ
		duplicate/o $minIEIname, minieiTemp
		duplicate/o $saName, sa
		
		concatenate/np=1 {ec}, allCycE
		concatenate/np=1 {qc}, allCycQ
		concatenate/np=1 {nnq},allNnq
		concatenate/np=1 {minieiTemp}, allIEI
		meanSA+=sa
		
		onRUn +=1
	while (onRUn < nI)
	
	
	/// sort into order
	variable i
	duplicate/o contrasts, sortContInd
	sortContInd = p
	sort contrasts, contrasts, sortContInd
	
	
	killwaves/z tempE, tempQ, sortedEperC, sortedQperC, tempNnq, sortedNnq
	
	for (i=0;i<nI;i+=1)
		duplicate/o/rmd=[][sortContInd[i]] allCycE, tempE
		duplicate/o/rmd=[][sortContInd[i]] allCycQ, tempQ
		duplicate/o/rmd=[][sortContInd[i]] allNnq,tempNnq
		concatenate/np=1 {tempE}, sortedEperC
		concatenate/np=1 {tempQ},sortedQperC
		concatenate/np=1 {tempNnq}, sortedNnq
	endfor
	
	killwaves/z test,tempQ,tempNnq,tempE,qc, nnQ,iei,allCycQ,allNnq,allCycE,ec
	
	
	duplicate/o sortedEperC, allCycE
	duplicate/o sortedQperC, allCycQ
	duplicate/o sortedNnq, allNnq
	
	killwaves/z sortedEperC, sortedQperC, sortedNnq, sortContInd
	
	variable binSizeVal = .02
	
	//compute probabilities matrices
	findpnQ(ALLNnQ,nI,binSizeVal,30) //p(nQ)
	string pnqName = "pnQ"
	findpSums(allCycQ,allCycE) //p(sumQ) and p(sumE)
	string pSumEname = "pSumE" 
	string psumQname = "pSumQ"
	
	removezeroRows(pSumEname)
	removezerorows(pSumQname)
	removezerorows(pnqName)
	DUPLICATE/O $pSumEname, pSumE
	duplicate/o $pSumQname, pSumQ
	duplicate/o $pnqName, pnQ
	//Information Theory
	mutualinfo(pnqName)
	mutualinfo(psumeName)
	mutualinfo(pSumQName)
	specInf(pnqName)
	specInf(pSumEName)
	specInf(pSumQName)
	
	string pnqi2N = "pnQ_i2"
	string pnqpr = pnqi2N + "_pR"
	
	wave pr = $pnqpr
	entropy(pr)
	expVes(pr,binSizeVal)
	killwaves/z pr
	
	duplicate/o $pnqi2N, pnQi2perves
	deletepoints 0,1, pnqi2perves
	duplicate/o, pnqi2perves, i2pervesp
	i2pervesp = p+1
	pnqi2perves/=i2pervesp
	string pnqi2PerVesn = "pnQ_i2_perVes"
	duplicate/o pnqi2perves, $pnqi2perVesN
 
	killwaves/z pnqi2PerVes, i2pervesp
	
	
	DUPLICATE/O $pSumEname, pSumE
	duplicate/o $pSumQname, pSumQ
	duplicate/o $pnqName, pnQ
	//Information Theory
	mutualinfo(pnqName)
	mutualinfo(psumeName)
	mutualinfo(pSumQName)
	specInf(pnqName)
	specInf(pSumEName)
	specInf(pSumQName)
	
	 pnqi2N = "pnQ_i2"
	 pnqpr = pnqi2N + "_pR"
	
	wave pr = $pnqpr
	entropy(pr)
	expVes(pr,binSizeVal)
	killwaves/z pr
	
	duplicate/o $pnqi2N, pnQi2perves
	deletepoints 0,1, pnqi2perves
	duplicate/o, pnqi2perves, i2pervesp
	i2pervesp = p+1
	pnqi2perves/=i2pervesp
	 pnqi2PerVesn = "pnQ_i2_perVes"
	duplicate/o pnqi2perves, $pnqi2perVesN
 
	killwaves/z pnqi2PerVes, i2pervesp,outQ,outE, ptst,pts,pt, sumSs, sumRs,pMat2
	
	
	
	
	//display
	killwaves/z meanQs, meanEs, varQ,varE, newQse, newEse,eSe,qSE
	matrixop meanQs = meancols(allCycQ)
	matrixop meanEs = meanCols(allCycE)
	matrixop varQ = varcols(allCycQ)
	matrixop varE = varcols(allCycE)

	
	make/o/n=(nI) qSe, eSE
	for (i=0;i<nI;i+=1)
		qSE[i] = sqrt(varQ[i])/sqrt(nI)
		eSE[i] = sqrt(varE[i])/sqrt(nI)
	endfor
		
		
	insertPoints 0,1, meanQs, meanEs, contrasts, qSE, eSE
	meanQs[0] = meanSa[1]
	meanEs[0] = meanSa[0]
	qSE[0] = 0
	eSE[0] = 0
	contrasts[0] = 0
	matrixop newQse = qSE^t
	matrixop newEse = ESE^t
	
	
	display meanQs vs contrasts
	appendtograph meanEs vs contrasts
	ErrorBars meanQs Y,wave=(newQse,newQse)
	ErrorBars meanEs Y,wave=(newEse,newEse)
	
	
	killwaves/z sa, outS, events,amps
end













///////////////////////////////////////////////////////////////
//////////////////// Complex Contrast////////////////////////////
///////////////////////////////////////////////////////////////
Function complexContrastButt(ba) : ButtonControl
	Struct WMButtonAction & ba
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
			variable stimNum=1
			variable  runNum=1
			prompt stimNum, "Which Stimulus Protocol?"
			doprompt "Choose Stimulus Protocol", stimNum
			if (v_flag==1)
				abort
			endif
			
			prompt runNum, "Which Run?"
			doprompt "Choose Run Number", runNum
			if (V_flag==1)
				abort
			endif
			
			//Choose 0 contrast
			variable zeroContrast = 0
			prompt zeroContrast, "Select base contrast"
			doprompt "Select Base Contrast", zeroContrast
			if (V_flag==1)
				abort
			endif
			
			
			// CHoose name to save all analysis under
			string baseName = "complexC"
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
			analComplex(baseName,eventName,ampName,runNum,stimNum,zeroContrast,baseFolder)
			break
		case -1:
			break
		endswitch
	
	return 0
end
		


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////Analyze Complex Protocol //////////////////////
/////////////////////////////////////////////////////////////////////////////////

function analComplex(baseName,eventName,ampName,runNum,stimNum,zeroContrast,baseFOlder)
// Inputs : baseName - string that all results will be saved as - makes collecting easier after
//			 : evetnName - string of event wve
//			 : ampName - string of clustered amplitude wave
//			 : runNum - run of the protocol
//        : stimNum - which stimulus set did you run
// Outputs : cycE - nStim (rows) x nCyc (cols) of events per cycle
//			  : cycQ - nStim (rows) x nCyc (cols) of quanta per cycle
//			  : nnQ - nStim (rows) x maxQ (cols) of n(nQuanta) 
//			  : outOrder - sorted stim order
// 		  : minIei - minimum iei value

	string baseName,eventName, ampName, baseFolder
	variable runNum,stimNum,zeroContrast
	
	
	duplicate/o $eventName, events
	duplicate/o $ampName, amps
	
	
	///// set folder
	string newPathA = getDataFolder(1,$baseName)
	movewave events, $newPatha
	movewave amps, $newPatha
	setdataFolder baseName
	
	// Set some hard parameters
	variable maxQ = 15
	variable nEv = dimsize(events,0)
	
	//gather protocol results
	getParams(stimNum,runNum)
	string rawOrder = "R" + num2str(runNum) + "_rawOrder"
	string paramOut = "R" + num2str(runNum) + "_params"
	duplicate/o $rawOrder, stimO
	duplicate/o $paramOut, params
	variable stimStart = params[0]
	variable nTrials = params[1]
	variable isi = params[2]
	variable stimLen = params[3]
	variable freq = params[4]
	
	///Initialize stuff for later storing data
	variable nCycPerTrial = freq * stimLen
	make/o/n=(nCycPerTrial, nTrials) qPerC, ePerC
	qPerC = 0
	ePerC = 0
	variable cycStart
	variable tStart
	
	
	//////////////////////// Univariate Measures ////////////////////////////////
	// Compute quanta and events per cycle for each contrast
	variable i,j
	for (i=0;i<nTrials;i+=1)
		for (j=0;j<nCycPerTrial;j+=1)
			cycStart = stimStart + (i * (isi+stimLen)) + (j/freq)
			qPerC[j][i] = qCounts(cycStart,cycStart+(1/freq), events,amps)
			ePerC[j][i] = eCounts(cycStart,cycStart+(1/freq), events,amps)
			print qPerC[j][i]
		endfor
	endfor
	
	//Find number of each nQuanta for each stim
	make/o/n=(maxQ,nTrials) nnQ
	nnQ = 0
	variable k
	for (i=0;i<nTrials;i+=1)
		tStart = stimStart + (i * (isi+stimLen))
		for (j=0;j<nEv;j+=1)
			for (k=1;k<maxQ+1;k+=1)
				if (events[j] > tStart && events[j] < tStart + stimLen && amps[j]==k)
					nnQ[k-1][i] +=1
				endif
			endfor
		endfor
	endfor

	
	/// Find minimum iei
	make/o/n=1 minIeiVal
	minIeiVal[0] = findMinIEI(events)
	
	
	
	/////////////////////////Multivariate Measures///////////////////////////////
	
	
	
	/////////////////////// Sort back into order //////////////////////////////
	stimO += zeroContrast
	make/o/n=(nTrials) sortInd
	sortInd = p
	sort stimO,stimO,sortInd
	
	
	killwaves/z tempE, tempQ, sortedEperC, sortedQperC, tempNnq, sortedNnq
	for (i=0;i<nTrials;i+=1)
		duplicate/o/rmd=[][sortInd[i]] ePerC, tempE
		duplicate/o/rmd=[][sortInd[i]] qPerC,tempQ
		duplicate/o/rmd=[][sortInd[i]] nnQ,tempNnq
		concatenate/np=1 {tempE}, sortedEperC
		concatenate/np=1 {tempQ},sortedQperC
		concatenate/np=1 {tempNnq}, sortedNnq
	endfor
	setscale/p y,0,1, sortedNNq, sortedEperC,sortedQperC
	
	
	
	
	/////////////////////////Define how to save shit ///////////////////////////////
	string ePerCName = "R" + num2str(runNum) + "_cycE"
	string qPerCname = "R" + num2str(runNum) + "_cycQ"
	string saveOutNnq = "R" + num2str(runNum) + "_nnQ"
	string saveOut = "R" + num2str(runNum) + "_outOrder"
	string minIEIout = "R" + num2str(runNum) + "_minIEI"
	duplicate/o sortedEperC, $ePerCname
	duplicate/o sortedQperC, $qPerCname
	duplicate/o sortedNnq, $saveOutNnq
	duplicate/o stimO, $saveOut
	duplicate/o minIEIVal, $minIEIout
	
	
	/////////////// Kill the leftovers /////////////
	killwaves/z qPerC, ePerC,nnQ,stimO, sortInd, tempE, tempQ,tempNnq,sortedEperC,sortedQperC,sortedNnq
	killwaves/z events, amps, minIEIVal

	/////////// Move back to base folder
	setdatafolder baseFolder
end




//////////////////////////////////////////////////////////////
/////////////////////////////////Gather Complex////////////////
/////////////////////////////////////////////////////////////////

// This function takes the outputs of analCOmplex (for each run of the protocols)
// and combines the data
// INPUTS : 
//			 : basename - the name of the folder data from analComplex is stored in
//			 : baseFolder - the root folder for the analysis
// OUTPUTS : 
//			 : allCycQ - nStim (row) x nCycle (column) off quanta per cycle
//        : allCycE - nStime (row) x nCycle (column) of events per cycle
//			 : nnq - nStim (row) x maqQ (column) of number of nQ events found in each stimBlock
//			 : pnQ -  nStim (row) x maxQ+1 (column) of the number of probability of finding
//					each event type (including zero events) in a window
//			 : i2 - maxQuanta + 1 vector of specific information between each event type and the stim
//			 : i2PerVes - Same as i2, but divided by the number of vesicles in each event
//			 : mutInfo - mutual information

function gatherComplex(baseName, baseFolder)
	string baseName, baseFolder
	setdatafolder baseName
	killwaves/z allCycQ, allCycE,allNnq, minIEIs
	
	string ecBase = "*" + "cycE"
	string ecList = wavelist(ecBase,";","")
	
	string qcBase = "*" + "cycQ"
	string qcList = wavelist(qcBase, ";","")
	
	string nnqBase = "*" + "nnQ"
	string nnQList = wavelist(nnqBase,";","")
	
	string minIEIbase = "*" + "_minIEI"
	string minIEIlist = wavelist(minIEIBase,";","")
	
	
	///loop over all runs
	variable nI = itemsinlist(ecList)
	variable onRun = 0
	
	
	
	do
		string ecName = stringfromlist(onRun, ecList)
		string qcName = stringfromlist(onRun, qcList)
		string nnqName = stringfromlist(onRun, nnqList)
		string minIEIName = stringfromlist(onRun, minIEIbase)
		
		duplicate/o $ecName, ec
		duplicate/o $qcName, qc
		duplicate/o $nnqName, nnQ
		duplicate/o $minIEIname, minieiTemp

		
		concatenate/np=0 {ec}, allCycE
		concatenate/np=0 {qc}, allCycQ
		concatenate/np=0 {minIEItemp}, minIeis
		if (onRun==0)
			duplicate/o nnQ, allNnq
		else
			allNNq += nnq
		endif
		onRun+=1
		
			
	while (onRun < nI)
	
	
	
	matrixop/o qMean = meancols(allCycQ)
	matrixop/o varQ = varCols(allCycQ)
	matrixop/o eMean = meanCols(allCycE)
	matrixop/o varE = varCols(allCycE)
	
	matrixop/o eVar = varE^t
	matrixop/o qVar = varQ^t
	variable nCyc = dimsize(allCycQ,0)
	variable rootN = sqrt(nCYc)
	variable nVars = dimsize(varQ,1)
	//comupte se
	make/o/n=(nVars) qSe, eSe
	variable i
	for (i=0;i<nVars;i+=1)
		qSe[i] = sqrt(qvar[i])/rootN
		eSe[i] = sqrt(evar[i])/rootN
	endfor
	
		
	
	
	
	
	
	///////////// INformation Theory Stuff//////////////////////////
	//find bin size for single variable metrics
	make/o/n=1 globalMinIei
	globalMinIEI[0] = wavemin(minIeis)
	variable binSizeVal = findBinSize(globalMinIEI[0])
	make/o/n=1 binSizeParam
	binSizeParam[0] = binSizeVal
	
	
	//find triple joint
	tripleJoint(allCycQ,allCycE)
	wave pQESdis
	//compute probabilities matrices
	findpnQ(ALLNnQ,nI,binSizeVal,2) //p(nQ)
	string pnqName = "pnQ"
	findpSums(allCycQ,allCycE) //p(sumQ) and p(sumE)
	string pSumEname = "pSumE" 
	string psumQname = "pSumQ"
	
	removezeroRows(pSumEname)
	removezerorows(pSumQname)
	removezerorows(pnqName)
	
	DUPLICATE/O $pSumEname, pSumE
	duplicate/o $pSumQname, pSumQ
	duplicate/o $pnqName, pnQ
	//Information Theory
	mutualinfo(pnqName)
	mutualinfo(psumeName)
	mutualinfo(pSumQName)
	specInf(pnqName)
	specInf(pSumEName)
	specInf(pSumQName)
	tripleMI(pQESdis)
	condMutInf(pQESdis)
	string pnqi2N = "pnQ_i2"
	string pnqpr = pnqi2N + "_pR"
	
	wave pr = $pnqpr
	entropy(pr)
	expVes(pr,binSizeVal)
	killwaves/z pr
	
	duplicate/o $pnqi2N, pnQi2perves
	deletepoints 0,1, pnqi2perves
	duplicate/o, pnqi2perves, i2pervesp
	i2pervesp = p+1
	pnqi2perves/=i2pervesp
	string pnqi2PerVesn = "pnQ_i2_perVes"
	duplicate/o pnqi2perves, $pnqi2perVesN
 
	killwaves/z pnqi2PerVes, i2pervesp
	
	
	DUPLICATE/O $pSumEname, pSumE
	duplicate/o $pSumQname, pSumQ
	duplicate/o $pnqName, pnQ
	//Information Theory
	mutualinfo(pnqName)
	mutualinfo(psumeName)
	mutualinfo(pSumQName)
	specInf(pnqName)
	specInf(pSumEName)
	specInf(pSumQName)
	
	 pnqi2N = "pnQ_i2"
	 pnqpr = pnqi2N + "_pR"
	
	wave pr = $pnqpr
	entropy(pr)
	expVes(pr,binSizeVal)
	killwaves/z pr
	
	duplicate/o $pnqi2N, pnQi2perves
	deletepoints 0,1, pnqi2perves
	duplicate/o, pnqi2perves, i2pervesp
	i2pervesp = p+1
	pnqi2perves/=i2pervesp
	 pnqi2PerVesn = "pnQ_i2_perVes"
	duplicate/o pnqi2perves, $pnqi2perVesN
 
	killwaves/z pnqi2PerVes, i2pervesp
	
	

end



end


///////////////////////////////////////////////////////
////////////////////////// Find p(sumQ) ////////////
/////////////////////////////////////////////////////
function FindpSums(allCycQ,allCycE)
wave allCycQ, allCycE

variable nCyc = dimsize(allCycQ,0)
variable nStim = dimsize(allCycQ,1)

variable maxQSum = wavemax(allCycQ)
variable maxW = maxQSum + 1

make/o/n=(maxW,nStim) pSumQ, pSumE

variable i,j,k
pSumQ = 0
pSuME = 0
for (i=0;i<nCYc;i+=1)
	for (j=0;j<nStim;j+=1)
		for (k=0;k<maxW;k+=1)
			if (allCycQ[i][j]==k)
				pSumQ[k][j]+=1
			endif
			if (allCycE[i][j]==k)
				pSumE[k][j]+=1
			endif
		endfor
	endfor
endfor
pSumE/=nCyc
pSumE/=nStim
pSumQ/=nCyc
pSumQ/=nStim


duplicate/o pSume, outE
duplicate/o pSumQ, outQ
killwaves/z pSume,psumq

string pSumEname = "pSumE"
string pSumQName = "pSumQ"
duplicate/o outE, $pSumEName
duplicate/o ouTQ,$pSumQName
end


/////////////////////////////////////////////////////
/////////////////////////Find pnQ ////////////////
/////////////////////////////////////////////////
function findpnQ(allNNq,nI,binSize,stimLen)
wave allNNq
variable nI, binSize
variable stimLen

variable nBins = nI * stimLen/binSize
variable nStims = dimsize(allnnQ,1)
make/o/n=(nStims) nZero
nZero=nBins
matrixop/o nqColSum = sumcols(allNNq)
nZero-=nqColSum

killwaves/z pnQ,pnq22

matrixop/o nnqT = allNNq^t
matrixop/o nZt = nZero^t


concatenate/np=0 {nZt, allNNQ}, pnq22

pnq22 /= nBins 
pnq22/=nStims
duplicate/o pnq22, pnQ
killwaves/z nZero, nZt, nnqT, nqColsum, pnq22



end

//////////////////////////////////////////////////////////
//////////////////

function removeZeroRows(pmat2Name)
string pmat2Name
duplicate/o $pmat2Name, pMat2

matrixop/o pt = pMat2^t
matrixop/o pts = sumcols(pt)
matrixop/o ptst = pts^t
findlevel/p/q ptst,0
if (v_flag<1)
	deletepoints/m=0 v_levelx, dimsize(pts,1)-v_levelx,pmat2
endif
duplicate/o pmat2,$pmat2Name
end

///////////////////////////////////////////////////////////////////////
//////////////////////////// Find min IEI /////////////////////
//////////////////////////////////////////////////////////////

function findMinIEI(events)
wave events
variable nEv = dimsize(events,0)
variable i
make/o/n=(nEv) ieid

ieid[0] = events[0]
for (i=1;i<nEv;i+=1)
	ieid[i] = ieid[i]-ieid[i-1]
endfor
variable minIEI = wavemin(ieid)
return minIEI

end
	

///////////////////////////////////////////////////////////////
//////////////////////// Select bin size ///////////////////////
///////////////////////////////////////////////////////////////

function findBinSize(minIeiVal)
variable minIeiVal
variable binSize = .02

do
	if (binSize>minIeiVal)
		binSize/=2
	endif
while (binSize > minIeiVal)

return binSize
end


	
	/////////////////////////////////////////////////////////////////////////
	////////////////////////// Protocol Stuff //////////////////
	//////////////////////////////////////////////////////////////////////
function getParams(stimNum,runNum)
	variable stimNum,runNum
	variable stimStart, nTrials, isi, stimLen, freq
	if (stimNum < 3)
		stimStart = 5
		nTrials = 10
		isi = 2
		stimLen = 2
		freq = 5
	elseif (stimNum==6)
		stimStart = 5
		nTrials = 11
		isi = 2
		stimLen = 2
		freq = 5
	else
		stimStart = 9
		nTrials = 11
		isi = 2
		stimLen = 2
		freq = 5
	
	endif
	
	make/o/n=5 params
	params = {stimStart, nTrials, isi, stimLen, freq}
	
	make/o/n=(nTrials) stimO
	if (stimNum==1) // protocol A
		stimO ={2,14,6,18,10,12,4,16,8,20}
	elseif (stimNum==2) 
		stimO= {10,6,8,2,4,20,16,18,12,14}
	elseif (stimNum==3)
		stimO = {0,2,14,6,18,10,12,4,16,8,20}
	elseif (stimNum==4)
		stimO= {10,6,8,2,4,20,16,18,12,14,0}
	elseif (stimNum==5)
		stimO = {50,20,60,40,80,30,70,10,90,0,100}
	elseif (stimNum==6)
		stimO = {50,20,60,40,80,30,70,10,90,0,100}
	endif
	
	string rawOrder = "R" + num2str(runNum) + "_rawOrder"
	string paramOut = "R" + num2str(runNum) + "_params"
	duplicate/o stimO, $rawOrder
	duplicate/o params, $paramOut
	
	killwaves/z params, stimO
end




//////////////////////////////////////////////////////////////////////////////////	
/////////////////////////////////Count Quanta from each Cycle/////////////////////////
//////////////////////////////////////////////////////////////////////////////////	

function qCounts(winStart,winStop,events,amps)
	variable winStart, winStop
	wave events, amps
	
	variable qInCyc = 0
	
	variable nEv = dimsize(events,0)
	variable i
	
	for (i=0;i<nEv;i+=1)
		if (events[i] > winStart && events[i] < winStop)
			qInCyc+= amps[i]
		endif
	endfor
	
	return qInCyc
end



//////////////////////////////////////////////////////////////////////////////////	
/////////////////////////////Count Events from Each Cycle///////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//////////////////////////////////////////////////////////////////////////////////	

function eCounts(winStart,winStop,events,amps)
	variable winStart, winStop
	wave events, amps
	
	variable eInCyc = 0
	
	variable nEv = dimsize(events,0)
	variable i
	
	for (i=0;i<nEv;i+=1)
		if (events[i] > winStart && events[i] < winStop)
			eInCyc+= 1
		endif
	endfor
	
	return eInCyc
end



function spontaneousStuff(events,amps,contrastNum)

wave events,amps
variable contrastNum

variable stiMStart = 10
variable rollStart = 5
variable cycLen = .2
variable dt= .001

variable nSamps = (stiMstart-rollStart)/dt -(cycLen/dt)

make/o/n=(nSamps) spontQperC, spontEperC
variable i, cycStart

for (i=0;i<nSamps;i+=1)
	cycStart = rollStart + (i*dt)
	spontQPerC[i] = qCounts(cycStart,cycStart+cycLen,events,amps)
	spontEPerC[i] = eCounts(cycStart,cycStart+cycLen,events,amps)
endfor

variable meanSpontQ = mean(spontQperC)
variable meanSpontE = mean(spontEperC)

string outSpont = "C" + num2str(contrastNum) + "_SA"
make/o/n=2 outS
outS = {meanSpontE,meanSpontQ}
duplicate/o outS, $outSpont

killwaves/z  spontQperC, spontEPerC
end







///////////////////////////////////////////////// Find p(Q,E,S)
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

function tripleJoint(allCycQ,allCycE)

wave allCycQ, allCycE



variable nCyc = dimsize(allCycQ,0)
variable nStim = dimsize(allCycQ,1)

variable maxQ = wavemax(allCycQ)+1
variable maxEv = wavemax(allCycE)+1


make/o/n=(maxQ,nStim,maxEv) pQES

variable i,j,k,l


pQES = 0



for (i=0;i<nCyc;i+=1)
	for (j=0;j<nStim;j+=1)
		for (k=0;k<maxQ;k+=1)
			for (l=0;l<maxEv;l+=1)
				if (allCycQ[i][j]==k && allCycE[i][j]==l)
					pQES[k][j][l]+=1
				endif
			endfor
		endfor
	endfor
endfor
variable sumP = sum(pQES)
pQES/=sumP



string QESname = "pQESdis"

duplicate/o pQES, $QESname

killwaves/z pQES
end


end