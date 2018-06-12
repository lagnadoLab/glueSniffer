#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


///////////////////////////////////////////////////////////////////
///////////////////////// GTA GUI//////////////////////////////////
//////////////////////////////////////////////////////////////////

Window analGTA() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(600,0,950,171) as "GTA Analysis"
	Button button1,pos={17.00,14.00},size={120.00,20.00},proc=analInd,title="Analyze"
	Button button1,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button1a,pos={159.00,14.00},size={120.00,20.00},proc=compButt,title="Compile"
	Button button1a,fSize=14,fStyle=1,fColor=(65015,39905,23740)
end




/////////////////////////////////////////////////////////////////////
///////////////////////////// Compile Butt/////////////////
////////////////////////////////////////////////////////////////////
function compButt(ba) : ButtonControl
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
			compileFilters(baseName,baseFolder)
			break
		case -1:
			break
		endswitch
	return 0
end



/////////////////////////////////////////////////////////////
/////////////////////Analyze ///////////////////////////////
///////////////////////////////////////////////////////////

function analInd(ba) : ButtonControl
	Struct WMButtonAction & ba
	switch(ba.eventcode)
		case 2:
			// select events wave
			string eList = wavelist("*E",";","")
			string eventName
			prompt eventname, "Select event wave", popup, eList
			doprompt "Pick Event Wave", eventName
			if (V_flag==1)
				Abort
			Endif
			
			// Use events string to grab amps string
			string ampName = eventName[0,strlen(eventName)-2] + "AQ"
			
			//select stimwave
			string stimList = wavelist("*stim*",";","")
			string stimName
			prompt stimName, "Pick stim wave", popup, stimList
			doprompt "Pick stim wave", stimname
			if (V_flag==1)
				abort
			endif

			// CHoose name to save all analysis under
			string baseName = "GTA"
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
			
			
			// Choose stim Number and run number
			
			variable  runNum =0
			prompt runNum, "Which run?"
			doprompt "Choose run", runNum
			if (v_flag==1)
				abort
			endif
			
			qGTA(baseName,baseFolder, eventName,ampName,stimName,runNum)
			break
		case -1:
			break
		endswitch
	
	return 0
end





//////////////////////////////////////////////////////////////
/////////////////////////// qGTA ////////////////////////////
//////////////////////////////////////////////////////////////


function qGTA(baseName,baseFolder, eventName,ampName,stimName,runNum)

	string baseName, baseFolder, eventName, ampName, stimName
	variable runNum
	duplicate/o $eventName, events
	duplicate/o $ampName, amps
	duplicate/o $stimName, stims
	
	//set folder
	string newPathA = getDataFolder(1,$baseName)
	moveWave events, $newPathA
	moveWave amps, $newPathA
	setDataFolder baseName
	
	
	// params
	variable maxQ = 15
	variable stimStart = 10
	variable stimDur = 120
	variable stimStop = stimStart + stimDur
	variable drawWin = .5
	variable nEv = dimsize(events,0)
	variable dt = .001
	variable nPtsWin = drawWin/dt
	variable i,j
	make/o/n=(maxQ) nQ = 0
	killwaves/z qFilters
	
	
	for (i=1;i<maxQ+1;i+=1)
		make/o/n=(nPtsWin+1) tempFilt = 0
		setscale/p x,0,dt,tempFilt
		for (j=0;j<nEv;j+=1)
			if (events[j] > stimStart + drawWin+.001 && events[j] < stimStop)
				if (amps[j]==i)
					nQ[i-1] +=1
					duplicate/o/r=(events[j]-stimStart-drawWin,events[j]-stimStart) stims, tempStim
					setscale/p x,0,dt,tempStim
					tempFilt+=tempStim
				endif
			endif
		endfor
		duplicate/o tempFilt, blah
		concatenate/np=1 {blah}, qFilters
		killwaves/z blah
	endfor
	
	string qFiltName = "R" + num2str(runNum) + "_rawQfilt"
	string nEvFiltName = "R" + num2str(runNum) + "_qN"
	duplicate/o qFilters, $qFiltName
	duplicate/o nQ, $nEvFiltName
	killwaves/z stims
	killwaves/z qFilters, tempFilt, tempStim, nQ, events, amps
	setdatafolder baseFolder
end





////////////////////////////////////////////////////////////////
//////////////////////// Compile ///////////////////////
///////////////////////////////////////////////////////////////

function compileFilters(baseName,baseFolder)

	string baseName, baseFolder
	setdatafolder baseName
	
	string filtBase = "*" + "_rawQfilt"
	string filtList = wavelist(filtBase,";","")
	
	string nBase = "*" + "_qN"
	string nList = wavelist(nBase,";","")
	
	variable nI = itemsInList(filtList)
	variable onRun = 0
	killwaves/z nQ
	// loop over all analyzed
	do
		string filtName = stringfromlist(onRun,filtList)
		string nName = stringfromlist(onRun,nList)
		
		duplicate/o $filtName, tempFilters
		duplicate/o $nName, tempN
		if (onRun==0)
			duplicate/o tempFilters, qFilters
			duplicate/o tempN, nQ
		elseif (onRun< nI)
			duplicate/o tempFilters, tempQFilters
			qFilters+=tempQFilters
			duplicate/o tempN, tempNQ
			nQ +=tempNQ
		endif
	
		onRUn +=1
	while (onRUn < nI)
	
	killwaves/z finalQFilts
	variable maxQ = 15
	variable i
	for (i=0;i<maxQ;i+=1)
		
		duplicate/o/rmd=[][i] qFilters, tempQFilt
		if (nQ[i]>0)
		tempQFilt /= nQ[i]
		endif
		concatenate/np=1 {tempQFilt}, finalQFilts
		
	endfor
	
	
	
	matrixop/o singleFilt = sumRows(qFilters)
	variable sumEv = sum(nQ)
	singleFilt/=sumEv


	setscale/p x,-.5,.001, finalQFIlts, singleFIlt
	display/k=1 finalQFilts[][0]
	for (i=1;i<maxQ;i+=1)
		appendtograph finalQFilts[][i]
	endfor
	ModifyGraph rgb(finalQFilts)=(61166,61166,61166)
	ModifyGraph rgb(finalQFilts#1)=(56797,56797,56797)
	ModifyGraph rgb(finalQFilts#2)=(52428,52428,52428)
	ModifyGraph rgb(finalQFilts#3)=(48059,48059,48059)
	ModifyGraph rgb(finalQFilts#4)=(43690,43690,43690)
	ModifyGraph rgb(finalQFilts#5)=(39321,39321,39321)
	ModifyGraph rgb(finalQFilts#6)=(30583,30583,30583)
	ModifyGraph rgb(finalQFilts#7)=(26214,26214,26214)
	ModifyGraph rgb(finalQFilts#8)=(21845,21845,21845)
	ModifyGraph rgb(finalQFilts#9)=(17476,17476,17476)
	
	
	display/k=1 singleFIlt
	michelsonContrast(finalQFilts)
	killwaves/z tempQFilt, tempN, tempFilters, qFilters, events, amps
	end
	
	
	
	
	function michelsonContrast(finalQFilts)
	
	wave finalQFilts
	
	variable maxQ = 15
	make/o/n=(maxQ)  contrast
	make/o/n=(maxQ,2) ttp, peaks
	variable i
	
	for (i=0;i<maxQ;i+=1)
		duplicate/o/r=[][i] finalQFilts, tempFilt
		variable minP = wavemin(tempFilt)
		variable maxP = wavemax(tempFilt)
		peaks[i][0] = maxP
		peaks[i][1] = minP
		contrast[i] = (maxP-minP)/(maxP+minP)
		findlevel/r=(-.3,0) tempFilt, maxP
		variable ttpMaxT = v_levelx
		findlevel/r=(-.3,0) tempFilt, minP
		variable ttpMinP = v_levelx
		ttp[i][0] = ttpMaxT
		ttp[i][1] = ttpMinP
	endfor
	
	display/k=1 contrast
	label left "Michelson Contrast"
	label bottom "nQuanta"
	duplicate/o/rmd=[][0] ttp, ttpMaxs
	duplicate/o/rmd=[][1] ttp,ttpMins
	display/k=1 ttpMaxs
	appendtograph ttpMins
	label left "Time To Peak"
	ModifyGraph rgb(ttpMaxs)=(0,0,0)
	Legend/C/N=text0/J/A=MC "\\s(ttpMaxs) ttpMax\r\\s(ttpMins) ttpMin"
	label bottom "nQuanta"

	end