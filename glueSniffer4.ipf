#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "Sarfia"
#include "Advanced ROI tools"
#include "Linescan"
#include "matrix"
#include "EventDetectionnew"
#include "parseoutse"
#include "Ch2LineRes"
#include "expDiff"
#include "emGMM"
#include "testGMM"
#include "tempPrec"
#include "threshGUI"
#include "timingStuff"
#include "analCon"
#include "infoTheory"
#include "statTests"
#include "lifModel"
#include "gtaAnal"
#include "analFreq"


////////////////////////////////////////// Panel /////////////////////////////////////////////////////
Window glueSniffer4() : Panel
	PauseUpdate; Silent 1		// buildiang window...
	NewPanel /W=(2551,-135,2852,171) as "Events Analysis New"
	ModifyPanel cbRGB=(11454,30202,51877)
	SetDrawLayer UserBack
	SetDrawEnv linethick= 0,fillfgc= (64824,27308,21496)
	DrawRRect 5,50,290,205
	Button button1,pos={17.00,14.00},size={120.00,20.00},proc=loadMoviesButt,title="Load Movies"
	Button button1,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button1a,pos={159.00,14.00},size={120.00,20.00},proc=enterParamsButt,title="Enter Params"
	Button button1a,fSize=14,fStyle=1,fColor=(65015,39905,23740)
	Button button2,pos={17.00,56.00},size={120.00,20.00},proc=splitChannelButt,title="Split channels"
	Button button2,fColor=(65015,39905,23740)
	Button button3,pos={159.00,57.00},size={120.00,20.00},proc=lineScanButt,title="Linescan"
	Button button3,fColor=(65015,39905,23740)
	Button button4,pos={16.00,87.00},size={120.00,20.00},proc=azButt,title="New Define AZs"
	Button button4,fColor=(65015,39905,23740)
	Button button5,pos={159.00,89.00},size={120.00,20.00},proc=azButtold,title="Old Define AZs"
	Button button5,fColor=(65015,39905,23740)
	Button button6,pos={17.00,120.00},size={120.00,20.00},proc=temporalProfileButt,title="New Temporal Profile"
	Button button6,fColor=(65015,39905,23740)
	Button button7,pos={159.00,120.00},size={120.00,20.00},proc=temporalProfileButtOld,title="Old Temporal Profile"
	Button button7,fColor=(65015,39905,23740)
	Button button12,pos={17.00,150.00},size={120.00,20.00},proc=threshButt,title="Threshold"
	Button button12,fColor=(65015,39905,23740)
	Button button13,pos={159.00,150.00},size={120.00,20.00},proc=eventDetButt,title="Event Detection"
	Button button13,fColor=(65015,39905,23740)
	Button button9,pos={17.00,180.00},size={120.00,20.00},proc=clusterButt,title="Cluster"
	Button button9,fColor=(65015,39905,23740)
	
	Button button20,pos = {17,230}, size = {120,20}, proc = contrastAnalButt, title ="Contrast"
	Button button20,fColor=(65015,39905,23740)
	
	Button button21, pos = {17,260}, size = {120,20}, proc = LIFButt, title ="LIF"
	Button button21,fColor=(65015,39905,23740)
	
	Button button22, pos = {159,230}, size = {120,20}, proc = freqButt, title ="Freq"
	Button button22,fColor=(65015,39905,23740)
	
	Button button23, pos = {159,260}, size = {120,20}, proc = GTA, title = "GTA"
	Button button23,fColor=(65015,39905,23740)
	Button button24, pos = {17,290}, size = {120,20}, proc = IEI, title = "IEI"
	Button button24,fColor=(65015,39905,23740)
	Button button25, pos = {159,290}, size = {120,20}, proc = batch, title = "batch"
	Button button25, fColor=(65015,39905,23740)
end














/////////////////////////////////////////////
////////////////////////////////////////// Load Movies ////////////////////////////////////////
Function batch(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			Variable refNum
			String message = "Select One or More Files"
			String outputPaths
			String fileFilters = "Data Files (*.tif):.tif;"
			fileFilters +="All Files:.*;"
			Open /D /R /Mult=1 /F=fileFilters /M=message refNum
			outputPaths = S_fileName
			
			String filenamestr, loadedFilenameStr, FluorFilenameStr, LSFilenameStr, DFFileNameStr
			Variable filenameStart
			String pathstr
			Variable nChannels=2, stimChannel=2
			if (strlen(outputPaths) == 0)
				Print "Cancelled"
			else
				Variable numFilesSelected = ItemsInList(outputPaths, "\r")
				Variable i
				for(i=0; i<numFilesSelected; i+=1)
					String path = StringFromList(i, outputPaths, "\r")
					filenameStart = strsearch(path, ":", strlen(path), 1)
					filenamestr=path[filenameStart+1, strlen(path)]
					pathstr=path[0,filenameStart-1]
					AutoLoadScanImage_L(pathstr, filenamestr)
					loadedFilenameStr=filenamestr[0, strlen(filenamestr)-5]   //remove the ".tif" at the end of the filename
					splitchannels_L(loadedFilenameStr, 2, 2)
					string fName = loadedFilenamestr + "_F"
					LSimage($fName)
					string lName = loadedFilenamestr + "_L"
					newDefine(lName)
					string muName = (lName[0, (strlen(lName)-3)])+"_muHat"
					string AName = (lName[0, (strlen(lName)-3)])+"_AHat"
					string sigmaName = (lName[0,strlen(lName)-3])+"_sigmaHat"
					TemporalProfile(lName,aName,muName,sigmaName)
					string dmName = loadedFileNameStr + "_deconMat"
					threshGui(dmName)
					duplicate/o $dmName, deconMat
					variable nROIs = dimsize(deconMat,1)
					variable nPoints = dimsize(deconMat,0)
					string threshName = (dmName[0,strlen(dmName)-10]) + "_thresh"
					duplicate/o $threshName, thr
					string baseName = (dmName[0,strlen(dmName)-10]) + "_AZ"
					string e1 = basename + "1_E"
					string a1 = basename + "1_A"
					//wave ampMat,evMat
					variable roi
					variable nEvent = 0
					for (roi=0;roi<nROIs;roi+=1)
				duplicate/o/RMD=[][roi] deconMat, decWave, tempA,tempE
				tempE = 0
					tempA = 0
					nEvent = 0
				for (i=1;i<nPoints-1;i+=1)
					
					if (decWave[i] > decWave[i-1] && decWave[i]>decWave[i+1] && decWave[i] > thr[roi])
						tempA[nEvent] = decWave[i]
						tempE[nEvent] = pnt2x(decWave,i)
						nEvent+=1
					endif
				endfor
				if (roi==0)
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a1
					duplicate/o tempE, $e1
				elseif (roi==1)
					string e2 = basename + "2_E"
					string a2 = basename + "2_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a2
					duplicate/o tempE, $e2
				elseif (roi==2)
					string e3 = basename + "3_E"
					string a3 = basename + "3_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a3
					duplicate/o tempE, $e3
				elseif (roi==3)
					string e4 = basename + "4_E"
					string a4 = basename + "4_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a4
					duplicate/o tempE, $e4
				elseif (roi ==4)
					string e5 = basename + "5_E"
					string a5 = basename + "5_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a5
					duplicate/o tempE, $e5
				elseif (roi==5)
					string e6 = basename + "6_E"
					string a6 = basename + "6_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a6
					duplicate/o tempE, $e6
				endif
			endfor	
		
		
		killwaves/z decWave,deconMat,data,aHat, decWave, tempE,tempA
				endfor
			endif	
			break
		case -1:
			break
		endswitch
		
	return 0
End



end





//////////////////////////////////////////////////////////////////

Function GTA(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			string executeGTA = "analGTA()"
			execute executeGTA
		break
		case -1:
		break
	endswitch
	return 0
end


///////////////////////////////////////////////////////////////////////////
Function LIFButt(ba) : ButtonControl
	Struct WMButtonAction &ba
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
			lifModel(eventName,ampName)
				case -1:
			break
		endswitch
	
	return 0
end
			
///////////////////// Contrast Analysis Button /////////////////
Function freqButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			string executeFreqAnal = "analFreqWin()"
			execute executeFreqAnal
			break
		case -1:
			break
	endswitch
	return 0
	
End


///////////////////// Contrast Analysis Button /////////////////
Function contrastAnalButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			string executeConAnal = "analFreqWin)"
			execute executeConAnal
			break
		case -1:
			break
	endswitch
	return 0
	
End

	

////////////////////////////////////////// Load Movies ////////////////////////////////////////
Function loadMoviesButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			Variable refNum
			String message = "Select One or More Files"
			String outputPaths
			String fileFilters = "Data Files (*.tif):.tif;"
			fileFilters +="All Files:.*;"
			Open /D /R /Mult=1 /F=fileFilters /M=message refNum
			outputPaths = S_fileName
			
			String filenamestr, loadedFilenameStr, FluorFilenameStr, LSFilenameStr, DFFileNameStr
			Variable filenameStart
			String pathstr
			Variable nChannels=2, stimChannel=2
			if (strlen(outputPaths) == 0)
				Print "Cancelled"
			else
				Variable numFilesSelected = ItemsInList(outputPaths, "\r")
				Variable i
				for(i=0; i<numFilesSelected; i+=1)
					String path = StringFromList(i, outputPaths, "\r")
					filenameStart = strsearch(path, ":", strlen(path), 1)
					filenamestr=path[filenameStart+1, strlen(path)]
					pathstr=path[0,filenameStart-1]
					AutoLoadScanImage_L(pathstr, filenamestr)
					loadedFilenameStr=filenamestr[0, strlen(filenamestr)-5]   //remove the ".tif" at the end of the filename
				endfor
			endif	
			break
		case -1:
			break
		endswitch
		
	return 0
End



/////////////////////Enter Params///////////////
Function enterParamsButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			Variable tauRiseTemp, tauTemp, thresholdTemp
			
			
			tauRiseTemp = 0.001
			prompt tauRiseTemp, "Enter Tau Rise"
			doprompt "Enter Tau Rise", tauRiseTemp
			if(V_flag==1)
					Abort
			endif
					case -1: // control being killed
			
			
			tauTemp = 0.06
			prompt tauTemp, "Enter Tau"
			doprompt "Enter Tau", tauTemp
			if(V_flag==1)
					Abort
			endif
					case -1: // control being killed
			

			
			Variable /G tauRise = tauRiseTemp
			Variable /G tau=tauTemp 
			break
	endswitch
	return 0	
end


	////////////////////////////////////////// Split Channels ////////////////////////////////////////
Function splitChannelButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			string list=wavelist("*",";","DIMS:3")
			string movieName
			variable nChannels, stimChannel
			
			prompt movieName, "Wave to split", popup,list
			doprompt "Pick your movie", MovieName
			if(V_flag==1)
					Abort
			endif
			
			prompt nChannels, "How many channels?"
			doprompt "Choose the number of channels", nChannels
			if(V_flag==1)	
					Abort
			endif
			
			prompt stimChannel, "Which channel contains stimulus (1, 2, 3 or 4)?"
			doprompt "Choose the number of channels", stimChannel
			if(V_flag==1)	
					Abort
			endif
			
			splitchannels_L(movieName, nChannels, stimChannel)
			break
		
		case -1: // control being killed
			break
	endswitch
	return 0	
end






//////////////////////////////////////////////////////////////////////////////////
Function SplitChannels_L(waveName, nChannels, stimChannel)// split channels
		
	string waveName
	Variable nChannels
	Variable stimChannel

	Duplicate/O $waveName PicWave
	variable nFrames = DimSize(PicWave,2), FramesPerChannel, Rest, ii
	String wvName
	FramesPerChannel=nFrames/nChannels
	Rest=FramesPerChannel-trunc(FramesPerChannel)
	
	if(Rest)
		Print "WARNING: unequal number of frames per channel."
		FramesPerChannel=trunc(FramesPerChannel)
	endif
	
	For(ii=0;ii<nChannels;ii+=1)
		if(ii!=(stimChannel-1))
			wvName=NameOfWave($waveName)+"_F"
			Duplicate /o PicWave $wvName
			wave w=$wvName
			Redimension/n=(-1,-1,FramesPerChannel) w
			MultiThread w=PicWave[p][q][r*nChannels+ii]
		else
			wvName=NameOfWave($waveName)+"_stim"
			Duplicate /o PicWave $wvName
			wave w=$wvName
			Redimension/n=(-1,-1,FramesPerChannel) w
			MultiThread w=PicWave[p][q][r*nChannels+ii]	
		endif			
	EndFor
	KillWaves/Z PicWave
	
	return 0
			
end
//////////////////////////////////////////////////////////////////////////////




Function /wave AutoLoadScanImage_L(pathstr, filenamestr)			//Does not open a file dialogue, takes path and filename as input parameters
	String pathstr, filenamestr

	String ImgWaveName, FirstWave
	string header, s_info = "No header info available\r"
	Variable PointPos
	
	NewPath/o/q path, pathstr
	
	ImageLoad /Q /O /C=-1/p=path filenamestr
	
	if (v_flag == 0)
		Abort
	endif
	
		header = s_info
		PointPos = strsearch(S_Filename, ".tif", 0)
		ImgWaveName = S_FileName[0,PointPos-1]
		ImgWaveName = ReplaceString("-", ImgWaveName, "_")
		
		PointPos = strsearch(S_Wavenames, ";", 0)
		FirstWave =S_Wavenames[0,PointPos-1]
		
	if (waveexists($ImgWaveName))
		killwaves /z $ImgWaveName
	endif
	
	
	duplicate /o $FirstWave, $ImgWaveName
	Killwaves /z $FirstWave
	
	redimension /d $ImgWaveName		//convert to double precision floating point
	
	Note $ImgWaveName, header
	Note $ImgWaveName, "file.path="+s_path
	Note $ImgWaveName, "file.name="+s_filename
	
	Wave ReturnWv =  $ImgWaveName
	Return ReturnWv
End

//////////////////////////////////////////////////////


////////////////////////////////////////// Linescan ////////////////////////////////////////
Function linescanButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			string list=wavelist("*_F",";","DIMS:3")
			string name
			prompt name, "wave to analyse", popup,list
			doprompt "pick your movie", name
			
			LSimage($name)
			
			if(V_flag==1)	
					Abort
			endif
			
		break
		
		case -1: // control being killed
			break
	endswitch

	return 0
End

/////////////////////////////////////////////////////////////////////////////////

Function newDefine(wavename)
	string wavename
	variable autoAbortSecs = 0
	variable mu1,mu2,mu3,mu4,mu5,mu6
	variable A1,A2,A3,A4,A5,A6
	string muName = (waveName[0, (strlen(waveName)-3)])+"_muHat"
	string AName = (waveName[0, (strlen(waveName)-3)])+"_AHat"
	string out=(waveName[0, (strlen(waveName)-3)])+"_profile"
	string coefName = (waveName[0,strlen(waveName)-3])+"_Coef"
	string fitName =(waveName[0,strlen(waveName)-3])+"_fit"
	string spatFilt = (waveName[0,strlen(waveName)-3])+"_SF"
	string sigmaHatName = (waveName[0,strlen(waveName)-3])+"_sigmaHat"
	
	duplicate/o $waveName $out
	duplicate/o $waveName w1
         
   make/o/n=(dimsize(w1,1)) w2
   matrixop/o w2=sumCols(w1)
   make/o/n=(dimsize(w2,1)) w3
   w3=w2[0][p]
   variable sm = mean(w3)
   variable sw = WAVEMAX(w3)
   //w3-=sm
   w3/=sw
   variable i
   for (i=0;i<dimsize(w3,0);i+=1)
   	if (w3[i]<0)
   		w3[i] =0
   	endif
   endfor
   duplicate/o w3 $out
   
   
   display/k=1/N=Profile $out
   Legend/C/N=text0/J/F=0/A=RT "\\f01\\Z14Define AZ means with cursors and close window"
   ShowInfo/CP={0,1,2}/W=Profile
   MoveWindow/W=Profile 0, 400, 600, 600
   
   
   ///////////Cursor Input ///////////////////////////////////////////
   if (UserCursorAdjust("Profile",autoAbortSecs)!=0)
   	return -1
   endif
   
   mu1 = pcsr(A)
   A1= vcsr(A)
   /// CSR 2, if there
   string csrResult
   csrResult=CsrInfo(B,"Profile")
   if (cmpstr(csrResult,"")==0)
   	mu2 = 0
   	A2 = 0
   else
   	mu2 = pcsr(B)
   	A2 = vCsr(B)
   endif
   
   //Csr 3, if there
   csrResult=CsrInfo(C,"Profile")
    if (cmpstr(csrResult,"")==0)
   	mu3 = 0
   	A3 = 0
   else
   	mu3 = pcsr(C)
   	A3 = vCsr(C)
   endif
   
   //Csr 4, if there
   csrResult=CsrInfo(D,"Profile")
    if (cmpstr(csrResult,"")==0)
   	mu4 = 0
   	A4 = 0
   else
   	mu4 = pcsr(D)
   	A4 = vCsr(D)
   endif
   
   //Csr 5, if there
   csrResult=CsrInfo(E,"Profile")
    if (cmpstr(csrResult,"")==0)
   	mu5 = 0
   	A5 = 0
   else
   	mu5 = pcsr(E)
   	A5 = vCsr(E)
   endif
   
   //Csr 6, if there
   csrResult=CsrInfo(E,"Profile")
    if (cmpstr(csrResult,"")==0)
   	mu6 = 0
   	A6 =  0
   else
   	mu6 = pcsr(E)
   	A6 = vcsr(E)
   endif
   
   make/o/n=6 tmpMu, tmpA
   tmpMu = {mu1,mu2,mu3,mu4,mu5,mu6}
   tmpA = {A1,A2,A3,A4,A5,A6}
   findValue/z/v=0 tmpMu
   variable firstZero = V_Value
   print firstZero
   deletePoints firstZero, 6-firstZero, tmpMu
   deletePoints firstZero, 6-firstZero, tmpA
	DoWindow/K Profile
	Killwaves w1,w2,w3
	duplicate/o tmpMu, $muName
	duplicate/o tmpA, $AName
	killwaves tmpMu, tmpA
	fitGauss(out,muName,aName)
	Display/k=1/n=fitGraph $out as "Gaussian Mixture Fit"
	print fitName
	print spatFilt
	AppendToGraph $fitName
	Label left "Normalized Intensity";DelayUpdate
	Label bottom "Location (pixels)"
	ModifyGraph rgb($fitName)=(0,0,0)
	Legend/C/N=text1/A=MC

	
	
   end
////////////////////////////////////////// Define AZs ////////////////////////////////////////
Function azButt(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Detect AZs")
			string list=wavelist("*_L",";","DIMS:2")
			string waveName
			prompt waveName, "Line scan to analyse", popup,list
			doprompt "Pick your Line scan", waveName
			newDefine(waveName)	   	
		break
				
		case -1: // control being killed
			break
	endswitch

	return 0
End

////////////////////////////////////////////////////////////

Function UserCursorAdjust(graphName,autoAbortSecs)
String graphName
Variable autoAbortSecs
DoWindow/F $graphName // Bring graph to front
if (V_Flag == 0) // Verify that graph exists
Abort "UserCursorAdjust: No such graph."
return -1
endif
NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
DoWindow/C tmp_PauseforCursor // Set to an unlikely name
AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
DrawText 21,20,"Adjust the cursors and then"
DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdjust_ContButtonProc
Variable didAbort= 0
if( autoAbortSecs == 0 )
PauseForUser tmp_PauseforCursor,$graphName
else
SetDrawEnv textyjust= 1
DrawText 162,103,"sec"
SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
Variable td= 10,newTd
Variable t0= ticks
Do
newTd= autoAbortSecs - round((ticks-t0)/60)
if( td != newTd )
td= newTd
SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
if( td <= 10 )
SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
endif
endif
if( td <= 0 )
DoWindow/K tmp_PauseforCursor
didAbort= 1
break
endif
PauseForUser/C tmp_PauseforCursor,$graphName
while(V_flag)
endif
return didAbort
End


/////////////////////////////////////////////////////////

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
String ctrlName
DoWindow/K tmp_PauseforCursor // Kill panel
DoWindow/K Linescan
End


/////////////////////////////////////////////////////////////////////////////////
Function DefineAZsOld(waveName)
	String waveName
	
	Variable AZleftP_1, AZrightP_1, AZleftP_2, AZrightP_2, AZleftP_3, AZrightP_3
	//At this point, the base name has "_L" added to it.  We now remove this and add the suffix "_AZs"
	//The ".._AZs" wave will contain the left and right limits of each AZ 
	String AZname=(waveName[0, (strlen(waveName)-3)])+"_AZs"		
	string out=(waveName[0, (strlen(waveName)-3)])+"_profile"
	
	duplicate/o $waveName $out
	duplicate/o $waveName w1
         
   make/o/n=(dimsize(w1,1)) w2
   matrixop/o w2=sumCols(w1)
   make/o/n=(dimsize(w2,1)) w3
   w3=w2[0][p]
   duplicate/o w3 $out
         
   display/k=1/N=Profile $out
   Legend/C/N=text0/J/F=0/A=RT "\\f01\\Z14Define AZ's using cursor pair(s) and close window"
   ShowInfo/CP={0,1,2}/W=Profile
   MoveWindow/W=Profile 0, 400, 600, 600
   CursorInput(0, "Profile")
         
   AZleftP_1=pcsr(A)
   AZrightP_1=pcsr(B)
   
   string csrResult
   csrResult=CsrInfo(C,"Profile")
   
   if(cmpstr(csrResult, "")==0)
   	AZleftP_2=0
   	AZrightP_2=0
   	AZleftP_3=0
   	AZrightP_3=0
   else
   	AZleftP_2=pcsr(C)
   	AZrightP_2=pcsr(D)
   endif
   
   csrResult=CsrInfo(E,"Profile")
   if(cmpstr(csrResult, "")==0)
   	AZleftP_3=0
   	AZrightP_3=0
   else
   	AZleftP_3=pcsr(E)
   	AZrightP_3=pcsr(F)
   endif
 
        
   Make/O/N=6 $AZname, temp
   temp= {AZleftP_1,AZrightP_1,AZleftP_2,AZrightP_2, AZleftP_3,AZrightP_3}
   Duplicate/O temp $AZname 
   KillWindow Profile 
             
	KillWaves temp, w1, w2, w3, $out	 
    
End          
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////// Define AZs ////////////////////////////////////////
Function azButtOld(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Detect AZs")
			string list=wavelist("*_L",";","DIMS:2")
			string waveName
			prompt waveName, "Line scan to analyse", popup,list
			doprompt "Pick your Line scan", waveName
			DefineAZsOld(waveName)	   	
		break
				
		case -1: // control being killed
			break
	endswitch

	return 0
End





/////////////////////////////////////////////////////////////////////////
Function CursorInput(autoAbortSecs, graphName)
	Variable autoAbortSecs
	String graphName
	
	if (UserCursorAdjust_L(graphName,autoAbortSecs) != 0)
		return -1
	endif

	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
	endif
End
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
Function UserCursorAdjust_ButtonProc(ctrlName) : ButtonControl
	String ctrlName

	DoWindow/K tmp_PauseforCursor				// Kill self
	DoWindow/K LineScan							//Kill Profile window too
End
////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////
Function UserCursorAdjust_L(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs

	DoWindow/F $graphName							// Bring graph to front
	if (V_Flag == 0)									// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName			// Put panel near the graph

	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ButtonProc
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
				
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	
	return didAbort
End
////////////////////////////////////////////////////////////////////





/////////////////////////
Function temporalProfileButt(ba): ButtonControl

Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Temporal Profile")
			string list=wavelist("*_L",";","DIMS:2")
			string waveName
			prompt waveName, "Line scan wave to analyse (ends in _L)", popup,list
			doprompt "Pick your wave", waveName
     		if(V_flag==1)
				Abort
			endif  	
			
			string muName = (waveName[0, (strlen(waveName)-3)])+"_muHat"
	string AName = (waveName[0, (strlen(waveName)-3)])+"_AHat"
	string out=(waveName[0, (strlen(waveName)-3)])+"_profile"
	string coefName = (waveName[0,strlen(waveName)-3])+"_Coef"
	string fitName =(waveName[0,strlen(waveName)-3])+"_fit"
	string spatFilt = (waveName[0,strlen(waveName)-3])+"_SF"
	string sigmaName = (waveName[0,strlen(waveName)-3])+"_sigmaHat"
			
			
			
						
			string list1=wavelist("*_SF",";","DIMS:2")
			string AZName
			prompt AZName, "Corresponding AZ wave (ends in _SF)", popup,list1
			doprompt "Pick your wave", AZName
     		if(V_flag==1)
				Abort
			endif
			
			TemporalProfile(WaveName,aName,muName,sigmaName)
		break
				
		case -1: // control being killed
			break
	endswitch

	return 0
End


//////////////////////////Temporal Profile
Function TemporalProfile(WaveName,aName,muName,sigmaName)
	string WaveName, aName, muName, sigmaName
	duplicate/o $WaveName, lData
	duplicate/o $aName, aHat
	duplicate/o $muName, muHat
	duplicate/o $sigmaName, sigmaHat
	
	variable dt = dimdelta(lData,0)
	variable nPoints = dimsize(lData,0)
	variable nX = dimsize(lData,1)
	make/o/n=(nX) xS
	variable i
	for (i=0;i<nX;i+=1)
		xS[i] = i
	endfor
	variable nROI = dimsize(aHat,0)
	
	
	/// AZs
	make/o/n=(nPoints) tempRoi
	killwaves/z deconMat
	//wave deconMat
	//wave roiDat
	variable j, k
	for (i=0;i<nRoi;i+=1)
	tempRoi =0
		for (j=0;j<nPoints;j+=1)
			for (k=0;k<nX;k+=1)
				tempRoi[j] += lData[j][k] *aHat[i] * (2 * pi * sigmaHat[i])^-.5 * exp(-(xS[k]-muHat[i])^2/(2 * sigmaHat[i]))
			endfor
		endfor
		concatenate/np=1 {tempRoi}, roiDat
	endfor
	
	string azName1 = (waveName[0, (strlen(waveName)-3)])+"_AZ1"
	string decon1 = azName1 + "_D"
	duplicate/o/RMD=[][0] roiDat, $azName1
	setscale/P x,0,dt, $azName1
	BaseCorrection(azName1)
	DeltaF(azName1)
	newWiener(azName1)
	concatenate/np=1 {$decon1}, deconMat
	display/k=1 $decon1 as "Temporally Deconvolved AZ1"
	label left "Deconvolved Amplitude"
	label bottom "Time (s)"
	if (nRoi>=2)
		string azName2 = waveName[0, (strlen(waveName)-3)]+"_AZ2"
		duplicate/o/RMD=[][1] roiDat, $azName2
		setscale/P x,0,dt, $azName2
		BaseCorrection(azName2)
		DeltaF(azName2)
		newWiener(azName2)
		string decon2 = azName2 + "_D"
		concatenate/np=1 {$decon2}, deconMat
		print decon2
			display/k=1 $decon2 as "Temporally Deconvolved AZ2"
		label left "Deconvolved Amplitude"
		label bottom "Time (s)"
	endif
	if (nRoi>=3)
		string azName3 = waveName[0, (strlen(waveName)-3)]+"_AZ3"
		
		duplicate/o/RMD=[][2] roiDat, $azName3
		setscale/P x,0,dt, $azName3
		BaseCorrection(azName3)
		
	DeltaF(azName3)
	newWiener(azName3)
	
	string decon3 = azName3 + "_D"
		concatenate/np=2 {$decon3}, deconMat
		display/k=1 $decon3 as "Temporally Deconvolved AZ3"
	label left "Deconvolved Amplitude"
	label bottom "Time (s)"
	endif
	if (nRoi >=4)
		string azName4 = waveName[0, (strlen(waveName)-3)]+"_AZ4"
		
		duplicate/o/RMD=[][3] roiDat, $azName4
		setscale/P x,0,dt, $azName4
		BaseCorrection(azName4)
	DeltaF(azName4)
	newWiener(azName4)
	string decon4 = azName4 + "_D"
		concatenate/np=2 {$decon4}, deconMat

	display/k=1 $decon4 as "Temporally Deconvolved AZ4"
	label left "Deconvolved Amplitude"
	label bottom "Time (s)"
	endif
	if (nROI >=5)
		string azName5 = waveName[0, (strlen(waveName)-3)]+"_AZ5"
	
		duplicate/o/RMD=[][4] roiDat, $azName5
		setscale/P x,0,dt, $azName5
		BaseCorrection(azName5)
	DeltaF(azName5)
	newWiener(azName5)
	string decon5 = azName5 + "_D"
		concatenate/np=2 {$decon5}, deconMat
		display/k=1 $decon5 as "Temporally Deconvolved AZ5"
	label left "Deconvolved Amplitude"
	label bottom "Time (s)"
	endif
	
	if (nRoi>=6)
		string azName6 = waveName[0, (strlen(waveName)-3)]+"_AZ6"
		duplicate/o/RMD=[][5] roiDat, $azName6
		setscale/P x,0,dt, $azName6
		BaseCorrection(azName6)
	DeltaF(azName6)
	newWiener(azName6)
	string decon6 = azName6 + "_D"
		concatenate/np=2 {$decon6}, deconMat
		display/k=1 $decon6 as "Temporally Deconvolved AZ6"
	label left "Deconvolved Amplitude"
	label bottom "Time (s)"
	endif
	string roiMat = (waveName[0,strlen(waveName)-3]) + "_roiDatMat"
	duplicate/o roiDat, $roiMat
	string baseWaveName = (waveName[0,strlen(wavename)-3])
	killwaves aHat,lData,muHat,sigmaHat,tempRoi,roiDat
	string deconMatName = (waveName[0,strlen(wavename)-3]) + "_deconMat"
	duplicate/o deconMat, $deconMatName
	killwaves/z deconMat
	killwaves/z xS
	
	
	
	end
	

////////////////////////////////////////////////////////////////////////////////////////
//BaseCorrection only does a linear correction, ignoring any initial faster component
////////////////////////////////////////////////////////////////////////////////////////
Function BaseCorrection(waveName)
	String waveName
	
	Duplicate/O $waveName, w
	variable deltat=dimdelta($waveName,0)
	variable nPoints=dimsize($waveName,0)
	duplicate/o w, fit
	
	
	variable endX = pnt2x(w,nPoints)
			
//Choosing here 7 s and 45 s as being roughtly in middle of baselines before and after the stimuus period			
	variable position1=7
	variable position2=45
	Variable m 
	m = (mean($waveName, endX-3, endX-1)- mean($waveName, 5, 9)) / (endX-3-7)
	fit = m*x
	w-=fit
	
	
	
	string out=waveName
	duplicate /o w,$out
	KillWaves/z w, fit
	
End



/////////////////////////// Non Linear Baselines
Function nonLinearBase(wavename)
string wavename
duplicate/o $wavename, w, fit

variable stimStart = 5
variable recLen = pnt2x(w,numpnts(w))

variable mean1 = mean(w,stimStart-4,stiMStart-1)
variable mean2 = mean(w,recLen-4,recLen-1)


variable ratio = mean1/mean2
fit = 0
fit = (ratio-1)*(p) / numpnts(w)+1
w*=fit
string out = waveName
duplicate/o w, $out
killwaves/z w, fit
end






////////////////////////////////////////////////////////////////////////////////////////
	Function/S DeltaF(waveName)
	String waveName

	Variable lengthWave=dimsize($waveName, 0) * dimdelta($waveName, 0)
	Variable x0=0.05*lengthWave, x1=0.95*lengthWave
	string AZname=waveName[(strlen(wavename)-4),strlen(wavename)] 
	string out=waveName[0, (strlen(wavename)-4)] +"DF" + AZname
			
	duplicate/o $waveName tempwave
			
	Make/N=100/O Hist;DelayUpdate
	Histogram/B=1 $waveName, Hist;DelayUpdate
	Duplicate/O Hist, Hist_smth;DelayUpdate
	Smooth 1, Hist_smth;
			
	FindPeak/Q Hist_smth
	tempwave=(tempwave/V_PeakLoc)-1        
	Duplicate/O tempwave $waveName
	KillWaves/z Hist, Hist_smth, tempwave
	return out      
End



////////////////////////////////////////// Temporal Profile ///////////////////////////////////////
Function temporalProfileButtOld(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Temporal Profile")
			string list=wavelist("*_L",";","DIMS:2")
			string waveName
			prompt waveName, "Line scan wave to analyse (ends in _L)", popup,list
			doprompt "Pick your wave", waveName
     		if(V_flag==1)
				Abort
			endif
			
//Get the wave containing the positions of the active zones, made by the DefineAZs function.
// NB As currently configured, you can only hav a maximum of three AZ's
			
			string list1=wavelist("*_AZs",";","DIMS:1")
			string AZName
			prompt AZName, "Corresponding AZ wave (ends in _AZs)", popup,list1
			doprompt "Pick your wave", AZName
     		if(V_flag==1)
				Abort
			endif
	
			variable deltat=dimdelta($waveName,0)
			variable nPoints=dimsize($waveName,0)
			Duplicate/o $AZName tempAZ
			
			
/////////////////////////////////////////
	TemporalProfileOld(waveName)             
//////////////////////////////////////////
		
		case -1: // control being killed
			break
	
	endswitch
	//KillWaves tempAZ, tempAZ1, tempAZ2, tempAZ3

	return 0
	
End


Function TemporalProfileOld(waveName)

	string waveName
	
//Takes the line scan 2D wave, which will end with "_L".  
//Then opens the positions of the active zones, which will end with "_AZs".
//made by the DefineAZs function. Corresponding definitions of active zones will end with "_AZs" but without the "_L".

// NB As curently configured, you can only have a maximum of three AZ's
	
	String AZname=(waveName[0, (strlen(waveName)-3)])+"_AZs"	
	variable deltat=dimdelta($waveName,0)
	variable nPoints=dimsize($waveName,0)
	Duplicate/o $AZName tempAZ

//Get the response in AZ1			
	string outName1=waveName[0, (strlen(waveName)-3)]+"_AZ1"
	Variable leftPoint = tempAZ[0]
	Variable rightPoint = tempAZ[1]
	Variable nPointsAZ1 = tempAZ[1]-tempAZ[0]
	killwaves/z decMat
	Make/O/N=(nPoints) $outName1
	Duplicate/O/R=[][leftPoint,rightPoint]  $waveName tempAZ1
	matrixop/o $outName1=sumRows(tempAZ1)
	duplicate/O $outName1, tempAZ1
	tempAZ1/=nPointsAZ1
	Redimension/N=(nPoints) $outName1, tempAZ1
	Duplicate/O tempAZ1, $outName1
	setscale/P x, 0, deltat, $outName1, tempAZ1 
   BaseCorrection(outName1)
   
   DeltaF(outName1)
   //nonLinearBase(outName1)
   //EventD(outName1)
   newWiener(outName1)
   //setscale/P x, 0, deltat, $outName1 
   duplicate/o $outName1, tempAZ1
   string deconMat1 = outName1 + "_D"
   duplicate/o $deconMat1, dec1
   concatenate/NP=1 {dec1}, decMat

//Get the response in AZ2, if it exists			
	if(tempAZ[2]!=0)
		string outName2=waveName[0, (strlen(waveName)-3)]+"_AZ2"
		leftPoint = tempAZ[2]
		rightPoint = tempAZ[3]	
		Variable nPointsAZ2 = tempAZ[3]-tempAZ[2]
		
		Make/O/N=(nPoints) $outName2
		Duplicate/O/R=[][leftPoint,rightPoint]  $waveName tempAZ2
		matrixop/o $outName2=sumRows(tempAZ2)
		duplicate/O $outName2, tempAZ2
		tempAZ2/=nPointsAZ2
		Duplicate/O tempAZ2, $outName2
		Redimension/N=(nPoints) $outName2, tempAZ2
		setscale/P x, 0, deltat, $outName2 
      BaseCorrection(outName2)
   	DeltaF(outName2)
   	newWiener(outName2)
   
   	//EventD(outName2)
   	//setscale/P x, 0, deltat, $outName1 
   	duplicate/o $outName2, tempAZ2
   	 duplicate/o $outName1, tempAZ1
   string deconMat2 = outName2 + "_D"
   duplicate/o $deconMat2, dec2
   concatenate/NP=1 {dec2}, decMat


	endif
 
 //Get the response in AZ3, if it exists        
    if(tempAZ[4]!=0)
		string outName3=waveName[0, (strlen(waveName)-3)]+"_AZ3"
		leftPoint = tempAZ[4]
		rightPoint = tempAZ[5]
		Variable nPointsAZ3 = tempAZ[5]-tempAZ[4]
					
		Make/O/N=(nPoints) $outName3
		Duplicate/O/R=[][leftPoint,rightPoint]  $waveName tempAZ3
		matrixop/o $outName3=sumRows(tempAZ3)
		duplicate/O $outName3, tempAZ3
		tempAZ3/=nPointsAZ3
		Redimension/N=(nPoints) $outName3
		setscale/P x, 0, deltat, $outName3
		BaseCorrection(outName3)
   	DeltaF(outName3)
   	newWiener(outName3)
   	//EventD(outName3) 
   	//setscale/P x, 0, deltat, $outName1
   	duplicate/o $outName3, tempAZ3
   	 duplicate/o $outName1, tempAZ1
   string deconMat3 = outName3 + "_D"
   duplicate/o $deconMat3, dec3
   concatenate/NP=1 {dec3}, decMat

 
	endif                   

	KillWaves tempAZ, tempAZ1, tempAZ2, tempAZ3
	Execute "TileWindows/O=1/C"
	string deconMatName = (waveName[0,strlen(wavename)-3]) + "_deconMat"
	duplicate/o decMat, $deconMatName
	
	
	
	return 0

End



Function eventDetBUTT(ba): ButtonControl
Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Temporal Profile")
			string list=wavelist("*deconMat",";","DIMS:2")
			string waveName
			prompt waveName, "Select ROI Data Matrix Wave", popup,list
			doprompt "Pick your wave", waveName
     		if(V_flag==1)
				Abort
			endif
			duplicate/o $waveName, deconMat
			variable nROIs = dimsize(deconMat,1)
			variable nPoints = dimsize(deconMat,0)
			string threshName = (wavename[0,strlen(wavename)-10]) + "_thresh"
			duplicate/o $threshName, thr
			string baseName = (wavename[0,strlen(wavename)-10]) + "_AZ"
			string e1 = basename + "1_E"
			string a1 = basename + "1_A"
			//wave ampMat,evMat
			variable roi,i
			variable nEvent = 0
			print nRois
			for (roi=0;roi<nROIs;roi+=1)
				duplicate/o/RMD=[][roi] deconMat, decWave, tempA,tempE
				tempE = 0
					tempA = 0
					nEvent = 0
				for (i=1;i<nPoints-1;i+=1)
					
					if (decWave[i] > decWave[i-1] && decWave[i]>decWave[i+1] && decWave[i] > thr[roi])
						tempA[nEvent] = decWave[i]
						tempE[nEvent] = pnt2x(decWave,i)
						nEvent+=1
					endif
				endfor
				if (roi==0)
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a1
					duplicate/o tempE, $e1
				elseif (roi==1)
					string e2 = basename + "2_E"
					string a2 = basename + "2_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a2
					duplicate/o tempE, $e2
				elseif (roi==2)
					string e3 = basename + "3_E"
					string a3 = basename + "3_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a3
					duplicate/o tempE, $e3
				elseif (roi==3)
					string e4 = basename + "4_E"
					string a4 = basename + "4_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a4
					duplicate/o tempE, $e4
				elseif (roi ==4)
					string e5 = basename + "5_E"
					string a5 = basename + "5_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a5
					duplicate/o tempE, $e5
				elseif (roi==5)
					string e6 = basename + "6_E"
					string a6 = basename + "6_A"
					deletepoints nEvent,nPoints, tempA
					deletepoints nEvent,nPoints, tempE
					duplicate/o tempA, $a6
					duplicate/o tempE, $e6
				endif
			endfor	
		endswitch
		
		killwaves/z decWave,deconMat,data,aHat, decWave, tempE,tempA
end




Function clusterButt(ba): ButtonControl
Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Temporal Profile")
			string list=wavelist("*_A",";","DIMS:2")
			string waveName
			prompt waveName, "Select Amplitude Wave", popup,list
			doprompt "Pick your wave", waveName
     		if(V_flag==1)
				Abort
			endif
			emGMM(waveName)
		endswitch
	end
	





Function threshButt(ba): ButtonControl
Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2: 
			print ("You booped Temporal Profile")
			string list=wavelist("*deconMat",";","DIMS:2")
			string waveName
			prompt waveName, "Select ROI Data Matrix Wave", popup,list
			doprompt "Pick your wave", waveName
     		if(V_flag==1)
				Abort
			endif
			
			threshGui(waveName)
			
			

				
			
		endswitch
end



/////////////////////////////////////////////////////////////////
///////////////////////////IEI /////////////////////////
///////////////////////////////////////////////////////////
function IEI(ba): ButtonControl
Struct WMButtonAction &ba
	switch(ba.eventCode)
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
				
				
			VARIABLE stimStart = 10
			variable stimStop = 40
			prompt stimStart, "Enter stim Start"
			doprompt "Enter stimStart", stimStart
			
			
			if(V_flag==1)
			abort
			endif
			prompt stimStop, "Enter stim Stop"
			doprompt "Enter stimstop", stimStop
			if(V_flag==1)
			abort
			endif
			ieiAnal(eventname, ampname,stimStart, stimStop)
		endswitch
end


///////////////////////////////////////////////////////
////////////////////// IEI Anal ////////////////////
/////////////////////////////////////////////////////
function ieiAnal(eventname, ampname,stimStart, stimStop)

	string eventName, ampName
	variable stiMStart, StimStop
	duplicate/o $eventname,events
	duplicate/o $ampName, amps
	
	
	variable nEv = dimsize(events,0)
	make/o/n=(15,nEv) ieiN = NaN
	MAKE/O/N=15 NNQ=0
	variable i,j
	
	for (i=1;i<16;i+=1)
		duplicate/o events, nnEvents

		VARIABLE ONwHICH = 0
		for (j=0;j<nEv;j+=1)
			if (AMPS[J]==I && events[j] > stimStart && events[j] < stimStop)
				NNeVENTS[ONwHICH]= EVENTS[J]
				ONwHICH+=1
			ENDIF
		endfor
		variable nnEv = dimsize(nnEvents,0)
		for (j=1;j<ONwHICH;j+=1)
			ieiN[i-1][j-1] = nnEvents[j]-nnEvents[j-1]
		endfor
	endfor
	KILLWAVES/Z NNeVENTS, AMPS, EVENTS, NNQ
	
	string ieiName = eventName + "_IEI"
	duplicate/o ieiN, $ieiName
	killwaves/z ieiN
end
	
			