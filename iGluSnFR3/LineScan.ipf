#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////////////////////
Function LSimage(w)
	
	// this code will take a linescan file from scanimage and replot it so that you can see the whole record with time on the 
	// x axis and position on the y. 
	
	wave w		// input of the line scan image
	variable i, lx=dimsize(w,0),ly=dimsize(w,1),lz=dimsize(w,2)  // counter and length of x, y and z
	
	duplicate/o/FREE w, temp		// duplicate it
	
	redimension/N=(lx*ly*lz) temp	// these 2 redimension commands will transform it to look pretty
	redimension/N=(lx,ly*lz) temp
	
	matrixtranspose temp 		// flip the data so that time is on the x axis
	
	string notes=note(w)		// get the info from the header to scale the image properly
	variable msPline = str2num(notes[(strsearch(notes,"state.acq.msPerLine=",0)+20),(strsearch(notes,"state.acq.fillFraction=",0)-2) ])
	variable zoom = str2num(notes[(strsearch(notes,"state.acq.zoomFactor=",0)+21),(strsearch(notes,"state.acq.scanAngleMultiplierFast=",0)-2) ])
	
	setscale/I y,0,(610/zoom), temp			// scale the data
	setscale/P x,0,(msPline/1000), temp
	
	
	
	string name=nameofwave(w)		//new name for wave
	string outName=(name[0, (strlen(name)-3)])+"_L"
	duplicate/o temp, $outName
	
	// code below is for making a plot of the data
	
	display/K=1/N=Linescan
	appendimage $outName
	Label bottom "Time (s)"
	Label left "Position (µm)"
	imagestats temp
	ModifyImage $outName ctab= {*,(V_max*0.7),YellowHot256,0}
	ModifyGraph mirror(left)=0,standoff(left)=0
	ModifyGraph mirror(bottom)=0,standoff(bottom)=0	
	
end

//////////////////////////////////////////////////////////////////