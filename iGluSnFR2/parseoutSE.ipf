#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function slowevents(w,spikes_test,thr)

wave w, spikes_test
variable thr
variable nl=dimsize(w,0),i,k,s=sum(spikes_test),j,r,v, maxi,l, delta=dimdelta(w, 0)
make/o/n=(100,s) event_before
k=0
r=0
 
 //count the number of peaks
maxi=0
	for(i=50;i<nl-60;i+=1)
	  if(spikes_test(i)==1)
	     duplicate/o/free/r=[i-10,i+60] w temp 
	  
	     duplicate/o/free temp tempdif
	  
	    tempdif-=tempdif[0]
	     wavestats/q tempdif 
	    l=area(tempdif)
	   event_before[][r]=temp[p]
	   r+=1
	   

	            if(l>thr)
	            
	       
	       	if (maxi != V_max)
	       	 wavestats/q temp
	               maxi=V_max
	               k+=1
	       	endif
	       endif
         endif
         l=0
         endfor
         //print(k)
         //create a new wave with peaks
          //make/o/n=(31,k) sevents 
          make/o/n=(71,k) sevents 
          make/o/n=(k,1) sevents_int 
          make/o/n=(nl,0) finallocationevent=0
          make/o/n=(nl,0) test=0
    	    make/o/n=(k,1) finalamplitude=0
    	    make/o/n=(k,1) locationwave
    	   // finallocationevent[p][q]=0
          v=0
          
          maxi=0
	   for(j=10;j<nl-60;j+=1)
	      if(spikes_test(j)==1)
	       	
	    //  duplicate/o/free/r=[j-5,j+25] w temp 
	    duplicate/o/free/r=[j-10,j+60] w temp 
	      duplicate/o/free temp tempdif
	  
	        tempdif-=tempdif[0]
	        wavestats/q tempdif 
	        l=area(tempdif)
	     
	     
		if(l>thr)
		          
		  	 if (maxi!=V_max)
		  	 wavestats/q temp
		    	  maxi=V_max
		    	 
		    	  sevents_int[v][0]=l
		    	  duplicate/o/r=[V_maxrowloc-10,V_max] temp baseline
		    	  wavestats/q baseline
		    	
		    	    locationwave[v][1]=j
		    	    //finalamplitude[v][]= maxi-V_min
		    	   // string noncor=nameofwave(w)-"_cor"
		    	    //wave fivehz030_sum
		    	  //  wave w=fivehz030_sum
	           //    duplicate/o/R=[x2pnt(w,locationwave[v])-5,x2pnt(w,locationwave[v])+25 ] w, tempamp
	           //duplicate/o/R=[locationwave[v]-5,locationwave[v]+25 ] w, tempamp
	           duplicate/o/R=[locationwave[v]-10,locationwave[v]+50 ] w, tempamp
	                //print(tempamp)	// the "-5" and "+15" will need changed for different sample rates

		//wavestats/q/R=[39000,39500] w
		//wavestats/q/R=[0,3] tempamp
		wavestats/q/R=[0,6] tempamp
		 tempamp-=V_min
		//tempamp/=V_min						// converts waves to ÆF/F
		duplicate/FREE/o tempamp, tempampdif
		Differentiate tempamp/D=tempampdif			// differentiate to find peak by the location of inflection
		//wavestats/q/R=[6,30] tempamp
		wavestats/q/R=[10,60] tempamp
	//wavestats/q/R=[11,20] tempamp
			finalamplitude[v][]= tempamp[V_maxrowloc]
	       	test[j+V_maxRowLoc-5]= tempamp[V_maxRowLoc]
	       	    finallocationevent[j+V_maxRowLoc-5]=1
	       	    
	       	 
	              sevents[][v]=tempamp[p]
	               v+=1
	               
	               endif
	           endif
	        endif
	    endfor
	
	wave ws
	wave eventwave
	wave ones
	//print(v)
	//print(k)
	setscale/P  x,0, delta, finallocationevent
	locationwave*=delta
	
	 	display/K=1/W=(35,276,772,481) ws					//display data and detected events
	appendtograph/R ones vs eventwave
	ModifyGraph mode(ones)=3,rgb(ones)=(0,0,0)
	SetAxis right 0,3
	ModifyGraph noLabel(right)=2,axThick(right)=0
	
	appendtograph/l  finalamplitude vs locationwave
	ModifyGraph rgb(finalamplitude)=(16385,49025,65535)
	ModifyGraph mode(finalamplitude)=1
	ModifyGraph lsize(finalamplitude)=2
	k=0

//print "threshold was = "+num2str(V_levelX)
	print "found "+num2str(sum(finallocationevent))+" events"
	
end