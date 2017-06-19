Function UserCursorAdjust(nextlooptest,autoAbortSecs)
//String graphName
Variable autoAbortSecs
string nextlooptest
DoWindow/f $nextlooptest            // Bring graph to front

if (V_Flag == 0)                       // Verify that graph exists
Abort "UserCursorAdjust: No such graph."
return -1
endif
NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
DoWindow/C tmp_PauseforCursor          // Set to an unlikely name
AutoPositionWindow/E/M=0/R=$nextlooptest   // Put panel near the graph
DrawText 21,20,"Adjust the cursors and then"
DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdjust_ContButtonProc
Variable didAbort= 0
if( autoAbortSecs == 0 )
PauseForUser tmp_PauseforCursor,$nextlooptest 
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
PauseForUser/C tmp_PauseforCursor,$nextlooptest 
while(V_flag)
endif
return didAbort
End
Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
String ctrlName
DoWindow/K tmp_PauseforCursor             // Kill panel
End
Function Demo(autoAbortSecs)
Variable autoAbortSecs
string name
//Newimage/N=nextlooptest fivehz_057_Ch1_LSim_sy1
DoWindow nextlooptest
ShowInfo
//print(v_flag)
if( V_Flag==0 )
Newimage/N=nextlooptest $ name
ModifyGraph height={Aspect,0.2}, expand=1
ShowInfo
endif

if (UserCursorAdjust("nextlooptest",autoAbortSecs) != 0)
return -1
endif
//place your code
End
