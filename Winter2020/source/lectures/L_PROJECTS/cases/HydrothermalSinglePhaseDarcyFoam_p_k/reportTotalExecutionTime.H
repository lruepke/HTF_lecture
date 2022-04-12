// report progress
const scalar progress = runTime.time().value()/endTime;
Info<< "Time = " << runTime.timeName() <<" s = "<<runTime.time().value()/sec_year<<" year "
<<" Progress = "
<<progress*100.0<<"%"
<<nl<<endl;

// Report total execution time
const scalar totalTime_clock = runTime.elapsedCpuTime();
const scalar totalTime_cpu = runTime.elapsedClockTime();
const scalar totalTimeEstimate = totalTime_clock/progress;
if(totalTime_clock<=60)
{
    Info<<"ExecutionTime = "<<totalTime_clock<<" s"
        <<" ClockTime = "<<totalTime_cpu<<" s"
        << nl << endl;
}else if(totalTime_clock>60 && totalTime_clock<=3600)
{
    Info<<"ExecutionTime = "<<totalTime_clock/60<<" m"
        <<" ClockTime = "<<totalTime_cpu/60<<" m"
        << nl << endl;
}else if(totalTime_clock>3600 && totalTime_clock<=86400)
{
    Info<<"ExecutionTime = "<<totalTime_clock/3600<<" h"
        <<" ClockTime = "<<totalTime_cpu/3600<<" h"
        << nl << endl;
}else if(totalTime_clock>86400 && totalTime_clock<=2592000)
{
    Info<<"ExecutionTime = "<<totalTime_clock/86400<<" d"
        <<" ClockTime = "<<totalTime_cpu/86400<<" d"
        << nl << endl;
}else
{
    Info<<"ExecutionTime = "<<totalTime_clock/2592000<<" months"
        <<" ClockTime = "<<totalTime_cpu/2592000<<" months"
        << nl << endl;
}
