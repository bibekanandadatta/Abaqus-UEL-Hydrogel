# to run PowerShell scipt, type: .\script_name on PowerShell terminal

$INPUTFILE  = "U4_hydrogel_1x1_cylinder_elem_64" 
$JOBNAME   	= $INPUTFILE
$SRC        = "../src/uel_hydrogel.for"
$NPROC      = 1

clear

echo "ABAQUS JOB RUNNING: $JOBNAME"
echo "UEL SUBROUTINE: $SRC"

abaqus interactive double analysis ask_delete=off job=$JOBNAME input=$INPUTFILE user=$SRC cpus=$NPROC

echo "ABAQUS JOB $JOBNAME IS COMPLETED"