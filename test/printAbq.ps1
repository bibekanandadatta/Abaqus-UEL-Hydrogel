# to run PowerShell scipt, type: .\script_name on PowerShell terminal

# Get-Location == pwd (prints out working directory)
# Get-Location

# Deckare variables with $ in front of it
$JOBNAME    = "U4_hydrogel_1x1_cylinder_elem_64"
$STAFILE   	= "$JOBNAME.sta"

clear

echo "PRINING: $STAFILE"

# Set-Location == cd (change directory)
# Set-Location -path $WORKDIR

# wait until the status file is created
while (!(Test-Path $STAFILE)) 
{
    # check every 10 seconds
    Start-Sleep 5
}

# print out the Abaqus job status file
Get-Content -wait $STAFILE