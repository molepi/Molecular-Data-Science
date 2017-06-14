:: windows software install script for FOS course 2017
:: Maarten van Iterson
:: Tue Jun 13 13:49:08 2017
:: useful resource: http://steve-jansen.github.io/guides/windows-batch-scripting/

:: set variables
SET SRC=
SET DST=C:\FOS2017

:: create folder
rmdir /S %DST%
mkdir %DST%


if not exist Cats goto driveE
XCOPY /E/Y . %DST%
goto copied

:driveE
if not exist E:Cats goto driveF
XCOPY /E/Y E: %DST%
goto copied

:driveF
if not exist F:Cats goto driveG
XCOPY /E/Y F: %DST%
goto copied

:driveG
pause 
exit

:copied

cd %DST%
C:


%DST%\.exe


C:\Progra~1\R\R-3.4.0\bin\r --vanilla < %DST%\packages2017.r

pause alles geinstalleerd!

exit


