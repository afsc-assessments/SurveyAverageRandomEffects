@echo off
rw -ind data\%1.dat -nox 
copy rwout.rep arc\%1.rep
echo %1
