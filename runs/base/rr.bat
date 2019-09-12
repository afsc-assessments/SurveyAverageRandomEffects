@echo off
rw -ind sims\sim_%1.dat -nox >NULL
cat rw.rep >>simsetout.rep
copy rwout.rep arc\simout_%1.rep
echo %1
