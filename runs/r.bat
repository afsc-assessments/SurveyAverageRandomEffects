@echo off
%mod% -ind %1\sim_%2.dat -nox >NUL
cat %mod%.rep >>simsetout.rep
:: copy rwout.rep arc\simout_%1_%2.rep
echo %1 %2 %mod%
