del simsetout.rep
FOR /L %%y IN (1,1,3600) DO call r.bat %1 %%y   
copy simsetout.rep %1.rep
