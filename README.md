# toolpackage5_v2
This is a ROOT6-based program for TOF or Alpha/Beta-TOF analysis.
To parse the MRTOF.lst file or do the drift correction, launch the root --> .L MassAnaGUI4.cc+ ( '+' at the end means load and compile. If does not work, quit and launch ROOT again, then use ".L MassAnaGUI4.cc" to load only without compile, although compile is always recommend.)  --> execute MassAnaGUI4() to open the control panel

Main program for mass and Alpha/Beta-TOF analysis is "preview5.C". Using above way to load and compile it. Then use the preview5() or autopreview() function to open spectrum. Documents in "demo" folder introduce in detail.

If just want to have a look the spectrum of Alpha/Beta.lst files or just want to convert data from .lst to a tree. Beta_Beta_coin.cc can be used independently.
".L Beta_Beta_coin.cc+" -> Read_Beta_lst_batch("path of .lst","name of .lst file with or without .lst");   set "MakeBetaRawTree=false = true" before executing Read_Beta_lst_batch() if wanna to convert data from .lst to a tree. --> [Option: execute Compare_Beta_Beta() to get beta-beta coincident candicates.] --> ShowBeta_Beta_Coin() to view histogram;


Main improvement: 
TOF:
1. Decoder and drift correction control window can combine .lst; select specific channel for decode, and/or select specific channel and veto time window for drift correction. Sort data in sequence.
2. preview5.C can select data from specific MCS channel and set veto time window for analysis

Beta/Alpha
1. Compatible to .txt and big-endian binary beta.lst.
2. Switch between different data read modes to handle large beta data set.
3. Switch between alpha and beta analysis

For detail introduction, please refer to documents in "demo" and videos below.

Demo video：
https://1drv.ms/f/c/107fa6fc1cbfed8d/EkB1vSMLjp9IujlS4tcKGQMBETyCDRof4BUFXrWUEbnRgQ?e=cy6wk1
