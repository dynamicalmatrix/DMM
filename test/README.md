How to test:
1. Remove/rename the previous my-output.dmm if exists;
2. start OOMMF;
3. load (in Oxsii) test.mif;
4. send the Magnetization to mmDisp;
5. run;
6. in mmDisp save the magnetization as test.omf (Text format);
7. run the main program (e.g. ../micro3d inputfile.template)
8. find the output in my-output.dmm.

Then the output can be converted to mode profiles with the utility program "dmmtosilo" (directory "uti").
