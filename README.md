Interface for VASP to Cc4s
==========================


This repository contains a bash script that patches the VASP source to support the Cc4s code and adds a number of additional features.

How to patch your VASP version?
-----------------------------------------

Find the correct patch in this repository depending on your VASP version (e.g. `./patch/vasp.6.3.1/patch_vasp.6.3.1_cc4s.sh`)
Copy the patch into the `./src` directory of your VASP copy.
Run the patch from the `./src` directory of your VASP version by typing, e.g. `./patch_vasp.6.3.1_cc4s.sh`.

If everything works, you will see something like the following in your terminal
```
vasp.6.3.1/src$ ./patch_vasp.6.3.1_cc4s.sh 
Trying to patch VASP version 6.3.1 to include Cc4s interface. 
Patching file  chi.F 
patching file chi.F
Patching file  chi_glb.F 
patching file chi_glb.F
Patching file  main.F 
patching file main.F
Patching file  ump2.F 
patching file ump2.F
Patching file  ump2no.F 
patching file ump2no.F
Patching file  .objects 
patching file .objects
Adding file  bracketst.F 
Adding file  cc4s.F 
Adding file  ccsd.F 
```

Recompile VASP and run the conventional complex or gamma-only vasp binaries.
