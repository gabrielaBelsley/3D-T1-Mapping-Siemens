# 3D-T1-Mapping-Siemens

## 3D T1 MAPPING using VFA SPGR 
**SIEMENS Data acquired on a 3T Prisma scanner**

Scripts to compute whole-liver T1 maps using the widely available SPGR sequence with corrections for B1+ inhomogeneities.
B1+ map generated from a GRE EPI sequence with fat saturation using a dual-angle approach and modelling:
1. Slice profile effects
2. B0 gradient through slice de-phasing
3. Detection and elimination of Fat band from EPI chemical shift 
4. Correction for EPI distortion of GRE EPI images

An alternative B1+ mapping method-pre-conditioned RF pulse with a Turbo FLASH readout- which is available in some Siemens 3T scanners (e.g. Prisma) is also used to correct the nominal FAs in the VFA acquisition. 

T1 Maps corrected for incomplete spoiling and with spatial saturation turned off.

This code is distributed under the MIT license. If you use it, please cite the code: TO DO: include zenodo doi when I make the code public

Author: Gabriela Belsley, University of Oxford, gabriela.belsley@stx.ox.ac.uk OR gabi.belsley@gmail.com

Please contact me, if you have any questions. Happy T1 Mapping!
————————————————————————————————————————————

**General Instructions:**

Run MainScript.m: it has several sections to run. 

It can be modified to run everything at once and more efficiently, if one wishes to for batch processing. 

The data needs to be extracted and saved in the following directory tree:

Main Folder: Prisma/Run1/

Subfolders

1. 2DGRE_EPI/FS/FA65 + 2DGRE_EPI/FS/FA130: place GRE EPI
2. 3DVIBE_VFA_Dixon/FA2 + 3DVIBE_VFA_Dixon/FA2R + 3DVIBE_VFA_Dixon/FA15 +3DVIBE_VFA_Dixon/FA15R: place water only images, PrescanOn is the second set (first set is PrescanOff). Use Prescan On for T1 Mapping.
3. Siemens2DMultiSlice_B1Map/
    - FA80: place Siemens B1+ map images computed directly at the scanner
    - Magnitude/FA80: place Siemens Pre-RF turbo FLASH magnitude only images 
4. B0_2DGREMultiSlice/UnwrapPrelude_B0Fugue/
    - Phase/TE1/FA15 + Phase/TE2/FA15 
    - Magnitude/TE1/FA15 + Magnitude/TE2/FA15 
    - PhaseTE1  + PhaseTE2/
    - MagTE1/
    - greEPI/FS/FA65 + greEPI/FS/FA130 + greEPI/FS/NoRegist


Dependencies:
Uses the Bloch simulator written by Professor Brian Hargreaves: download it from http://mrsrl.stanford.edu/~brian/blochsim/
