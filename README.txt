  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHORS:
        
        Taewon Cho, Department of Mathematics, Virginia Tech

        Julianne Chung, Department of Mathematics, Virginia Tech
       
        Hodjat Pendar, Department of Biomedical Engineering and Mechanics, Virginia Tech
   
   REFERENCE:

       "Computational Tools for Inversion and Uncertainty Estimation in Respirometry". 
            2021.

   SOFTWARE LANGUAGE:

       MATLAB 9.6 (R2019a)


=====================================================================
SOFTWARE
=====================================================================
The MainDrivers require the following packages:

   (1) RestoreTools by James Nagy, Katrina Palmer, and Lisa Perrone
	http://www.mathcs.emory.edu/~nagy/RestoreTools/index.html

   (2) IRTools by James Nagy, Sivia Gazzola, and Per Christian Hansen
	https://github.com/silviagazzola/IRtools

   (3) genHyBR by by Julianne Chung and Arvind Saibaba
	https://github.com/juliannechung/genHyBR

   (4) Codes from the SIAM Book by Johnathan M. Bardsley
      "Computational Uncertainty Quantification for Inverse Problems"
        https://github.com/bardsleyj/SIAMBookCodes  

   (5) FISTA, a MATLAB implementation by Tiep Vu
	https://github.com/tiepvupsu/FISTA

First run startup_Respirometry.m for setting paths.

MainDrivers for each numerical experiments of inverse problems

  EX_1_MainDriver_sim.m           Sets up and runs a linear inverse
                                  problem corresponding to Sections  
                                  "Linear Respirometry Reconstruction"
                                  and "Nonlinear Respirometry
                                  Reconstruction" with simulations.
                                  Generates Figures 1-10

  EX_2_MainDriver_real.m          Sets up and runs real blind
                                  deconvolution corresponding to
                                  "Nonlinear case study: Abdominal
                                  pumping and CO2 emission in a darkling
                                  beetle" and generating Figures 13, 14

  Ex_3_MainDriver_real_true.m     Sets up and runs real linear case study
                                  generating Figures 11, 12
