********************************
*  SIDIS Cross Section Model   *
*  -- Zhihong Ye, 06/10/2014   *
********************************

Note:
1, This model is extracted fro Xin Qian's SIDIS generator "collider.C" 
   which is included as a reference.
2, The XS model is coded in "SIDIS.h". 
  1) A generator is given in "GetSIDIS.C" (a bug is found, please use collider.C as the generator for now-Z. Ye,07/17/2014).
  2) A simpler example of using the model is given in "example.C"
3, LHAPDF is needed to be installed. Unpar "lhapdf.5.8.8.tar.gz", 
   and follow the instruction to install it. Specify the path in "SIDIS.h"
4, A older version of the model wrotten in FORTRAN is also given in "sidis_fortran/", as a reference.
5, I held no responsiability to this code but you are welcome to discuss questions with me :->
