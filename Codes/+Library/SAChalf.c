/*********************************************************************
* SAChalf.c 
* Created by: SP 21 Jun, 2017
* 
* This is a MEX-file for MATLAB.
*
* Computes Shuffled Auto Correlogram for positive delays from a Spike Matrix (symmetric around 0 delay)
* 
* Call:  [SAC,delay] = InnerSACmex(SpikeMAT1,NUMspikes,TOTALspikes);
* 		SpikesMAT1: (Kmax x NUMreps) matrix with spike trains (NaNs fill end of each column for shorter spike trains)
*		NUMspikes: (1x NUMreps) vector with number of spikes in each train
*		TOTALspikes: needed to allocate SAC
*
*		SAC: Shuffled Auto Correlogram (for only positive delays)
*		delay: delay corresponding to SAC. 
*		
*********************************************************************/

#include "mex.h"

//InnerShufAutoCorrMEX(SpikeMAT1,NUMspikes,NUMspikeREPS,Kmax,ints,delay);
void SAChalf(double *SpikeMAT1, double *NUMspikes, long NUMspikeREPS, long Kmax, double *DelayBinWidth, long nzmax, double *ints)
{
   long REPindREF,SpikeIND,REPindCOMP,COMPSpikeIND,intIND=0;
   
   for (REPindREF=0; REPindREF<NUMspikeREPS; REPindREF++) {
	   for (SpikeIND=0; SpikeIND<NUMspikes[REPindREF]; SpikeIND++) {
		   for (REPindCOMP=0; REPindCOMP<NUMspikeREPS; REPindCOMP++) {
			   if (REPindCOMP!=REPindREF) {
				   for (COMPSpikeIND=0; COMPSpikeIND<NUMspikes[REPindCOMP]; COMPSpikeIND++) {
					   intIND=floor((SpikeMAT1[REPindCOMP*Kmax+COMPSpikeIND]-SpikeMAT1[REPindREF*Kmax+SpikeIND])/DelayBinWidth[0]);
					   if (intIND>=0 && intIND<=nzmax){
					   ints[intIND]+=1;
						}
      }  }  }	}  }
}


// the gateway function
// [SAC,delay] = SAChalf(SpikeMAT1,NUMspikes,TOTALspikes);
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   double   *SpikeMAT1,*SAC,*NUMspikes, *DelayMax, *DelayBinWidth; // input-output arrays
	long     TOTALspikes; // input scalar
   double   *delay;  // output scalar
   long     NUMspikeREPS,Kmax, nzmax;

	int      i;

//	__asm int 3h ; // Setting a breakpoint for debugging (compile with the -g flag!)
							// For "Nel_matlab", run 'doallmex' with '(0,1)' flags.
							// Run Nel software as usual, and Matlab will automatically
							//  break at this line and open this file in Visual Studio
							//  debugger environment.

   /*  check for proper number of arguments */
   if(nrhs!=5)
    mexErrMsgTxt("Five inputs required.");
   if(nlhs!=2)     mexErrMsgTxt("One output required.");
  
  	/*  create a pointer to the input matrix SpikeMAT1 */
   // SpikeMAT1 is stored in MATLAB as: each row is one spike train
	//   in MEX/C, linear indexing goes down columns first, so we pass the transpose to 
	//   allow easier indexing
   SpikeMAT1 = mxGetPr(prhs[0]);  // each column is one spike train
	/*  get the dimensions of the matrix input SpikeMAT1 */
   NUMspikeREPS = mxGetN(prhs[0]);
   Kmax = mxGetM(prhs[0]);
  
	/*  create a pointer to the input vector NUMspikes */
   NUMspikes = mxGetPr(prhs[1]);

   /*  get the scalar input TOTALspikes */
   TOTALspikes = mxGetScalar(prhs[2]);

   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
      mexErrMsgTxt("Input TOTALspikes must be a scalar.");
   }
   
   /*  get the scalar input DelayMax */
   DelayMax = mxGetPr(prhs[3]);
   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
      mexErrMsgTxt("Input DelayMax must be a scalar.");
   }

    /*  get the scalar input DelayBinWidth  */
   DelayBinWidth = mxGetPr(prhs[4]);
   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
      mexErrMsgTxt("Input DelayBinWidth must be a scalar.");
   }
   
   
   /*  set the output pointer to the output matrix SAC */
    nzmax = 1+floor(DelayMax[0]/DelayBinWidth[0]);
	//printf("%d is the size \n", nzmax);
	plhs[0] = mxCreateDoubleMatrix(1,nzmax, mxREAL);
   /*  create a C pointer to a copy of the output matrix SAC */
   SAC = mxGetPr(plhs[0]);
  
   /*  set the output pointer to the output scalar delay */
   plhs[1] = mxCreateDoubleMatrix(1,nzmax, mxREAL);
   /*  create a C pointer to a copy of the output variable TOALints */
   delay = mxGetPr(plhs[1]);
    
	for (i=0; i<=nzmax; i++)
   {
	   delay[i]=DelayBinWidth[0]*i;
   }

   /*  call the C subroutine */
   SAChalf(SpikeMAT1,NUMspikes,NUMspikeREPS,Kmax,DelayBinWidth, nzmax, SAC);
	
	}
