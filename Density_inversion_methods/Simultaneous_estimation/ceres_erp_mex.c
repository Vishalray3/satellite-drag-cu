// Box-wing Earth radiation pressure model using CERES data
// Doornbos et al. air density models derived from multiple satellite observations
// Written by Vishal Ray, Oct 26, 2022
#include "mex.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define m_max  180
#define n_max  360
#define k_max  20
#define pi  3.14159265358979323846
#define delp 100.0
////////////////////////////
//////// Outputs /////////
/////////////////////////

double *Accel_ptr, *Jacobian_ptr;

////////////////////////////
//////// Inputs /////////
/////////////////////////

double *sunvector_ptr,  *nrm_eci_ptr, *x_state_ptr;
double *plt_area_ptr, *refl_vis_spec_ptr, *refl_vis_diff_ptr, *refl_ir_spec_ptr, *refl_ir_diff_ptr;
double *ECEFtoJ2000_ptr;
double *Psw_ptr, *Plw_ptr, *lam_ptr, *phi_ptr;
double scmass,  Cr_erp, Re;
double nrm_eci[3][k_max];
double plt_areas[k_max], refl_vis_spec[k_max], refl_vis_diff[k_max], refl_ir_spec[k_max], refl_ir_diff[k_max];
double ECEFtoJ2000[3][3];
double Psw[n_max+1][m_max+1], Plw[n_max+1][m_max+1], lam[n_max], phi[m_max];
int Mphi, Nlam, Kplt;

////////////////////////////////
/////Intermediate variables/////
////////////////////////////////
int  mm, nn, kk, pp, qq, ii;
double m;
double lamk, phik, Ak;
double x_state[3];
double rk_ecef[3], rksat[3], rksun[3];
double rk[3];
double rk_norm, rksat_norm, rksun_norm;
double uk[3], uksat[3], uksun[3];
double sinEksat, sinEksun, Pkalb, Pkir;
double cos_theta, e_coeff, n_coeff;
double CfA_alb[3]; 
double CfA_ir[3];
double a_erpx, a_erpy, a_erpz;
double a_erpx_plus, a_erpy_plus, a_erpz_plus;
double a_erpx_minus, a_erpy_minus, a_erpz_minus;

/////// Functions /////////////
void calc_earthrad(void);


/******** MAIN PROGRAM *****************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {

    /* make sure there are six inputs to the function */
    if(nrhs!=19) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "Nineteen inputs required.");
    }  
     /* make sure there are six inputs to the function */
    if(nlhs!=2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "Two outputs required.");
    }  

    /* make sure the first input argument is a row vector of 3 components of type double*/
    // sunvector_j2000
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetM(prhs[0])!=1 ||
         mxGetN(prhs[0])!=3) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "First input must be a row vector of 3 components of type double.");
    }

    /* make sure the second input argument is of type double*/
    // nrm_eci
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Second input matrix must be of type double.");
    }

    /* make sure the third input argument is a row vector of type double*/
    // plt_area
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ||
         mxGetM(prhs[2])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Third input must be a row vector of type double.");
    }    

     /* make sure the fourth input argument is a scalar of type double*/
     // scmass
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3]) ||
         mxGetNumberOfElements(prhs[3])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Fourth input must be a scalar of type double.");
    }    

    /* make sure the fifth input argument is a row vector of type double*/
    // refl_vis_spec
    if( !mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4]) ||
         mxGetM(prhs[4])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Fifth input must be a row vector of type double.");
    }   

    /* make sure the sixth input argument is a row vector of type double*/
    // refl_vis_diff
    if( !mxIsDouble(prhs[5]) || 
         mxIsComplex(prhs[5]) ||
         mxGetM(prhs[5])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Sixth input must be a row vector of type double.");
    }  

    /* make sure the seventh input argument is a row vector of type double*/
    // refl_ir_spec
    if( !mxIsDouble(prhs[6]) || 
         mxIsComplex(prhs[6]) ||
         mxGetM(prhs[6])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Seventh input must be a row vector of type double.");
    }  

    /* make sure the eighth input argument is a row vector of type double*/
    // refl_ir_diff
    if( !mxIsDouble(prhs[7]) || 
         mxIsComplex(prhs[7]) ||
         mxGetM(prhs[7])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Eighth input must be a row vector of type double.");
    }  

    /* make sure the ninth input argument is a 3x3 matrix of type double*/
    // ECEFtoJ2000
    if( !mxIsDouble(prhs[8]) || 
         mxIsComplex(prhs[8]) ||
         mxGetM(prhs[8])!=3 || 
         mxGetN(prhs[8])!=3 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Ninth input must be a 3x3 matrix of type double.");
    } 

    /* make sure the tenth input argument is a matrix of type double*/
    // Psw
    if( !mxIsDouble(prhs[9]) || 
         mxIsComplex(prhs[9]) )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Tenth input matrix must be of type double.");
    } 

    /* make sure the eleventh input argument is a matrix of type double*/
    // Plw
    if( !mxIsDouble(prhs[10]) || 
         mxIsComplex(prhs[10]) ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Eleventh input matrix must be of type double.");
    } 

    /* make sure the twelfth input argument is a row vector of type double*/
    // lam
    if( !mxIsDouble(prhs[11]) || 
         mxIsComplex(prhs[11]) ||
         mxGetM(prhs[11])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Twelfth input must be a row vector of type double.");
    } 

    /* make sure the thirteenth input argument is a row vector of type double*/
    // phi
    if( !mxIsDouble(prhs[12]) || 
         mxIsComplex(prhs[12]) ||
         mxGetM(prhs[12])!=1 ) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Thirteenth input must be a row vector of type double.");
    } 

    /* make sure the fourteenth input argument is a scalar of type double*/
    // Cr_erp
    if( !mxIsDouble(prhs[13]) || 
         mxIsComplex(prhs[13]) ||
         mxGetNumberOfElements(prhs[13])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Fourteenth input must be a scalar of type double.");
    }  

    /* make sure the fifteenth input argument is a scalar of type double*/
    // Re
    if( !mxIsDouble(prhs[14]) || 
         mxIsComplex(prhs[14]) ||
         mxGetNumberOfElements(prhs[14])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Fifteenth input must be a scalar of type double.");
    } 
    /* make sure the sixteenth input argument is a row vector of type double*/
    // X_state
    if( !mxIsDouble(prhs[15]) || 
         mxIsComplex(prhs[15]) ||
         mxGetM(prhs[15])!=1  ||
         mxGetN(prhs[15])!=3 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Sixteenth input must be a row vector of dimension 3 of type double.");
    } 

    /* make sure the seventeenth input argument is a scalar of type double*/
    // Mphi 
    if( !mxIsDouble(prhs[16]) || 
         mxIsComplex(prhs[16]) ||
         mxGetNumberOfElements(prhs[16])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Seventeenth input must be a scalar of type double.");
    } 
    
    /* make sure the eighteenth input argument is a scalar of type double*/
    // Nlam
    if( !mxIsDouble(prhs[17]) || 
         mxIsComplex(prhs[17]) ||
         mxGetNumberOfElements(prhs[17])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Eighteenth input must be a scalar of type double.");
    }  
    
    /* make sure the nineteenth input argument is a scalar of type double*/
    // Kplt
    if( !mxIsDouble(prhs[18]) || 
         mxIsComplex(prhs[18]) ||
         mxGetNumberOfElements(prhs[18])!=1 )  {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                           "Nineteenth input must be a scalar of type double.");
    }       


    /* Reallocating the inputs */
    sunvector_ptr = mxGetPr(prhs[0]);
    nrm_eci_ptr   = mxGetPr(prhs[1]);
    plt_area_ptr  = mxGetPr(prhs[2]);
    scmass        = mxGetScalar(prhs[3]);
    refl_vis_spec_ptr = mxGetPr(prhs[4]);
    refl_vis_diff_ptr = mxGetPr(prhs[5]);
    refl_ir_spec_ptr  = mxGetPr(prhs[6]);
    refl_ir_diff_ptr  = mxGetPr(prhs[7]);
    ECEFtoJ2000_ptr    = mxGetPr(prhs[8]);
    Psw_ptr           = mxGetPr(prhs[9]);
    Plw_ptr           = mxGetPr(prhs[10]);
    lam_ptr           = mxGetPr(prhs[11]);
    phi_ptr           = mxGetPr(prhs[12]);
    Cr_erp            = mxGetScalar(prhs[13]);
    Re                = mxGetScalar(prhs[14]);
    x_state_ptr       = mxGetPr(prhs[15]);
    Mphi              = mxGetScalar(prhs[16]);
    Nlam              = mxGetScalar(prhs[17]);
    Kplt              = mxGetScalar(prhs[18]);

     /*///////////////////
    // -- Outputs -- //
    //////////////////*/
    
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);     // [m/s^2] ERP Acceleration
    plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);     // [m/s^2] Jacobian

    /*Acceleration*/
    
    Accel_ptr   = mxGetPr(plhs[0]);   

    /*Jacobian*/

    Jacobian_ptr = mxGetPr(plhs[1]);
    
    /*Allocation*/
    for (mm = 0; mm < Mphi; mm++) {
        for (nn = 0; nn < Nlam; nn++){
            Psw[nn][mm] = Psw_ptr[Nlam*mm+nn];
            Plw[nn][mm] = Plw_ptr[Nlam*mm+nn];
        }
    }

    for (kk = 0; kk< Kplt; kk++){
        for (mm = 0; mm<3; mm++){
        nrm_eci[mm][kk] = nrm_eci_ptr[3*kk+mm];
        }
    }

    for (mm = 0; mm <3; mm++){
        Accel_ptr[mm] = 0.0;
        x_state[mm] = x_state_ptr[mm];
        
        for (nn = 0; nn <3; nn++){
            ECEFtoJ2000[nn][mm] = ECEFtoJ2000_ptr[3*mm+nn];
            Jacobian_ptr[3*nn+mm] = 0.0;
        }
    }

    
    calc_earthrad();

    ///////////////// Assign the output ////////////////////////////
    Accel_ptr[0] = a_erpx;
    Accel_ptr[1] = a_erpy;
    Accel_ptr[2] = a_erpz;


    ///////////// Calculate Jacobian /////////////////////////////

    for (ii = 0; ii < 3; ii++){
        x_state[ii] = x_state_ptr[ii] + delp;
        calc_earthrad();

        a_erpx_plus = a_erpx;
        a_erpy_plus = a_erpy;
        a_erpz_plus = a_erpz;

        x_state[ii] = x_state_ptr[ii] - delp;
        calc_earthrad();

        a_erpx_minus = a_erpx;
        a_erpy_minus = a_erpy;
        a_erpz_minus = a_erpz;

        Jacobian_ptr[3*ii+0] = (a_erpx_plus - a_erpx_minus)/(2*delp);
        Jacobian_ptr[3*ii+1] = (a_erpy_plus - a_erpy_minus)/(2*delp);
        Jacobian_ptr[3*ii+2] = (a_erpz_plus - a_erpz_minus)/(2*delp);

    }
    
    

return;
}  /* End of main function */


////////////////////////////////////////////////////////////////
void calc_earthrad(void)
{
    
    /*/////////////////////////
    // -- Pre-allocation -- //
    /////////////////////////*/
    
    a_erpx = 0.0; a_erpy = 0.0; a_erpz = 0.0; 
    

    //////////////////// Calculating the EARTH RADIATION PRESSURE acceleration/////////////////////////////
    // loop through latitude elements
    for (mm = 0; mm < Mphi; mm++){  
        m = (double) mm;
        m = m + 1.0;
        phik = phi_ptr[mm];
        Ak = 4.0*pi*Re*Re/Nlam*sin (pi/2/Mphi)*sin ((m-0.5)*pi/Mphi);

        /* Loop through the longitude elements */
        for (nn = 0; nn < Nlam; nn++){
            lamk = lam_ptr[nn];
            rk_ecef[0] = Re*cos (lamk)*cos (phik);          // earth grid elements in ecef 
            rk_ecef[1] = Re*sin (lamk)*cos (phik);
            rk_ecef[2] = Re*sin (phik);

            for (pp = 0; pp < 3; pp++){
                rk[pp] = 0.0;
                for (qq = 0; qq < 3; qq++){
                    rk[pp] += ECEFtoJ2000[pp][qq]*rk_ecef[qq];     // convert earth grid elements to eci
                }
            }

            for (pp = 0; pp < 3; pp++){                // vector from earth grid element to satellite and sun
                rksat[pp] = x_state[pp] - rk[pp];
                rksun[pp] = sunvector_ptr[pp] - rk[pp];
            }
            
            rk_norm = sqrt (rk[0]*rk[0] + rk[1]*rk[1] + rk[2]*rk[2]);
            rksat_norm = sqrt (rksat[0]*rksat[0] + rksat[1]*rksat[1] + rksat[2]*rksat[2]);
            rksun_norm = sqrt (rksun[0]*rksun[0] + rksun[1]*rksun[1] + rksun[2]*rksun[2]);

            for (pp = 0; pp < 3; pp++){                   // unit vector elements
                uk[pp] = rk[pp]/rk_norm;
                uksat[pp] = rksat[pp]/rksat_norm;
                uksun[pp] = rksun[pp]/rksun_norm;
            }

            sinEksun = uk[0]*uksun[0] + uk[1]*uksun[1] + uk[2]*uksun[2];
            sinEksat = uk[0]*uksat[0] + uk[1]*uksat[1] + uk[2]*uksat[2];
            
            if (sinEksat < 0){ 
                sinEksat = 0.0;
            }
            if (sinEksun <= 0){
                sinEksun = 0.0;
            } else{
                sinEksun = 1.0;
            }

            Pkalb = Psw[nn][mm]*Ak*sinEksat*sinEksun/(rksat_norm*rksat_norm)/pi;
            Pkir  = Plw[nn][mm]*Ak*sinEksat/(rksat_norm*rksat_norm)/pi;
            
            // radiation coefficient vector calculation
            if (Pkalb > 0.0 || Pkir > 0.0){
                for (kk = 0; kk < 3; kk++){
                    CfA_alb[kk] = 0.0;
                    CfA_ir[kk] = 0.0;
                }
                for (kk = 0; kk < Kplt; kk++){

                    cos_theta = - (uksat[0]* nrm_eci[0][kk] + uksat[1]* nrm_eci[1][kk]+ uksat[2]* nrm_eci[2][kk]);
                    if (cos_theta < 0){
                        cos_theta = 0;
                    }
                    e_coeff = plt_area_ptr[kk]*cos_theta*(1-refl_vis_spec_ptr[kk]);
                    n_coeff = plt_area_ptr[kk]*refl_vis_spec_ptr[kk]*cos_theta*cos_theta + plt_area_ptr[kk]*refl_vis_diff_ptr[kk]*cos_theta/3;
                    CfA_alb[0] += -e_coeff*uksat[0] + 2*n_coeff*nrm_eci[0][kk];
                    CfA_alb[1] += -e_coeff*uksat[1] + 2*n_coeff*nrm_eci[1][kk];
                    CfA_alb[2] += -e_coeff*uksat[2] + 2*n_coeff*nrm_eci[2][kk];

                    e_coeff = plt_area_ptr[kk]*cos_theta*(1-refl_ir_spec_ptr[kk]);
                    n_coeff = plt_area_ptr[kk]*refl_ir_spec_ptr[kk]*cos_theta*cos_theta + plt_area_ptr[kk]*refl_ir_diff_ptr[kk]*cos_theta/3;
                    CfA_ir[0] += -e_coeff*uksat[0] + 2*n_coeff*nrm_eci[0][kk];
                    CfA_ir[1] += -e_coeff*uksat[1] + 2*n_coeff*nrm_eci[1][kk];
                    CfA_ir[2] += -e_coeff*uksat[2] + 2*n_coeff*nrm_eci[2][kk];
                }
            }
            a_erpx -= Cr_erp*(CfA_alb[0]*Pkalb + CfA_ir[0]*Pkir)/scmass;
            a_erpy -= Cr_erp*(CfA_alb[1]*Pkalb + CfA_ir[1]*Pkir)/scmass;
            a_erpz -= Cr_erp*(CfA_alb[2]*Pkalb + CfA_ir[2]*Pkir)/scmass;
        }
    }

    return;

}




    







