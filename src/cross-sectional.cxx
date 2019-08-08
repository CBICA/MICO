/**
 * @file  cross-sectional.cxx
 * @brief Implements MICO for cross-sectional MR brain studies.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>

//#include <basis/assert.h> // assert(), ASSERT()
//#include <basis/except.h> // BASIS_THROW, std::invalid_argument

#include "mico/utilities.h"
#include "mico/cross-sectional.h"


// acceptable in .cxx file
using namespace std;


namespace mico {


// ===========================================================================
// API
// ===========================================================================

// ---------------------------------------------------------------------------
bool segment_and_bias_correct(Image*            image,
                              ImageMap*         memberships,
                              Image**           bcimage,
                              const Parameters& params,
                              int               verbosity)
{
    bool ok = true;

    // -------------------------------------------------------------------
    // check arguments
    if (image == NULL)
    {
      std::cerr << "Image is not defined, exiting.\n";
      abort();
    }
    
    if (memberships == NULL && bcimage == NULL) 
    {
      std::cerr << "segment_and_bias_correct(): "
                << " Neither output parameter is a non-NULL pointer.\n";
    }

    // -------------------------------------------------------------------
    // image attributes
    const int width  = image->hdr.dim[1];
    const int height = image->hdr.dim[2];
    const int slices = image->hdr.dim[3];
    const int total  = width * height * slices;

    // -------------------------------------------------------------------
    // reset output
    if (memberships) delete_images(*memberships);

    if (bcimage && *bcimage) {
        delete_image(*bcimage);
        *bcimage = NULL;
    }

    // -------------------------------------------------------------------
    // convert input image to float datatype (if required)
    Image* flimage = image;
    if (flimage->hdr.datatype != DT_FLOAT) {
        flimage = cast_image(image, DT_FLOAT);
        if (flimage == NULL) {
            fprintf(stderr, "Not enough memory available to convert image to DT_FLOAT!\n");
            return false;
        }
    }

    // -------------------------------------------------------------------
    // determine foreground region
    if (verbosity > 0) {
        printf("Determining region-of-interest...");
        fflush(stdout);
    }
    ImageRegion foreground = get_foreground_region(flimage, params.th_bg);
    if (verbosity > 0) printf(" - done");
    if (verbosity > 1) printf(": offset = [%d, %d, %d], size = [%d, %d, %d]",
                              foreground.ox, foreground.oy, foreground.oz,
                              foreground.nx, foreground.ny, foreground.nz);
    if (verbosity > 0) {
        printf("\n");
        fflush(stdout);
    }

    if (foreground.nx == 0 || foreground.ny == 0 || foreground.nz == 0) {
        fprintf(stderr, "WARNING: Image %s has no intensity value\n", image->name.c_str());
        fprintf(stderr, "WARNING: higher than the threshold %f! It will be segmented as background only!\n", params.th_bg);
    }

    // -------------------------------------------------------------------
    // calculate b, C, and M
    if (verbosity > 0) {
        printf("Segmenting brain image...");
        fflush(stdout);
    }

    // pointer to image data
    float* const img    = flimage->img.fl;
    float*       result = (float*)malloc(sizeof(float) * (total * 4 + 3));

    memset(result, 0, (total * 4 + 3) * sizeof(float));
    float* bias = get_basis_order(foreground.nx, foreground.ny, foreground.nz);
    updateBCM(foreground, bias, img, width, height, slices, params, result);

    free(bias);
    bias = NULL;

    // shield background for b
    for (int j = 0; j < total; j++) {
        if (img[j] <= params.th_bg) result[j] = 0.0;
    }
    // shield background for all M
    float* p = NULL;
    for (int j = total + 3; j < 4 * total + 3; j++, p++) {
        if (((j - 3) % total) == 0) p = img;
        if (*p <= params.th_bg) result[j] = 0.0;
    }

    // rescale
    float r_rescale = 200.0f / result[total + 2];
    for (int j = 0; j < total; ++ j) result[j] /= r_rescale;

    // extract membership functions from result
    if (memberships) {
        (*memberships)["CSF"] = new_image(width, height, slices, DT_FLOAT);
        (*memberships)["GM"]  = new_image(width, height, slices, DT_FLOAT);
        (*memberships)["WM"]  = new_image(width, height, slices, DT_FLOAT);

        for (ImageMap::iterator it = memberships->begin(); it != memberships->end(); ++it) {
            Image* membership = it->second;

            if (membership == NULL) {
                fprintf (stderr, "Failed to allocate memory for membership function of %s!\n", it->first);
                ok = false;
                break;
            }

            membership->hdr          = image->hdr;
            membership->hdr.bitpix   = 8 * sizeof(float);
            membership->hdr.datatype = DT_FLOAT;
            membership->compress     = image->compress;
            membership->nifti_type   = image->nifti_type;

            string descrip = "fuzzy ";
            descrip += it->first;
            descrip += " membership (MICO)";
            strncpy(membership->hdr.descrip, descrip.c_str(), 80);
        }

        if (ok) {
            float* csf = (*memberships)["CSF"]->img.fl;
            float* gm  = (*memberships)["GM"] ->img.fl;
            float* wm  = (*memberships)["WM"] ->img.fl;

            for (int j = 0; j < total; j++) {
                *(csf++) = result[j +     total + 3];
                *(gm++)  = result[j + 2 * total + 3];
                *(wm++)  = result[j + 3 * total + 3];
            }
        }
    }

    if (verbosity > 0) printf (ok ? " - done\n" : " - failed\n");

    // -------------------------------------------------------------------
    // calculate bias corrected image
    if (ok && bcimage) {
        if (verbosity > 0) {
            printf("Bias correcting image...");
            fflush(stdout);
        }

        // initialize image structure
        Image* bcflimage = new_image(width, height, slices, DT_FLOAT);
        if (bcflimage == NULL) {
            fprintf(stderr, "Failed to allocate memory for bias corrected image!\n");
            ok = false;
        } else {
            // set meta-data of bias corrected image
            bcflimage->hdr        = flimage->hdr;
            bcflimage->compress   = flimage->compress;
            bcflimage->nifti_type = flimage->nifti_type;
            // generate bias corrected image
            float* bcimg = bcflimage->img.fl;
            for (int j = 0; j < total; j++) {
                if (result[j] != 0.0f && img[j] > params.th_bg) {
                    bcimg[j] = img[j] / result[j];
                } else {
                    bcimg[j] = 0.0f;
                }
            }
            // rescale bias corrected image and cast to datatype of input image
            float min, max;
            get_intensity_range(image, min, max);
            if (image->hdr.datatype == DT_FLOAT) {
                scale_image(bcflimage, min, max);
                (*bcimage) = bcflimage;
            } else {
                (*bcimage) = cast_image(bcflimage, image->hdr.datatype, min, max);
                if ((*bcimage) == NULL) {
                    fprintf(stderr, "Failed to allocate memory for bias corrected image!\n");
                    ok = false;
                }
                delete_image(bcflimage);
            }
            // overwrite data scaling information in image header
            (*bcimage)->hdr.scl_slope = image->hdr.scl_slope;
            (*bcimage)->hdr.scl_inter = image->hdr.scl_inter;
        }

        if (verbosity > 0) {
            printf(ok ? " - done\n" : " - failed\n");
            fflush(stdout);
        }
    }

    // -------------------------------------------------------------------
    // clean up
    if (flimage != image) delete_image(flimage);
    free(result);
    result = NULL;

    return ok;
}

// ---------------------------------------------------------------------------
Image* create_label_map(const Image*      image,
                        const ImageMap&   memberships,
                        const Parameters& params,
                        int               verbosity)
{
    if (memberships.size () != 3) 
    {
      std::cerr << "create_label_map(): Expected three membership functions!\n";
      abort();
    }

    const unsigned int total = image->region.nx * image->region.ny * image->region.nz;
    ImageMap::const_iterator it_csf = memberships.find("CSF");
    ImageMap::const_iterator it_gm  = memberships.find("GM");
    ImageMap::const_iterator it_wm  = memberships.find("WM");

    if (it_csf == memberships.end() || it_csf->second == NULL) 
    {
      std::cerr << "create_label_map(): Missing membership function for CSF!\n";
      abort();
    }
    if (it_wm == memberships.end() || it_wm->second == NULL) 
    {
      std::cerr << "create_label_map(): Missing membership function for WM!\n";
      abort();
    }
    if (it_gm == memberships.end() || it_gm->second == NULL) 
    {
      std::cerr << "create_label_map(): Missing membership function for GM!\n";
      abort();
    }

    Image* segmentation = new_image(image->region.nx,
                                    image->region.ny,
                                    image->region.nz,
                                    DT_UNSIGNED_CHAR);
    if (segmentation == NULL) return NULL;

    segmentation->hdr          = image->hdr;
    segmentation->hdr.bitpix   = 8;
    segmentation->hdr.datatype = DT_UNSIGNED_CHAR;
    segmentation->compress     = image->compress;

    float*         csf   = it_csf->second->img.fl;
    float*         gm    = it_gm ->second->img.fl;
    float*         wm    = it_wm ->second->img.fl;
    unsigned char* label = segmentation->img.uc;

    for (unsigned int j = 0; j < total; j++, csf++, gm++, wm++, label++) {
        if ((*csf > *gm) && (*csf > *wm)) {
            *label = params.labelCSF;
        } else if((*gm > *csf) && (*gm > *wm)) {
            *label = params.labelGM;
        } else if((*wm > *csf) && (*wm > *gm)) {
            *label = params.labelWM;
        }
    }

    return segmentation;
}

// ---------------------------------------------------------------------------
Image* segment(Image*            image,
               const Parameters& params,
               int               verbosity)
{
    Image*   segmentation= NULL;
    ImageMap membership;
    if (segment_and_bias_correct(image, &membership, NULL, params, verbosity)) {
        segmentation = create_label_map(image, membership, params, verbosity);
    }
    delete_images(membership);
    return segmentation;
}

// ---------------------------------------------------------------------------
Image* bias_correct(Image*            image,
                    const Parameters& params,
                    int               verbosity)
{
    Image* bcimage = NULL;
    segment_and_bias_correct(image, NULL, &bcimage, params, verbosity);
    return bcimage;
}

// ===========================================================================
// orthogonal basis function
// ===========================================================================

// ---------------------------------------------------------------------------
float* get_basis_order(int width, int height, int slices)
{
	 int i, j, k, i1, i3;

	 float* x1 = (float*)calloc(sizeof(float), width);
	 float* y1 = (float*)calloc(sizeof(float), height);
	 float* z1 = (float*)calloc(sizeof(float), slices);
	 for (i=0;i<width;i++)
	 {
		 *(x1+i)=-1.0+i*(2.0/(width-1.0));
	 }
	 for (j=0;j<height;j++)
	 {
		 *(y1+j)=-1.0+j*(2.0/(height-1.0));
	 }
	 for (k=0;k<slices;k++)
	 {
		 *(z1+k)=-1.0+k*(2.0/(slices-1.0));
	 }

	 float* x = (float *)calloc(sizeof(float), width*height*slices);
	 float* y = (float *)calloc(sizeof(float), width*height*slices);
	 float* z = (float *)calloc(sizeof(float), width*height*slices);

	 int h,l,r;
	 for (i1=0;i1<width*height*slices;i1++)
	 {
		 l=(i1/width)%height;
		 h=i1/(width*height);
		 r=i1%width;
		 *(z+i1)=*(z1+h);
		 *(x+i1)=*(x1+r);
		 *(y+i1)=*(y1+l);

	 }

	 float *bias=NULL;
	 bias=(float *)malloc(20*width*height*slices*sizeof(float));

     for(i3=0;i3<width*height*slices;i3++)
	 {
		 bias[i3]                                = 1.0; 
		 bias[ 1 * width * height * slices + i3] = x[i3]; 
		 bias[ 2 * width * height * slices + i3] = (3.0*x[i3]*x[i3]-1.0)/2.0; 
		 bias[ 3 * width * height * slices + i3] = y[i3]; 
         bias[ 4 * width * height * slices + i3] = x[i3]*y[i3];  
		 bias[ 5 * width * height * slices + i3] = (3.0*y[i3]*y[i3]-1.0)/2.0;               
		 bias[ 6 * width * height * slices + i3] = z[i3];
		 bias[ 7 * width * height * slices + i3] = (3.0*z[i3]*z[i3]-1.0)/2.0; 
		 bias[ 8 * width * height * slices + i3] = x[i3]*z[i3];
		 bias[ 9 * width * height * slices + i3] = y[i3]*z[i3]; 
		 bias[10 * width * height * slices + i3] = (5.0*z[i3]*z[i3]*z[i3] - 3.0*z[i3])/2.0;  
         bias[11 * width * height * slices + i3] = y[i3]*(3.0*z[i3]*z[i3] - 1.0)/2.0;  
         bias[12 * width * height * slices + i3] = (5.0*y[i3]*y[i3]*y[i3] - 3.0*y[i3])/2.0;  
		 bias[13 * width * height * slices + i3] = x[i3]*(3.0*z[i3]*z[i3]-1.0)/2.0;
		 bias[14 * width * height * slices + i3] = x[i3]*y[i3]*z[i3];
		 bias[15 * width * height * slices + i3] = x[i3]*(3.0*y[i3]*y[i3]-1.0)/2.0;
         bias[16 * width * height * slices + i3] = (3.0*x[i3]*x[i3]-1.0)/2.0*z[i3]; 
		 bias[17 * width * height * slices + i3] = (3.0*x[i3]*x[i3]-1.0)/2.0*y[i3];
         bias[18 * width * height * slices + i3] = (5.0*x[i3]*x[i3]*x[i3]-3.0*x[i3])/2.0; 
		 bias[19 * width * height * slices + i3] = (3.0*y[i3]*y[i3]-1.0)/2.0*z[i3];  
	 }
	
	 // normalization
	 float sum=0.0f;
	 for (i=0;i<20;i++)
	 {
		 sum=0.0f;
		 for (j=0;j<width*height*slices;j++)
		 {
			sum+=pow(bias[i*width*height*slices+j],2);
		 }
		 sum=sqrt(sum);
		 for (j=0;j<width*height*slices;j++)
		 {
			 bias[i*width*height*slices+j]=bias[i*width*height*slices+j]/sum;
		 }

	 }

	 free(x1);
     free(y1);
	 free(z1);
	 free(x);
	 free(y);
	 free(z);
	
	 return bias;
}

// ===========================================================================
// update functions
// ===========================================================================

// ---------------------------------------------------------------------------
void updateMCB(float *img, float *ROI, float *M, float *C, float *b, float *bias, int total, const Parameters &params)
{
	int n,kk;
	for (n=1;n<=params.iterMCB;n++)
	{
		for (kk=1;kk<=params.iterMC;kk++)
		{
			updateC(img, ROI, M, C, b, total, params);
			updateM(img, ROI, M, C, b, total, params);
		}
		updateB(img, ROI, M, C, b, bias, total, params);
	}
}
 
// ---------------------------------------------------------------------------
float *updateB(float *img, float * /*ROI*/, float *M, float *C, float *b, float *bias, int total, const Parameters &params) 
{
	int kk,j,i;
	float *PC2=NULL,*PC=NULL;
	PC2=(float *)malloc(total*sizeof(float));
	memset(PC2,0,total*sizeof(float));
	PC=(float *)malloc(total*sizeof(float));
	memset(PC,0,total*sizeof(float));
	for(kk=0;kk<3;kk++)
	{
		for (j=0;j<total;j++)
		{
			PC2[j]+=pow(C[kk],2)*pow(M[kk*total+j],params.q)*params.tissueWeight[kk];
			PC[j]+=C[kk]*pow(M[kk*total+j],params.q)*params.tissueWeight[kk];
		}
	}
	 float *V=NULL;
	 V=(float *)malloc(20*sizeof(float));
	 memset(V,0,20*sizeof(float));
     int ii,iii,jj,jjj;
	    
	 float *ImgG_PC=NULL;
	 ImgG_PC=(float *)malloc(total*sizeof(float));
	 float *A=NULL;
     A=(float *)malloc(20*20*sizeof(float));
	 memset(A,0,20*20*sizeof(float));
	 float *B_out=NULL;
	 B_out=(float *)malloc(total*sizeof(float));
	
	  for(ii=0;ii<20;ii++)
	  {
		  for(iii=0;iii<total;iii++)
		  {
			  ImgG_PC[iii]=img[iii]*bias[ii*total+iii]*PC[iii];	
			  V[ii]+=ImgG_PC[iii];
		  }

		 for (jj=ii;jj<20;jj++)
		 {
 			 for(jjj=0;jjj<total;jjj++)
			 {
                  B_out[jjj]= bias[ii*total+jjj]*bias[jj*total+jjj]*PC2[jjj];
				  A[ii*20+jj]+=B_out[jjj];
			 }
			 A[jj*20+ii]=A[ii*20+jj];
		 }
	  }
	 free (PC);
	 free (PC2);
	 free (B_out);
	 free (ImgG_PC);

	 float * coeff=NULL;
	 coeff=(float *)malloc(20*sizeof(float));
	 memset(coeff,0,20*sizeof(float));
     invert_matrix(A, A, 20);	
     for (i=0;i<20;i++)
	 {
		 for (j=0;j<20;j++)
		 {
              coeff[i]+=A[i*20+j]*V[j];
		 }
	 }

	 memset(b,0,total*sizeof(float));

	 for(i=0;i<total;i++)
	 {
	      for (j=0;j<20;j++)
	      {
			  b[i]+=coeff[j]*bias[j*total+i];
	      }
	 }

	 free(coeff);
	 free (A);
	 free (V);
	 return b;
}

// ---------------------------------------------------------------------------
float *updateC(float *Img, float *ROI, float *M, float *C, float *b, int total, const Parameters &params)
{
	float * N=NULL;
    N=(float *)malloc(total*sizeof(float));
	float sN=0.0;
	float sD=0.0;
	float * D=NULL;
	D=(float *)malloc(total*sizeof(float));

	int N_class,i;
	for (N_class=0;N_class<3;N_class++)
	{
		for(i=0;i<total;i++)
		{
			N[i]=(b[i]+params.lambda)*Img[i]*pow(M[N_class*total+i],params.q)*ROI[i];
			D[i]=((b[i]*b[i])+params.lambda)*pow(M[N_class*total+i],params.q)*ROI[i]; 
			sN+=N[i];
			sD+=D[i];
		}

        if (sD==0)
		{
				C[N_class]=sN;
		}
		else
		{
				C[N_class]=sN/sD;
		}

	
		sN=0.0;
		sD=0.0;
	}
	free (N);
	free (D);

	return C;
}

// ---------------------------------------------------------------------------
float *updateM(float *Img, float *ROI, float *M, float *C, float *b, int total, const Parameters &params)
{
	 float * e=NULL;
	 e=(float *)malloc(3*total*sizeof(float));
	 int i,j;
     for(i=0;i<3;i++)
	 {
	     for (j=0;j<total;j++)
		 {
			 e[i*total+j]=((Img[j]-C[i]*b[j])*(Img[j]-C[i]*b[j])+params.lambda*(Img[j]-C[i])*(Img[j]-C[i]))*params.tissueWeight[i];
		 }

	 }

	 float * f_sum=NULL;
	 f_sum=(float *)malloc(total*sizeof(float));
	 memset(f_sum,0,total*sizeof(float));
	 
	 if(params.q>1)
	 {
         const float epsilon = 1.0e-12;
		 for (i=0;i<3*total;i++)
		 {
			 e[i]=e[i]+epsilon;  // avoid division by zero
		 }
		 float p = -1.0/(params.q-1.0);

		 for (j=0;j<3;j++)
		 {
			 for (i=0;i<total;i++)
			 {
                  f_sum[i]+=pow(e[j*total+i],p);
			 }
		 }

		 for (i=0;i<3;i++)
		 { 
			 for (j=0;j<total;j++)
			 {
				 
                  M[i*total+j]=pow(e[i*total+j],p)*ROI[j]/f_sum[j];
			 }
			 

		 }
		 
	 }
	
	 free (f_sum);
	 free (e);
	 return M;
}
 
// ---------------------------------------------------------------------------
void updateBCM (const ImageRegion &roi, float *bias, float *p, int width, int height, int slices, const Parameters &params, float *result)
{
	int i,j,k;

	float *b=NULL;
	float *C=NULL;
	float *M=NULL;
	
    //the size of croping image
	int width1=roi.nx;
	int height1=roi.ny;
	int slices1=roi.nz;

    //the number of croping image
    int total  = width  * height  * slices;  // number of voxels in original volume
	int total1 = width1 * height1 * slices1; // number of voxels in cropped volume

        //the data of croping image
	float *p_ROI = NULL;
	p_ROI=(float *)malloc(total1*sizeof(float));
	float *p_cut = NULL;
	p_cut=(float *)malloc(total1*sizeof(float));
	float *ROI_cut = NULL;							
	ROI_cut=(float *)malloc(total1*sizeof(float));

	for (k=0;k<slices1;k++)
	{
		for (i=0;i<width1;i++)
		{
			for (j=0;j<height1;j++)
			{
				p_cut[k*(width1*height1)+j*width1+i]=p[(k+roi.oz)*(width*height)+(j+roi.oy)*width+(i+roi.ox)];
			}
		}

	}

        //the background and forground of croping image
	for (i=0;i<total1;i++)
	{
		if (p_cut[i]<=params.th_bg) {
			ROI_cut[i]=0.0;
		} else {
			ROI_cut[i]=1.0;
        }
	}
	
	for (i=0;i<total1;i++)
	{
		p_ROI[i]=p_cut[i]*ROI_cut[i];
	}

	free(p_cut);
	
      
	 
	C=(float *)malloc(3*sizeof(float));
	M=(float *)malloc(3*total1*sizeof(float));
	memset(M,0,3*total1*sizeof(float));
	srand(params.seed);

        //Initializing C
        float A1;
        A1=	p_ROI[0];
        for(i=1;i<total1;i++)
        {
           if(p_ROI[i]>A1)
           {
               A1=p_ROI[i];
           }
        }
        C[0]=0.1*A1;
        C[1]=0.5*A1;
        C[2]=0.9*A1;

        //Initializing M
    for (j=0;j<total1*3;j++)
    {
          M[j]=(rand());
    }
	float *w=NULL;
	w=(float *)malloc(total1*sizeof(float));
	memset(w,0,total1*sizeof(float));

	 int l,l1;
	 for (j=0;j<total1;j++)
	 {
		 for(l=0;l<3;l++)
		 {
            w[j]=w[j]+M[l*total1+j];
		  }
	 }
	 for (j=0;j<total1*3;j++)
	 {
		 l1=j%total1;
		 M[j]=M[j]/w[l1];
	 }
	 free (w);
	 float *M_old=NULL;
	 M_old=(float *)malloc(3*total1*sizeof(float));
	
	 for (i=0;i<total1*3;i++)
	 {
          M_old[i]=M[i];
	 }
	
         //Initializing b
	 b=(float *)malloc(total1*sizeof(float));

	 for (i=0;i<total1;i++)
	 {
		 b[i]=1.0f;
	 }
	
         int n_iter;
	 
         //iteration for each 3D image
	 for (n_iter=1;n_iter<=params.iterNum;n_iter++)
	 {
		  updateMCB(p_ROI,ROI_cut,M,C,b,bias,total1,params);
          
        float *C_before=NULL;
        C_before=(float *)malloc(sizeof(float)*3);
         for(i=0;i<3;i++)
         {
             C_before[i]=C[i];
          }
         
         //order C 
	      sort(C,3);
      
	     float * M_sort1 =NULL;
	     M_sort1 =(float *)malloc(sizeof(float)*total1);
	     float * M_sort2 =NULL;
	     M_sort2 =(float *)malloc(sizeof(float)*total1);
	     float * M_sort3 =NULL;
	     M_sort3 =(float *)malloc(sizeof(float)*total1);
        
        //order M 
        for(k=0;k<total1;k++)
        {
             for(i=0;i<3;i++)
            {
                 if(C[0]==C_before[i])
                 {
                      M_sort1[k]=M[i*total1+k];
                  }
                 else if(C[1]==C_before[i])
                 {
                       M_sort2[k]=M[i*total1+k];

                 }
                 else if(C[2]==C_before[i])
                 {
                       M_sort3[k]=M[i*total1+k];

                 }

             }
         }
         
         for(i=0;i<total1;i++)
         {
             M[i]=M_sort1[i];
         }
         for(i=total1;i<2*total1;i++)
         {
             M[i]=M_sort2[i-total1];
         }
         for(i=2*total1;i<3*total1;i++)
         {
             M[i]=M_sort3[i-total1*2];
         }
         
         free(C_before);
         free(M_sort1);
         free(M_sort2);
         free(M_sort3);         
     }
     
     float * M1_sort =NULL;
	   M1_sort =(float *)malloc(sizeof(float)*total1);
	   float * M2_sort =NULL;
	   M2_sort =(float *)malloc(sizeof(float)*total1);
	   float * M3_sort =NULL;
	   M3_sort =(float *)malloc(sizeof(float)*total1);
     
     for(i=0;i<total1;i++)
     {
        M1_sort[i]=M[i];
        M2_sort[i]=M[i+total1];
        M3_sort[i]=M[i+total1*2];
     }
     	
	  float *b_out=NULL;
	  b_out=(float *)malloc(sizeof(float)*total);
         float *M1_out=NULL;
	  M1_out=(float *)malloc(sizeof(float)*total);
         float *M2_out=NULL;
	  M2_out=(float *)malloc(sizeof(float)*total);
         float *M3_out=NULL;
	  M3_out=(float *)malloc(sizeof(float)*total);
	  
	  memset(M1_out,0,total*sizeof(float));
	  memset(M2_out,0,total*sizeof(float));
	  memset(M3_out,0,total*sizeof(float));
	  memset(b_out,0,total*sizeof(float));
          
          //the bias filed b and membership function M1,M2,M3 of inital 3D image
	  for (k=0;k<slices1;k++)
	  {
		  for(i=0;i<width1;i++)
		  {
			  for(j=0;j<height1;j++)
			  {
				  b_out[(k+roi.oz)*(width*height)+(j+roi.oy)*width+(i+roi.ox)]=b[k*(width1*height1)+j*width1+i];
				  M1_out[(k+roi.oz)*(width*height)+(j+roi.oy)*width+(i+roi.ox)]= M1_sort [k*(width1*height1)+j*width1+i];
				  M2_out[(k+roi.oz)*(width*height)+(j+roi.oy)*width+(i+roi.ox)]= M2_sort [k*(width1*height1)+j*width1+i];
				  M3_out[(k+roi.oz)*(width*height)+(j+roi.oy)*width+(i+roi.ox)]= M3_sort [k*(width1*height1)+j*width1+i];
			  }
		  }
	  }

	  for(i=0;i<total;i++)
	  {
		  result[i]=b_out[i]; 
	  }
          //result[total]-result[total+2]:C--the center of CSF,GM and WM
	  for(i=total;i<(total+3);i++)
	  {
		  result[i]=C[i-total]; 
	  }
          //result[total+3]-result[total*2+2]: M1--the membership function of CSF
           for(i=total+3;i<(total*2+3);i++)
	  {
		  result[i]=M1_out[i-total-3]; 
	  }
          //result[total*2+3]-result[total*3+2]: M2--the membership function of GM
           for(i=total*2+3;i<(total*3+3);i++)
	  {
		  result[i]=M2_out[i-total*2-3]; 
	  }
          //result[total*3+3]-result[total*4+2]: M1--the membership function of WM
             for(i=total*3+3;i<(total*4+3);i++)
	  {
		  result[i]=M3_out[i-total*3-3]; 
	  }
			 
	
	  free(b);
	  free(b_out);
	  free(C);
	  free(M);
	  free(M_old);
	  free(p_ROI);
	  free(ROI_cut);
      free(M1_sort);
      free(M2_sort); 
      free(M3_sort);
      free(M1_out);
      free(M2_out);
      free(M3_out);
}


} // namespace mico
