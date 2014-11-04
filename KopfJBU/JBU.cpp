#include "head.h"
#include "function.h"

extern double Imagescale;   
extern int height,width;
extern int s_height,s_width;

inline double SpatialGauss(double x1, double y1, int x2, int y2, double sigma, double mu = 0.0)
{
	double dis = pow(x1 - x2, 2) + pow(y1 - y2, 2) - mu;     
	return exp(-1.0 * dis / (2 * sigma * sigma)); 
}

inline double RangeGauss(int x, double sigma, double mu = 0.0)
{
	double x_p = x - mu;
	return exp(-1.0 * (x_p * x_p) / (2 * sigma * sigma));
}

Mat& JBU(Mat& source, Mat &refIm, Mat &dest, int WinWidth)
{
	const double scale  = 1.0*s_width/width;
	const double sigmad = 0.50;
	const double sigmar = 25.5;
	int num_neighbors = WinWidth/2;   

	uchar  refPix = 0, srcPix = 0, neighborPix = 0;    	
	int    r_y    = 0,  r_ys = 0, r_x  = 0, r_xs = 0;
	double o_y    = 0.0, o_x = 0.0;                     
	double sgauss = 0.0, rgauss = 0.0, totalgauss = 0.0;
	double total_val = 0.0, normalizing_factor = 0.0;
	uchar *srcPixLoc = nullptr;
	uchar *neighborPixLoc = nullptr;

	for (int y = 0; y < height; y++)
	{
		uchar *refPixLoc = refIm.ptr<uchar>(y);
		uchar *destLoc   = dest.ptr<uchar>(y);

		for (int x = 0; x < width; x++)
		{
			refPix = refPixLoc[x];
			o_y    = y * scale;
			o_x    = x * scale;

			total_val = 0.0, normalizing_factor = 0.0;

			for (int j = -num_neighbors; j <= num_neighbors; j++)
			{
				// source
				r_y = o_y + j;                             
				r_y = (r_y > 0 ? (r_y < s_height ? r_y :s_height - 1) : 0) ;
				srcPixLoc = source.ptr<uchar>(r_y);     
				// refIm
				r_ys = y + j;
				r_ys = (r_ys > 0 ? (r_ys < height ? r_ys :height - 1) : 0) ;
				neighborPixLoc = refIm.ptr<uchar>(r_ys);

				for (int i = -num_neighbors; i <= num_neighbors; i++)
				{
					// source
				    r_x = o_x + i;
					r_x = (r_x > 0 ? (r_x < s_width ? r_x : s_width - 1) : 0);
					srcPix = srcPixLoc[r_x];            
                    // refIm
					r_xs = x + i;  
					r_xs = (r_xs > 0 ? (r_xs < width ? r_xs :width - 1) : 0) ;
					neighborPix = neighborPixLoc[r_xs]; 

					sgauss = SpatialGauss(o_x, o_y, r_x, r_y, sigmad);     
					rgauss = RangeGauss(abs(refPix - neighborPix), sigmar);
					totalgauss = sgauss * rgauss;
					normalizing_factor += totalgauss; 
					total_val += srcPix * totalgauss; 				

				}//end for i
			}//end for j
			
			destLoc[x] = ceil(total_val/normalizing_factor);

		}//end for x
	}//end for y

	return dest;
}