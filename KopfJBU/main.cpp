#include "head.h"
#include "function.h"

Mat RefImage,SrcImage;
int Imagescale;

Mat SrcDown, SrcUp;
int WinSize;

int height,width;
int s_height,s_width;

int main(int argc,char * argv[])
{
	//Input src image and verify 
	RefImage    = imread(argv[1], 0);
	SrcImage    = imread(argv[2], 0); 
    Imagescale  = atoi(argv[3]);      
	
	if (argc != 6)
	{
		cout<<"Error!"<<endl;
		cout<<"Please input the file again!"<<endl;
		return EXIT_FAILURE;
	}

	height = SrcImage.rows;
	width  = SrcImage.cols;
	s_height = height/Imagescale;
	s_width  = width /Imagescale;

	Mat DepthImage = SrcImage.clone();
	SrcDown = Mat::zeros(s_height, s_width, DepthImage.type());
	resize(DepthImage, SrcDown, SrcDown.size(), 0, 0, INTER_NEAREST);  

	/*
	*sigmad = 0.5
	*sigmar = 0.1 (Notice that if you use 8-bit image that has the 0~255 intensity values, varC must be multiplied by 255)
	*window = scale*scale +1 (in high resolution, ex. scale 2x --> window = 5)
	*so:
	*2x window = 5 ;
	*4x window = 17 ;
	*8x window = 65 ;
	*16x window = 257 ;
	*downsample: nearest neighbor
	*/

    WinSize = Imagescale * Imagescale + 1;
	SrcUp = Mat::zeros(DepthImage.size(), DepthImage.type());

	//time for joint bilateral upsampling
	double start_time, end_time; 
	start_time = clock();

	JBU(SrcDown, RefImage, SrcUp, WinSize);

	end_time = (clock() - start_time) / CLOCKS_PER_SEC; 
	cout<<"Imagescale #: "<<Imagescale<<endl;
	cout<<"The runtime of Kopf JBU is : "<<end_time<<"s."<<endl;

	imwrite(argv[4],SrcUp);

#pragma region Criterion BPR MSE RMSE PSNR

	Mat BadImage = Mat::zeros(SrcUp.size(),SrcUp.type());
	double BRP = 0.0,Mse = 0.0 , Rmse = 0.0, Psnr = 0.0;
	double blackcnt = 0.0,cnt = 0.0;
	int dv = 0;
	long sum = 0;

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if (SrcImage.at<uchar>(j,i) != 0)
			{
				//MSE RMSE PSNR
				dv = abs(SrcImage.at<uchar>(j,i) - SrcUp.at<uchar>(j,i));
				dv = pow(dv,2);
				sum += dv;
				//BPR
				if (dv > 1)
				{
					BadImage.at<uchar>(j,i) = 0;
					cnt++;
				}
				else
				{
					BadImage.at<uchar>(j,i) = 255;
				}
			} 
			else
			{
				BadImage.at<uchar>(j,i) = 255;
				blackcnt++;
			}
		}
	}

	BRP = 1.0 * cnt / (height * width - blackcnt)*100;
	Mse = 1.0 * sum / (height*width - blackcnt);
	Rmse = sqrt(Mse);
	Psnr = 10 * log10(255*255/(Mse));

	cout<<"Kopf JBU Imagescale "<<Imagescale<<" BPR  : "<<BRP<<"%"<<endl;
	cout<<"Kopf JBU Imagescale "<<Imagescale<<" MSE  : "<<Mse<<endl;
	cout<<"Kopf JBU Imagescale "<<Imagescale<<" RMSE : "<<Rmse<<endl;
	cout<<"Kopf JBU Imagescale "<<Imagescale<<" PSNR : "<<Psnr<<endl;
	cout<<endl<<endl;
	imwrite(argv[5],BadImage);
#pragma endregion

	return EXIT_SUCCESS;
}