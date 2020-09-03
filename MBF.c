#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
# define H 1080
# define W 1920
# define T 150
# define FilterLen 16
# define url_in "C:/Users/Cambricon/Desktop/yuvdataset/testdata_frame/HisiYUV_1920x1080_8bits_420sp_linear_600x_TR1_building.yuv"
# define url_out "C:/Users/Cambricon/Desktop/yuvdataset/output_600x_TR1_building.yuv"

double Lo_D[16] = {-0.00011747678412476953, 0.0006754494064505693, -0.00039174037337694705, -0.004870352993451574,
0.008746094047405777, 0.013981027917398282, -0.044088253930794755, -0.017369301001807547,
0.12874742662047847, 0.0004724845739132828, -0.2840155429615469, -0.015829105256349306,
0.5853546836542067, 0.6756307362972898, 0.31287159091429995, 0.05441584224310401};
double Hi_D[16] = {-0.05441584224310401, 0.31287159091429995, -0.6756307362972898, 0.5853546836542067,
0.015829105256349306, -0.2840155429615469, -0.0004724845739132828, 0.12874742662047847,
0.017369301001807547, -0.044088253930794755, -0.013981027917398282, 0.008746094047405777,
0.004870352993451574, -0.00039174037337694705, -0.0006754494064505693, -0.00011747678412476953};
double Lo_R[16] = {0.05441584224310401, 0.31287159091429995, 0.6756307362972898, 0.5853546836542067,
-0.015829105256349306, -0.2840155429615469, 0.0004724845739132828, 0.12874742662047847,
-0.017369301001807547, -0.044088253930794755, 0.013981027917398282, 0.008746094047405777,
-0.004870352993451574, -0.00039174037337694705, 0.0006754494064505693, -0.00011747678412476953};
double Hi_R[16] = {-0.00011747678412476953, -0.0006754494064505693, -0.00039174037337694705, 0.004870352993451574,
0.008746094047405777, -0.013981027917398282, -0.044088253930794755, 0.017369301001807547,
0.12874742662047847, -0.0004724845739132828, -0.2840155429615469, 0.015829105256349306,
0.5853546836542067, -0.6756307362972898, 0.31287159091429995, -0.05441584224310401};

int bilateral_filter(double **image, int img_height, int img_width, int winsize, double sigmaColor, double sigmaSpace, double **image_tran)
{
/**
 * @param image	      待滤波图像数组的二级指针
 * @param img_height  待滤波图像的高
 * @param img_width	  待滤波图像的宽
 * @param winsize 	  滤波窗口的半径
 * @param sigmaColor  色彩滤波器超参
 * @param sigmaSpace  空间滤波器超参
 * @param image_tran  双边滤波后图像数组的二级指针
 * @func              实现单通道图像的双边滤波功能
 */

    // 初始化超参数
    double color_coeff = -0.5 / (sigmaColor * sigmaColor);
    double space_coeff = -0.5 / (sigmaSpace * sigmaSpace);

    // 计算直径mask，初始化空间滤波器的存储数组weight_space，坐标数组weight_space_row、weight_space_col
    int mask = (2*winsize+1)*(2*winsize+1);
    double *weight_space = (double *)malloc(mask*sizeof(double));
    double *p = weight_space;
    int *weight_space_row = (int *)malloc(mask*sizeof(int));
    int *pr = weight_space_row;
    int *weight_space_col = (int *)malloc(mask*sizeof(int));
    int *pc = weight_space_col;

    // 计算高斯核的参数
    for (int i=-winsize; i<winsize+1; i++){
        for (int j=-winsize; j<winsize+1; j++){
            int r_square = i*i + j*j;
            *p = exp(r_square * space_coeff);
            *pr = i;
            *pc = j;
            p++;
            pr++;
            pc++;
        }
    }
    
    // 进行双边滤波
    for (int row=0; row<img_height; row++){
        for (int col=0; col<img_width; col++){
            double value = 0;
            double weight = 0;
            for (int i=0; i<mask; i++){
                int m = row + weight_space_row[i];
                int n = col + weight_space_col[i];
                double val = 0;
                double w = 0;
                if ((m<0) || (n<0) || (m >= img_height) || (n>=img_width))
                    val=0;
                else
                    val = image[m][n];
                w = weight_space[i] * exp(color_coeff * pow(abs(val - image[row][col]), 2.0));
                value += val * w;
                weight += w;
            }
            image_tran[row][col] = (double)(value/weight);
        }
    }
    free(weight_space);
    free(weight_space_row);
    free(weight_space_col);
    return 0;
}

int DWT(double *pSrcData,int srcLen,int filterLen,double *pDstCeof)
{
/**
 * @param pSrcData	  待小波变换数组的一级指针
 * @param srcLen      待小波变换数组的长度
 * @param filterLen	  滤波器的长度
 * @param pDstCeof    待返回的小波变换数组的一级指针
 * @func              实现一维数组的db8小波变换
 */
    int exLen = (srcLen + filterLen - 1) / 2;
	int k = 0;
	double tmp = 0.0;
	for (int i = 0; i < exLen; i++)
	{
		pDstCeof[i] = 0.0;
		pDstCeof[i + exLen] = 0.0;
        for (int j = 0; j < filterLen; j++)
        {
			k = 2 * i - j + 1;
			if ((k<0) && (k >= -filterLen + 1))//左边沿拓延
				tmp = pSrcData[-k - 1];
			else if ((k >= 0) && (k <= srcLen - 1))//保持不变
				tmp = pSrcData[k];
			else if ((k>srcLen - 1) && (k <= (srcLen + filterLen - 2)))//右边沿拓延
				tmp = pSrcData[2 * srcLen - k - 1];
			else
				tmp = 0.0;
			pDstCeof[i] += Lo_D[j] * tmp;
			pDstCeof[i + exLen] += Hi_D[j] * tmp;
        }
    }
    return 0;
}

int DWT_2D(double **pSrcImage,int height,int width,int filterLen,double **pDstCeof)
{
/**
 * @param pSrcImage	  待小波变换的单通道图像的二级指针
 * @param height      单通道图像高
 * @param width       单通道图像宽
 * @param filterLen	  滤波器的长度
 * @param pDstCeof    返回的小波变换单通道图像的二级指针
 * @func              实现二维数组的db8小波变换
 */
	int exwidth = (width + filterLen - 1) / 2 * 2;//pImagCeof的宽度
	int exheight = (height + filterLen - 1) / 2 * 2;//pImagCeof的高度

	double *tempImage = (double *) malloc(exwidth*height * sizeof(double)); 

	// 对每一行进行行变换
	double *tempAhang = (double *) malloc(width * sizeof(double)); // 临时存放每一行的未处理数据
	double *tempExhang = (double *) malloc(exwidth * sizeof(double));// 临时存放每一行的处理数据

	for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            tempAhang[j] = pSrcImage[i][j];//提取每一行的数据
        }

		DWT(tempAhang, width, filterLen, tempExhang);

		for (int j = 0; j < exwidth; j++)
			tempImage[i*exwidth + j] = tempExhang[j];
    }

	// 对每一列进行列变换
	double *tempAlie = (double *) malloc(height * sizeof(double)); // 临时存放每一列的转置数据
	double *tempEXlie = (double *) malloc(exheight * sizeof(double)); // 临时存放每一列的处理数据

	for (int i = 0; i < exwidth; i++)
	{
		// 列转置
		for (int j = 0; j < height; j++)
			tempAlie[j] = tempImage[j*exwidth + i];//提取每一列数据

		//执行变换
		DWT(tempAlie, height, filterLen, tempEXlie);

		// 反转置
		for (int j = 0; j < exheight; j++)
			pDstCeof[j][i] = tempEXlie[j];
	}

    free(tempImage);
    free(tempAhang);
    free(tempExhang);
    free(tempAlie);
    free(tempEXlie);

    return 0;
}

int IDWT(double *pSrcCoef, int dstLen, int filterLen, double *pDstData)
{
/**
 * @param pSrcCoef  待小波逆变换的一维数组指针
 * @param dstLen    重构后数组长度
 * @param filterLen 滤波器的长度
 * @param pDstData  返回的小波逆变换一维数组指针
 * @func            实现一维数组的db8小波逆变换
 */
	int p = 0;
	int caLen = (dstLen + filterLen - 1) / 2;
	for (int i = 0; i < dstLen; i++)
    {
		pDstData[i] = 0.0;
		for (int j = 0; j < caLen; j++)
        {
			p = i - 2 * j + filterLen - 2;
            // 信号重构
            if (p<0)
                break;
			if (p<filterLen)
				pDstData[i] += Lo_R[p] * pSrcCoef[j] + Hi_R[p] * pSrcCoef[j + caLen];
        }

    }
    return 0;
}

int IDWT_2D(double **pSrcCeof, int dstHeight, int dstWidth, int filterLen, double **pDstImage)
{
/**
 * @param pSrcCeof	  待小波逆变换的单通道图像的二级指针
 * @param dstHeight   重构后单通道图像高
 * @param dstWidth    重构后单通道图像宽
 * @param filterLen	  滤波器的长度
 * @param pDstImage   返回的小波逆变换单通道图像的二级指针
 * @func              实现二维数组的db8小波逆变换
 */
	int srcHeight = (dstHeight + filterLen - 1) / 2 * 2;//pSrcCeof的高度
	int srcWidth = (dstWidth + filterLen - 1) / 2 * 2;//pSrcCeof的宽度

    double *tempAline = (double *) malloc(srcHeight * sizeof(double)); // 临时存放每一列的数据
	double *tempdstline = (double *) malloc(dstHeight * sizeof(double)); // 临时存放每一列的重构结果

	double *pTmpImage = (double *) malloc(srcWidth*dstHeight * sizeof(double));

	// 列重构
	for (int i = 0; i < srcWidth; i++)//每一列
    {
		// 列转置
		for (int j = 0; j<srcHeight; j++)
			tempAline[j] = pSrcCeof[j][i];//提取每一列

		IDWT(tempAline, dstHeight, filterLen, tempdstline);

		// 反转置
		for (int j = 0; j < dstHeight; j++)
			pTmpImage[j*srcWidth + i] = tempdstline[j];

    }

    // 行重构
	double *tempAhang = (double *) malloc(srcWidth * sizeof(double));// 临时存放每一行的数据
	double *tempdsthang = (double *) malloc(dstWidth * sizeof(double));// 临时存放每一行的重构数据
	for (int i = 0; i < dstHeight; i++)
	{
		for (int j = 0; j < srcWidth; j++)
			tempAhang[j] = pTmpImage[i*srcWidth + j];//提取每一行的数据

		IDWT(tempAhang, dstWidth, filterLen, tempdsthang);

		for (int j = 0; j < dstWidth; j++)
			pDstImage[i][j] = tempdsthang[j];
	}

    free(tempAline);
    free(tempdstline);
    free(pTmpImage);
    free(tempAhang);
    free(tempdsthang);

    return 0;
}

double getThr(double **pDetCoef, int height,int width, int N)
{
/**
 * @param pDetCoef  待求取小波阈值的单通道图像的二级指针
 * @param height    信号的高
 * @param width     信号的宽
 * @param N	        原始信号的总数
 * @func            实现二维数组的小波阈值求取
 */
    double thr = 0.0;
	double sigma = 0.0;
    double *temp = (double *) malloc(height/2 * width/2 * sizeof(double));
    double sum = 0.0;
    int num = 0;

    // 对二维数组取绝对值
    for (int i=height/2; i<height; i++)
    {
        for (int j=width/2; j<width; j++)
        {
            temp[(i-height / 2) * width / 2 + j - width/2] = fabs(pDetCoef[i][j]);
            sum += fabs(pDetCoef[i][j]);
            num += 1;
        }
    }

    // 求取sigma
    sigma = (sum/num)/0.6745;
    thr = sigma * sqrt(2.0 * log(N));

    free(temp);

    return thr;
}

int Wthresh(double **pDstCoef, double thr, const int height, const int width, int model)
{
/**
 * @param pDetCoef  待小波阈值变换的单通道图像的二级指针
 * @param thr       小波阈值
 * @param height    图像的高
 * @param width     图像的宽
 * @param model	    1为硬阈值，0为软阈值
 * @func            实现二维数组的小波阈值变换
 */
    if (model)//硬阈值
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (((i>=height/2) || (j>=width/2)) && (fabs(pDstCoef[i][j]) < thr))
                    pDstCoef[i][j] = 0.0;
            }
        }
    }
    else//软阈值
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (((i>=height/2) || (j>=width/2)) && (fabs(pDstCoef[i][j]) < thr))
                {
                    pDstCoef[i][j] = 0.0;
                }
                else if (((i>=height/2) || (j>=width/2)))
                {
                    if (pDstCoef[i][j]<0)
                        pDstCoef[i][j] = thr - fabs(pDstCoef[i][j]);
                    else
                        pDstCoef[i][j] = fabs(pDstCoef[i][j]) - thr;
                }              
            }        
        }
    }
    return 0;
}

int main()
{
	//Setting YUV information
    FILE *input_fp;
    FILE *outputY_idw;
    int h = H, w = W, filterLen = FilterLen;
    const char *url = url_in;
    outputY_idw = fopen(url_out, "wb+");

    // 得到Y通道的数组指针。
	if ((input_fp = fopen(url, "rb")) == NULL)
	{
		printf("%s open error!\n", url);
	}
	else
	{
		printf("%s open.\n", url);
	}

    for (int t = 0; t <T; t++)
    {
        printf("Processing frame %d \n", t+1);
        unsigned char *pic =(unsigned char *) malloc(w * h * 3 / 2 * sizeof(unsigned char));

        double **input_Y_ui = (double**)malloc(sizeof(double*) * h);  
        for (int i = 0;i < h;i++)
        {
            input_Y_ui[i] = (double*)malloc(sizeof(double) * w);
        }

        fread(pic, sizeof(unsigned char), w * h * 3 / 2, input_fp);

        for (int i = 0; i<h; i++)
        {
            for (int j = 0; j<w; j++)
            {
                input_Y_ui[i][j] = (double)pic[i * w + j];
            }
        }

        // 小波一级分解
        int h1 = (h+filterLen-1)/2 * 2, w1 = (w+filterLen-1)/2 * 2;
        double **Y_dwt1 = (double**)malloc(sizeof(double*) * h1);
        for (int i = 0;i < h1;i++)
        {
            Y_dwt1[i] = (double*)malloc(sizeof(double) * w1);
        }

        DWT_2D(input_Y_ui, h, w, filterLen, Y_dwt1);

        //提取cA1
        double **cA1 = (double**)malloc(sizeof(double*) * h1/2);
        for (int i = 0;i < h1/2;i++)
        {
            cA1[i] = (double*)malloc(sizeof(double) * w1/2);
        }

        for (int i = 0;i< h1/2;i++)
        {
            for (int j=0;j< w1/2;j++)
            {
                cA1[i][j] = Y_dwt1[i][j];
            }
        }

        //二级分解
        int h2 = (h1/2+filterLen-1)/2 * 2, w2 = (w1/2+filterLen-1)/2 * 2;
        double **Y_dwt2 = (double**)malloc(sizeof(double*) * h2);
        for (int i = 0;i < h2;i++)
        {
            Y_dwt2[i] = (double*)malloc(sizeof(double) * w2);
        } 
        DWT_2D(cA1, h1/2, w1/2, 16, Y_dwt2);

        //提取CA2
        double **cA2 = (double**)malloc(sizeof(double*) * h2/2);
        for (int i = 0;i < h2/2;i++)
        {
            cA2[i] = (double*)malloc(sizeof(double) * w2/2);
        }

        for (int i = 0;i< h2/2;i++)
        {
            for (int j=0;j< w2/2;j++)
            {
                cA2[i][j] = Y_dwt2[i][j];
            }
        }
        
        //二层双边滤波
        double **cA2_tran = (double**)malloc(sizeof(double*) * h2/2);  
        for (int i = 0;i < h2/2;i++)
        {
            cA2_tran[i] = (double*)malloc(sizeof(double) * w2/2);
        }
        bilateral_filter(cA2, h2/2, w2/2, 4, 20, 1.8, cA2_tran);

        //cA赋值
        for (int i = 0;i<h2/2;i++)
        {
            for (int j=0;j<w2/2;j++)
            {
                Y_dwt2[i][j] = cA2_tran[i][j];
            }
        }

        //小波阈值去噪
        double thr2 = getThr(Y_dwt2, h2, w2, h*w);
        Wthresh(Y_dwt2, thr2, h2, w2, 1);

        //二级重构
        double **Y_idwt2 = (double**)malloc(sizeof(double*) * h1/2);
        for (int i = 0;i < h1/2;i++)
        {
            Y_idwt2[i] = (double*)malloc(sizeof(double) * w1/2);
        } 
        IDWT_2D(Y_dwt2, h1/2, w1/2, 16, Y_idwt2);

        //二层双边滤波
        double **cA1_tran = (double**)malloc(sizeof(double*) * h1/2);  
        for (int i = 0;i < h1/2;i++)
        {
            cA1_tran[i] = (double*)malloc(sizeof(double) * w1/2);
        }

        bilateral_filter(Y_idwt2, h1/2, w1/2, 4, 20, 1.8, cA1_tran);

        //cA赋值
        for (int i = 0;i<h1/2;i++)
        {
            for (int j=0;j<w1/2;j++)
            {
                Y_dwt1[i][j] = cA1_tran[i][j];
            }
        }

        //小波阈值去噪
        double thr1 = getThr(Y_dwt1, h1, w1, h*w);
        Wthresh(Y_dwt1, thr1, h1, w1, 1);  

        //一级重构
        double **Y_idwt1 = (double**)malloc(sizeof(double*) * h);
        for (int i = 0;i < h;i++)
        {
            Y_idwt1[i] = (double*)malloc(sizeof(double) * w);
        }
        IDWT_2D(Y_dwt1, h,w,16,Y_idwt1);

        // for (int i = 0;i< 20;i++)
        // {
        //     for (int j=0;j< 20;j++)
        //     {
        //         printf("%f ",Y_idwt1[i][j]);
        //     }
        //     printf("\n");
        // }

        // 图像double转为unsigned char
        unsigned char **output_Y = (unsigned char**)malloc(sizeof(unsigned char*) * h);  
        for (int i = 0;i < h;i++)
        {
            output_Y[i] = (unsigned char*)malloc(sizeof(unsigned char) * w);
        }

        for (int i=0; i<h; i++){
            for (int j=0; j<w; j++)
            {
                if (Y_idwt1[i][j]<0)
                    Y_idwt1[i][j] = 0;
                else if (Y_idwt1[i][j]>255)
                    Y_idwt1[i][j]=255;
                output_Y[i][j] = (unsigned char) Y_idwt1[i][j];
            }
        }

        // 将滤波的Y通道数据写回图像
        for (int i=0; i<h; i++){
            for (int j=0; j<w; j++){
                pic[i*w+j] = output_Y[i][j];
            }
        } 

        // 滤波图像输出文件
        fwrite(pic, sizeof(unsigned char), w * h * 3 / 2, outputY_idw);

        // 释放内存
        free(pic);

        for (int i = 0;i < h;i++)
            free(input_Y_ui[i]);
        free(input_Y_ui);

        for (int i = 0;i < h1;i++)
            free(Y_dwt1[i]);
        free(Y_dwt1);

        for (int i = 0;i < h1/2;i++)
            free(cA1[i]);
        free(cA1);

        for (int i = 0;i < h1/2;i++)
            free(cA1_tran[i]);
        free(cA1_tran);

        for (int i = 0;i < h2;i++)
            free(Y_dwt2[i]);
        free(Y_dwt2);

        for (int i = 0;i < h1/2;i++)
            free(Y_idwt2[i]);
        free(Y_idwt2);

        for (int i = 0;i < h;i++)
            free(Y_idwt1[i]);
        free(Y_idwt1);

        for (int i = 0;i < h2/2;i++)
            free(cA2[i]);
        free(cA2);

        for (int i = 0;i < h2/2;i++)
            free(cA2_tran[i]);
        free(cA2_tran);

        for (int i = 0;i < h;i++)
            free(output_Y[i]);
        free(output_Y);
    }

	return 0;
}