#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
# define H 1080
# define W 1920
# define T 250
# define FilterLen 16
# define url_in "C:/Users/Cambricon/Desktop/yuvdataset/testdata_frame/HisiYUV_1920x1080_8bits_420sp_linear_128x_TR26_pavement.yuv"
# define url_out "C:/Users/Cambricon/Desktop/yuvdataset/output_128x_TR26_pavement_fixed.yuv"

// 乘以了pow(2,20)
int Lo_D[16] = {-123, 708, -410, -5106, 9170, 14660, -46229, -18213, 135001, 495, -297811, -16598, 613788, 708450, 328069, 57059};
int Hi_D[16] = {-57059, 328069, -708450, 613788, 16598, -297811, -495, 135001, 18213, -46229, -14660, 9170, 5106, -410, -708, -123};
int Lo_R[16] = {57059, 328069, 708450, 613788, -16598, -297811, 495, 135001, -18213, -46229, 14660, 9170, -5106, -410, 708, -123};
int Hi_R[16] = {-123, -708, -410, 5106, 9170, -14660, -46229, 18213, 135001, -495, -297811, 16598, 613788, -708450, 328069, -57059};


int bilateral_filter(int **image, int img_height, int img_width, int winsize, int sigmaColor, int sigmaSpace, int **image_tran, int ** LUT_space, int ** LUT_color)
{
/**
 * @param image	      待滤波图像数组的二级指针
 * @param img_height  待滤波图像的高
 * @param img_width	  待滤波图像的宽
 * @param winsize 	  滤波窗口的半径
 * @param sigmaColor  色彩滤波器超参,int[10,20]
 * @param sigmaSpace  空间滤波器超参,double[1,3]->定点化后int[16,48]
 * @param image_tran  双边滤波后图像数组的二级指针
 * @param LUT         指数函数查找表
 * @func              实现单通道图像的双边滤波功能
 */

    // 初始化超参数,乘以pow(2, 14)
    int space_coeff = - 8192 / (sigmaSpace * sigmaSpace);//6位

    // 计算直径mask，初始化空间滤波器的存储数组weight_space，坐标数组weight_space_row、weight_space_col
    int mask = (2*winsize+1)*(2*winsize+1);
    int *weight_space = (int *)malloc(mask*sizeof(int));
    int *p = weight_space;
    int *weight_space_row = (int *)malloc(mask*sizeof(int));
    int *pr = weight_space_row;
    int *weight_space_col = (int *)malloc(mask*sizeof(int));
    int *pc = weight_space_col;

    // 计算高斯核的参数
    for (int i=-winsize; i<winsize+1; i++){
        for (int j=-winsize; j<winsize+1; j++){
            int r_square = i*i + j*j;
            // *p = exp(r_square * space_coeff);r_square为[0, 51),space_coeff为[-32, -3],乘以pow(2,14)
            *p = LUT_space[r_square][space_coeff+32];
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
            long long int value = 0;
            long long int weight = 0;
            for (int i=0; i<mask; i++){
                int m = row + weight_space_row[i];
                int n = col + weight_space_col[i];
                int val = 0;
                long long int w = 0;
                int temp1;
                int temp2;
                if ((m<0) || (n<0) || (m >= img_height) || (n>=img_width))
                    val=0;
                else
                    val = image[m][n];
                // w = weight_space[i] * exp(-0.5 / (sigmaColor * sigmaColor) * pow(abs(val - image[row][col]), 2.0));
                if (abs(val - image[row][col]) >= 130)
                    w = 0;
                else
                    w = weight_space[i] * LUT_color[sigmaColor-10][abs(val - image[row][col])];//29位,为了防止溢出用longlong
                value += val * w;
                weight += w;
            }
            image_tran[row][col] = (int)(value/weight);
        }
    }
    free(weight_space);
    free(weight_space_row);
    free(weight_space_col);
    return 0;
}

int DWT(int *pSrcData,int srcLen,int filterLen,int *pDstCeof)
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
	int tmp = 0;
	for (int i = 0; i < exLen; i++)
	{
		pDstCeof[i] = 0;
		pDstCeof[i + exLen] = 0;
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

int DWT_2D(int **pSrcImage, int height, int width, int filterLen, int **pDstCeof)
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

	int *tempImage = (int *) malloc(exwidth*height * sizeof(int)); 

	// 对每一行进行行变换
	int *tempAhang = (int *) malloc(width * sizeof(int)); // 临时存放每一行的未处理数据
	int *tempExhang = (int *) malloc(exwidth * sizeof(int));// 临时存放每一行的处理数据

	for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            tempAhang[j] = pSrcImage[i][j];//提取每一行的数据
        }

		DWT(tempAhang, width, filterLen, tempExhang);

		// 所有数右移20位
		for (int j = 0; j < exwidth; j++)
        {
            tempExhang[j] = tempExhang[j]>>20;
			tempImage[i*exwidth + j] = tempExhang[j];
        }
    }

	// 对每一列进行列变换
	int *tempAlie = (int *) malloc(height * sizeof(int)); // 临时存放每一列的转置数据
	int *tempEXlie = (int *) malloc(exheight * sizeof(int)); // 临时存放每一列的处理数据

	for (int i = 0; i < exwidth; i++)
	{
		// 列转置
		for (int j = 0; j < height; j++)
			tempAlie[j] = tempImage[j*exwidth + i];//提取每一列数据

		//执行变换
		DWT(tempAlie, height, filterLen, tempEXlie);

		// 所有数右移20位，反转置
		for (int j = 0; j < exheight; j++)
        {
            tempEXlie[j] = tempEXlie[j]>>20;
			pDstCeof[j][i] = tempEXlie[j];
        }
	}

    free(tempImage);
    free(tempAhang);
    free(tempExhang);
    free(tempAlie);
    free(tempEXlie);

    return 0;
}

int IDWT(int *pSrcCoef, int dstLen, int filterLen, int *pDstData)
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
		pDstData[i] = 0;
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

int IDWT_2D(int **pSrcCeof, int dstHeight, int dstWidth, int filterLen, int **pDstImage)
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

    int *tempAline = (int *) malloc(srcHeight * sizeof(int)); // 临时存放每一列的数据
	int *tempdstline = (int *) malloc(dstHeight * sizeof(int)); // 临时存放每一列的重构结果

	int *pTmpImage = (int *) malloc(srcWidth*dstHeight * sizeof(int));

	// 列重构
	for (int i = 0; i < srcWidth; i++)//每一列
    {
		// 列转置
		for (int j = 0; j<srcHeight; j++)
			tempAline[j] = pSrcCeof[j][i];//提取每一列

		IDWT(tempAline, dstHeight, filterLen, tempdstline);

		// 反转置
		for (int j = 0; j < dstHeight; j++)
        {
            tempdstline[j] = tempdstline[j]>>20;
			pTmpImage[j*srcWidth + i] = tempdstline[j];
        }

    }

    // 行重构
	int *tempAhang = (int *) malloc(srcWidth * sizeof(int));// 临时存放每一行的数据
	int *tempdsthang = (int *) malloc(dstWidth * sizeof(int));// 临时存放每一行的重构数据
	for (int i = 0; i < dstHeight; i++)
	{
		for (int j = 0; j < srcWidth; j++)
			tempAhang[j] = pTmpImage[i*srcWidth + j];//提取每一行的数据

		IDWT(tempAhang, dstWidth, filterLen, tempdsthang);

		for (int j = 0; j < dstWidth; j++)
        {
            tempdsthang[j] = tempdsthang[j]>>20;
			pDstImage[i][j] = tempdsthang[j];
        }
	}

    free(tempAline);
    free(tempdstline);
    free(pTmpImage);
    free(tempAhang);
    free(tempdsthang);

    return 0;
}

int getThr(int **pDetCoef, int height,int width, int N)
{
/**
 * @param pDetCoef  待求取小波阈值的单通道图像的二级指针
 * @param height    信号的高
 * @param width     信号的宽
 * @param N	        原始信号的总数
 * @func            实现二维数组的小波阈值求取
 */
    int thr = 0;
	int sigma = 0;
    int *temp = (int *) malloc(height/2 * width/2 * sizeof(int));
    int sum = 0.0;
    int num = 0;

    // 对二维数组取绝对值
    for (int i=height/2; i<height; i++)
    {
        for (int j=width/2; j<width; j++)
        {
            temp[(i-height / 2) * width / 2 + j - width/2] = abs(pDetCoef[i][j]);
            sum += fabs(pDetCoef[i][j]);
            num += 1;
        }
    }

    // 求取与之，131010为系数乘pow(2,14)的结果
    thr = (sum/num) * 131010;
    thr = thr >> 14;

    free(temp);

    return thr;
}

int Wthresh(int **pDstCoef, int thr, const int height, const int width, int model)
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
                if (((i>=height/2) || (j>=width/2)) && (abs(pDstCoef[i][j]) < thr))
                    pDstCoef[i][j] = 0;
            }
        }
    }
    else//软阈值
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (((i>=height/2) || (j>=width/2)) && (abs(pDstCoef[i][j]) < thr))
                {
                    pDstCoef[i][j] = 0;
                }
                else if (((i>=height/2) && (j>=width/2)))
                {
                    if (pDstCoef[i][j]<0)
                        pDstCoef[i][j] = thr - abs(pDstCoef[i][j]);
                    else
                        pDstCoef[i][j] = abs(pDstCoef[i][j]) - thr;
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

        int **input_Y_ui = (int**)malloc(sizeof(int*) * h);  
        for (int i = 0;i < h;i++)
        {
            input_Y_ui[i] = (int*)malloc(sizeof(int) * w);
        }

        fread(pic, sizeof(unsigned char), w * h * 3 / 2, input_fp);

        for (int i = 0; i<h; i++)
        {
            for (int j = 0; j<w; j++)
            {
                input_Y_ui[i][j] = (int)pic[i * w + j];
            }
        }

        //读取查找表
        int LUT_color_row = 0;
        int LUT_color_col = 0;
        FILE *fp_color = NULL;
        if((fp_color = fopen("C:/Users/Cambricon/Desktop/yuvdataset/LUT_color_exp_15bit.txt", "r")) == NULL)
        {
            printf("fopen failed\n");
        }
        fscanf(fp_color, "%d %d", &LUT_color_row, &LUT_color_col);

        int **LUT_color = (int**)malloc(sizeof(int*) * LUT_color_row);  
        for (int i = 0;i < LUT_color_row;i++)
        {
            LUT_color[i] = (int*)malloc(sizeof(int) * LUT_color_col);
        }

        for (int i = 0; i < LUT_color_row; i++)
        {
            for (int j = 0; j < LUT_color_col; j++)
            {
                fscanf(fp_color, "%d", &LUT_color[i][j]);
            }
        }
    
        int LUT_space_row = 0;
        int LUT_space_col = 0;
        FILE *fp_space = NULL;
        if((fp_space = fopen("C:/Users/Cambricon/Desktop/yuvdataset/LUT_space_exp_14bit.txt", "r")) == NULL)
        {
            printf("fopen failed\n");
        }
        fscanf(fp_space, "%d %d", &LUT_space_row, &LUT_space_col);

        int **LUT_space = (int**)malloc(sizeof(int*) * LUT_space_row);  
        for (int i = 0;i < LUT_space_row;i++)
        {
            LUT_space[i] = (int*)malloc(sizeof(int) * LUT_space_col);
        }

        for (int i = 0; i < LUT_space_row; i++)
        {
            for (int j = 0; j < LUT_space_col; j++)
            {
                fscanf(fp_space, "%d", &LUT_space[i][j]);
            }
        }

        // 小波一级分解
        int h1 = (h+filterLen-1)/2 * 2, w1 = (w+filterLen-1)/2 * 2;

        int **Y_dwt1 = (int**)malloc(sizeof(int*) * h1);
        for (int i = 0;i < h1;i++)
        {
            Y_dwt1[i] = (int*)malloc(sizeof(int) * w1);
        }

        DWT_2D(input_Y_ui, h, w, filterLen, Y_dwt1);

        //提取cA1
        int **cA1 = (int**)malloc(sizeof(int*) * h1/2);
        for (int i = 0;i < h1/2;i++)
        {
            cA1[i] = (int*)malloc(sizeof(int) * w1/2);
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
        int **Y_dwt2 = (int**)malloc(sizeof(int*) * h2);
        for (int i = 0;i < h2;i++)
        {
            Y_dwt2[i] = (int*)malloc(sizeof(int) * w2);
        } 
        DWT_2D(cA1, h1/2, w1/2, 16, Y_dwt2);

        //提取CA2
        int **cA2 = (int**)malloc(sizeof(int*) * h2/2);
        for (int i = 0;i < h2/2;i++)
        {
            cA2[i] = (int*)malloc(sizeof(int) * w2/2);
        }

        for (int i = 0;i< h2/2;i++)
        {
            for (int j=0;j< w2/2;j++)
            {
                cA2[i][j] = Y_dwt2[i][j];
            }
        }
        
        //二层双边滤波
        int **cA2_tran = (int**)malloc(sizeof(int*) * h2/2);  
        for (int i = 0;i < h2/2;i++)
        {
            cA2_tran[i] = (int*)malloc(sizeof(int) * w2/2);
        }
        // 高斯核的参数为乘以pow(2,4)后
        bilateral_filter(cA2, h2/2, w2/2, 4, 20, 28, cA2_tran, LUT_space, LUT_color);

        //cA赋值
        for (int i = 0;i<h2/2;i++)
        {
            for (int j=0;j<w2/2;j++)
            {
                Y_dwt2[i][j] = cA2_tran[i][j];
            }
        }

        //小波阈值去噪
        int thr2 = getThr(Y_dwt2, h2, w2, h*w);
        Wthresh(Y_dwt2, thr2, h2, w2, 1);

        //二级重构
        int **Y_idwt2 = (int**)malloc(sizeof(int*) * h1/2);
        for (int i = 0;i < h1/2;i++)
        {
            Y_idwt2[i] = (int*)malloc(sizeof(int) * w1/2);
        } 
        IDWT_2D(Y_dwt2, h1/2, w1/2, 16, Y_idwt2);

        //二层双边滤波
        int **cA1_tran = (int**)malloc(sizeof(int*) * h1/2);  
        for (int i = 0;i < h1/2;i++)
        {
            cA1_tran[i] = (int*)malloc(sizeof(int) * w1/2);
        }

        bilateral_filter(Y_idwt2, h1/2, w1/2, 4, 20, 28, cA1_tran, LUT_space, LUT_color);

        //cA赋值
        for (int i = 0;i<h1/2;i++)
        {
            for (int j=0;j<w1/2;j++)
            {
                Y_dwt1[i][j] = cA1_tran[i][j];
            }
        }

        //小波阈值去噪
        int thr1 = getThr(Y_dwt1, h1, w1, h*w);
        Wthresh(Y_dwt1, thr1, h1, w1, 1);  

        //一级重构
        int **Y_idwt1 = (int**)malloc(sizeof(int*) * h);
        for (int i = 0;i < h;i++)
        {
            Y_idwt1[i] = (int*)malloc(sizeof(int) * w);
        }
        IDWT_2D(Y_dwt1, h,w,16,Y_idwt1);

        // for (int i = 0;i< 20;i++)
        // {
        //     for (int j=0;j< 20;j++)
        //     {
        //         printf("%d ",Y_idwt1[i][j]);
        //     }
        //     printf("\n");
        // }

        // 图像int转为unsigned char
        unsigned char **output_Y = (unsigned char**)malloc(sizeof(unsigned char*) * h);  
        for (int i = 0;i < h;i++)
        {
            output_Y[i] = (unsigned char*)malloc(sizeof(unsigned char) * w);
        }

        for (int i=0; i<h; i++){
            for (int j=0; j<w; j++)
            {
                if (Y_idwt1[i][j]<0)
                {
                    Y_idwt1[i][j] = 0;
                }
                else if (Y_idwt1[i][j]>255)
                {
                    Y_idwt1[i][j]=255;
                }
                output_Y[i][j] = (unsigned char)Y_idwt1[i][j];
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

        for (int i = 0;i < LUT_space_row;i++)
            free(LUT_space[i]);
        free(LUT_space);

        for (int i = 0;i < LUT_color_row;i++)
            free(LUT_color[i]);
        free(LUT_color);
    }

	return 0;
}