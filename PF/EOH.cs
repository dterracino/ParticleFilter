using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    class EOH:Feature
    {
        const int BIN_NUM = 4;     /* 梯度区间个数 */

        double[] targetModel;

        public override void Initialize(Image<Bgr, byte> image, int x0, int y0, int wx, int hy)
        {
            targetModel = new double[BIN_NUM];
            /* 计算模板积分图 */
            double[][] integralImages;
            CalcEOHIntegralImage(image.Convert<Bgr, byte>(), x0 - wx, y0 - hy, x0 + wx, y0 + hy, out integralImages);

            GetEm(integralImages, targetModel);
            for (int i = 0; i < BIN_NUM; ++i)
            {
                targetModel[i] = integralImages[i][wx * 2 + hy * 2 - 1];
                targetModel[i] /= (wx * 2 * hy * 2);
            }
        }

        public override void CalcWeights(Image<Bgr, byte> image, SpaceState[] states, double[] weights)
        {
            // 找到左上和右下坐标
            int minX = int.MaxValue, minY = int.MaxValue, maxX = int.MinValue, maxY = int.MinValue;
            for (int i = 0; i < states.Length; ++i)
            {
                if (states[i].xt - states[i].Hxt < minX) minX = states[i].xt - states[i].Hxt;
                if (states[i].xt + states[i].Hxt > maxX) maxX = states[i].xt + states[i].Hxt;
                if (states[i].yt - states[i].Hyt < minY) minY = states[i].yt + states[i].Hyt;
                if (states[i].yt + states[i].Hyt > maxY) maxY = states[i].yt + states[i].Hyt;
            }
            // 计算积分图
            double[][] integralImages;
            double[] particleModel = new double[BIN_NUM];
            CalcEOHIntegralImage(image, minX, minY, maxX, maxY, out integralImages);
            // 计算粒子模型
            for (int j = 0; j < states.Length; ++j)
            {
                GetEm(integralImages, particleModel, ref states[j]);
                // Bhattacharyya系数
                double rho = CalcBhattacharyya(targetModel, particleModel);
                weights[j] = CalcWeightedPi(rho);
            }
        }

        public override int ModelUpdate(Image<Bgr, byte> image, ref SpaceState estState)
        {
            double[] estHist;
            double rho, Pi_E;

            estHist = new double[BIN_NUM];

            /* (1)在估计值处计算目标直方图 */
            double[][] integralImages;
            double[] particleModel = new double[BIN_NUM];
            CalcEOHIntegralImage(image, estState.xt - estState.Hxt, estState.yt - estState.Hyt,
                estState.Hxt * 2, estState.Hyt * 2, out integralImages);
            /* (2)计算Bhattacharyya系数 */
            rho = CalcBhattacharyya(estHist, modelHist);
            /* (3)计算概率权重 */
            Pi_E = CalcWeightedPi(rho);

            if (Pi_E > Pi_Thres)
            {
                for (int i = 0; i < BIN_NUM; i++)
                {
                    modelHist[i] = ((1.0 - ALPHA_COEFFICIENT) * modelHist[i] + ALPHA_COEFFICIENT * estHist[i]);
                }
                return 1;
            }

            return -1;
        }

        private void GetEm(double[][] integralImage, double[] model, ref SpaceState state)
        {
            int 
            int leftTop = (state.xt - state.Hxt)(state.yt - state.Hyt) * ;
            int rightBottom = state.yt - state.Hyt;
            for (int i = 0; i < BIN_NUM; ++i)
            {
                if (leftTop < 0 || rightBottom >= integralImage.Length)
                    model[i] = 0;
                else
                    model[i] = (integralImage[i][rightBottom] - integralImage[i][leftTop]) / (state.Hxt * 2 + state.Hyt * 2);
            }
        }

        private void GetEm(double[][] integralImage, double[] model, int x0, int y0, int wx, int hy)
        {
            int leftTop = x0 - wx;
            int rightBottom = state.yt - state.Hyt;
            for (int i = 0; i < BIN_NUM; ++i)
            {
                if (leftTop < 0 || rightBottom >= integralImage.Length)
                    model[i] = 0;
                else
                    model[i] = (integralImage[i][rightBottom] - integralImage[i][leftTop]) / (state.Hxt * 2 + state.Hyt * 2);
            }
        }

        private void CalcIntegralImage(double[] inputMatrix, int width, int height, out double[] outputMatrix)
        {
            outputMatrix = new double[width * height];
            double[] columnSum = new double[width]; // sum of each column
            // calculate integral of the first line 
            for (int i = 0; i < width; i++)
            {
                columnSum[i] = inputMatrix[i];
                outputMatrix[i] = inputMatrix[i];
                if (i > 0)
                {
                    outputMatrix[i] += outputMatrix[i - 1];
                }
            }
            for (int i = 1; i < height; i++)
            {
                int offset = i * width;
                // first column of each line  
                columnSum[0] += inputMatrix[offset];
                outputMatrix[offset] = columnSum[0];
                // other columns   
                for (int j = 1; j < width; j++)
                {
                    columnSum[j] += inputMatrix[offset + j];
                    outputMatrix[offset + j] = outputMatrix[offset + j - 1] + columnSum[j];
                }
            }  
        }

        private void CalcEOHIntegralImage(Image<Bgr, byte> image, int x0, int y0, int x1, int y1, 
            out double[][] integralImages)
        {
            integralImages = null;
            /* 计算实际高宽和区域起始点 */
            int x_begin = x0;
            int y_begin = y0;
            if (x_begin < 0) x_begin = 0;
            if (y_begin < 0) y_begin = 0;
            int x_end = x1;
            int y_end = y1;
            if (x_end >= image.Width) x_end = image.Width - 1;
            if (y_end >= image.Height) y_end = image.Height - 1;
            int width = x_end - x_begin;
            int height = y_end - y_begin;
            if (width <= 0 || height <= 0) return;
            /* 取出实际需要计算的子图，并转为灰度图 */
            var subImage = image.GetSubRect(new Rectangle(x_begin, y_begin, width, height)).Convert<Gray, byte>();
            /* Convolution */
            float[,] ky = new float[,] { { -1, 0, 1 }, 
                                         { -2, 0, 2 }, 
                                         { -1, 0, 1 } };
            float[,] kx = new float[,] { { -1, -2, -1 }, 
                                         {  0,  0, 0 }, 
                                         {  1, 2, 1 } };
            var gx = subImage.Convolution(new ConvolutionKernelF(kx));
            var gy = subImage.Convolution(new ConvolutionKernelF(ky));

            double[][] binMatrix = new double[BIN_NUM][];
            for (int i = 0; i < BIN_NUM; ++i)
            {
                binMatrix[i] = new double[height * width];
            }

            double theta, magnitude;

            /* calculate theta and magnitude */
            for (int y = 0; y < height; ++y)
            {
                for (int x = 0; x < width; ++x)
                {
                    theta = Math.Atan2(gy[y, x].Intensity, gx[y, x].Intensity); /* range(-pi,+pi]*/
                    theta += Math.PI; /* trun to range(0, +2pi] */
                    int bin_num = ((int)(theta * BIN_NUM / Math.PI + 0.5) - 1);
                    magnitude = Math.Sqrt(gx[y, x].Intensity * gx[y, x].Intensity +
                                          gy[y, x].Intensity * gy[y, x].Intensity);

                    binMatrix[bin_num][y * width + x] = magnitude;
                }
            }

            /* calculate imtegral image */
            integralImages = new double[BIN_NUM][];
            for (int i = 0; i < BIN_NUM; ++i)
            {
                CalcIntegralImage(binMatrix[i], width, height, out integralImages[i]);
            }
        }
    }
}
