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
        const int BIN_NUM = 8;                              /* 梯度区间个数 */
        const double Pi_Thres = 0.9;                        /* 权重阈值 */
        const double ALPHA_COEFFICIENT = 0.2;               /* 目标模型更新权重取0.1-0.3 */

        double[] targetModel;
        Image<Gray, double>[] integralImages;
        Image<Gray, double>[] eohImages;
        Rectangle ROI;

        public override void Initialize(Image<Bgr, byte> image, int x0, int y0, int wx, int hy)
        {
            this.ROI = new Rectangle(x0 - wx, y0 - hy, wx * 2, hy * 2);
            /* 初始化eog图和积分图 */
            eohImages = new Image<Gray, double>[BIN_NUM];
            integralImages = new Image<Gray, double>[BIN_NUM];
            /* 初始化目标模型 */
            targetModel = new double[BIN_NUM];
            /* 计算积分图 */
            CalcEOHIntegralImage(image);
            CalcModel(x0, y0, wx, hy, targetModel);
        }

        public override void CalcWeights(Image<Bgr, byte> image, SpaceState[] states, double[] weights)
        {
            /* 计算上下左右四个点 */
            int minX = int.MaxValue, 
                minY = int.MaxValue, 
                maxX = int.MinValue, 
                maxY = int.MinValue;
            int value;
            foreach (var state in states)
            {
                if ((value = state.xt - state.Hxt) < minX) minX = value;
                if ((value = state.yt - state.Hyt) < minY) minY = value;
                if ((value = state.xt + state.Hxt) > maxX) maxX = value;
                if ((value = state.yt + state.Hyt) > maxY) maxY = value;
            }
            this.ROI = new Rectangle(minX, minY, maxX - minX, maxY - minY);
            /* 计算该区域积分图 */
            CalcEOHIntegralImage(image);
            for (int i = 0; i < states.Length; ++i)
            {
                double[] particleModel = new double[BIN_NUM];
                double rho;
                /* 计算粒子模型 */
                CalcModel(states[i].xt, states[i].yt, states[i].Hxt, states[i].Hyt, particleModel);
                /* 计算距离 */
                rho = CalcBhattacharyya(particleModel, targetModel);
                /* 计算权重 */
                weights[i] = CalcWeightedPi(rho);
            }
        }

        public override int ModelUpdate(Image<Bgr, byte> image, ref SpaceState estState)
        {
            double[] estModel;
            double rho, Pi_E;

            estModel = new double[BIN_NUM];

            ///* 计算该区域积分图 */
            //CalcEOHIntegralImage(image, new Rectangle(estState.xt - estState.Hxt, estState.yt - estState.Hyt,
            //                                          estState.Hxt * 2, estState.Hyt * 2));

            /* 计算粒子模型 */
            CalcModel(estState.xt, estState.yt, estState.Hxt, estState.Hyt, estModel);
            /* 计算距离 */
            rho = CalcBhattacharyya(estModel, targetModel);
            /* 计算权重 */
            Pi_E = CalcWeightedPi(rho);

            if (Pi_E > Pi_Thres)
            {
                for (int i = 0; i < BIN_NUM; i++)
                {
                    targetModel[i] = ((1.0 - ALPHA_COEFFICIENT) * targetModel[i] + ALPHA_COEFFICIENT * estModel[i]);
                }
                return 1;
            }

            return -1;
        }

        private void CalcModel(int x0, int y0, int wx, int hy, double[] model)
        {
            int x = x0 + wx - ROI.X,
                y = y0 + hy - ROI.Y,
                x_1 = x0 - wx - ROI.X,
                y_1 = y0 - hy - ROI.Y;
            double sum = 0.0;
            /* 初始化 */
            for (int i = 0; i < model.Length; ++i)
                model[i] = 0.0;
            /* 验证有效性 */
            if (x >= integralImages[0].Width || y >= integralImages[0].Height || x_1 < 0 || y_1 < 0)
                return;

            for (int i = 0; i < BIN_NUM; ++i)
            {
                /* four table-ups */
                var integral = integralImages[i];
                var em = integral[y, x].Intensity - integral[y_1, x].Intensity - integral[y, x_1].Intensity +
                    integral[y_1, x_1].Intensity;
                em /= (wx << 2) * (hy << 2);
                model[i] = em;
                sum += em;
            }
            /* 归一化 */
            for (int i = 0; i < BIN_NUM; ++i)
                model[i] /= sum;
        }

        private void CalcEOHIntegralImage(Image<Bgr, byte> image)
        {
            double theta, magnitude;
            /* 验证ROI的有效性 */
            if (ROI.X < 0) ROI.X = 0;
            if (ROI.Y < 0) ROI.Y = 0;
            if (ROI.X + ROI.Width > image.Width) ROI.Width = image.Width - ROI.X;
            if (ROI.Y + ROI.Height > image.Height) ROI.Height = image.Height - ROI.Y;
            /* 初始化 */
            for (int i = 0; i < BIN_NUM; ++i)
                eohImages[i] = new Image<Gray, double>(ROI.Width, ROI.Height);
            /* 转为灰度图 */
            CvInvoke.cvSetImageROI(image, ROI);
            var tmp = image.Convert<Gray, byte>();
            CvInvoke.cvResetImageROI(image);
            /* 卷积，使用Sobel算子计算梯度 */
            var gx = tmp.Sobel(1, 0, 3);
            var gy = tmp.Sobel(0, 1, 3);
            /* 计算梯度方向和大小 */
            for (int y = 0; y < ROI.Height; ++y)
            {
                for (int x = 0; x < ROI.Width; ++x)
                {
                    theta = Math.Atan2(gy[y, x].Intensity, gx[y, x].Intensity); /* range(-pi,+pi] */
                    theta += Math.PI; /* trun to range(0, +2pi] */
                    int bin_num = (int)Math.Floor(theta * BIN_NUM / (Math.PI * 2));
                    if (bin_num == BIN_NUM) bin_num--;
                    magnitude = Math.Sqrt(gx[y, x].Intensity * gx[y, x].Intensity +
                                          gy[y, x].Intensity * gy[y, x].Intensity);
                    eohImages[bin_num].Data[y, x, 0] = magnitude;
                }
            }
            /* 计算积分图 */
            for (int i = 0; i < BIN_NUM; ++i)
                integralImages[i] = eohImages[i].Integral();
        }
    }
}
