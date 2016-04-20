using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    class ColorDistribution:Feature
    {
        /* 常量 */
        const int R_BINS = 8;                               /* R分量的bin数目 */
        const int G_BINS = 8;                               /* G分量的bin数目 */
        const int B_BINS = 8;                               /* B分量的bin数目 */
        const int R_SHIFT = 5;                              /* R分量对应的位移 */
        const int G_SHIFT = 5;                              /* G分量对应的位移 */
        const int B_SHIFT = 5;                              /* B分量对应的位移 */
        const int BIN_NUM = R_BINS * G_BINS * B_BINS;       /* 总bin的数目 */
        const double Pi_Thres = 0.9;                        /* 权重阈值 */
        const double ALPHA_COEFFICIENT = 0.2;               /* 目标模型更新权重取0.1-0.3 */

        /* 目标模型 */
        double[] modelHist;

        public override void Initialize(Image<Bgr, byte> image, int x0, int y0, int wx, int hy)
        {
            modelHist = new double[BIN_NUM];
            /* calculate target color distribution */
            CalcColorHistogram(image, x0, y0, wx, hy, modelHist);
        }

        public override void CalcWeights(Image<Bgr, byte> image, SpaceState[] states, double[] weights)
        {
            Parallel.For(0, states.Length, (i) =>
            {
                double[] colorHist = new double[BIN_NUM];
                double rho = 0.0;
                /* (1) 计算彩色直方图分布 */
                CalcColorHistogram(image, states[i].xt, states[i].yt, states[i].Hxt, states[i].Hyt, colorHist);
                /* (2) Bhattacharyya系数 */
                rho = CalcBhattacharyya(colorHist, modelHist);
                /* (3) 根据计算得的Bhattacharyya系数计算各个权重值 */
                weights[i] = CalcWeightedPi(rho);
            });
        }

        public override int ModelUpdate(Image<Bgr, byte> image, ref SpaceState estState)
        {
            double[] estHist;
            double rho, Pi_E;

            estHist = new double[BIN_NUM];

            /* (1)在估计值处计算目标直方图 */
            CalcColorHistogram(image, estState.xt, estState.yt, estState.Hxt, estState.Hyt, estHist);
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

        private void CalcColorHistogram(Image<Bgr, byte> image, int x0, int y0, int wx, int hy, double[] colorHist)
        {
            /* 指定图像区域的左上角坐标 */
            int x_begin, y_begin;
            int y_end, x_end;
            double a2;
            int r, g, b;
            double f, r2, kr;
            int index;

            /* 直方图各个值赋0 */
            for (int i = 0; i < BIN_NUM; ++i)
                colorHist[i] = 0.0;

            /* 考虑特殊情况：x0, y0在图像外面，或者，Wx<=0, Hy<=0 */
            /* 此时强制令彩色直方图为0 */
            if ((x0 < 0) || (x0 >= image.Width) ||
                (y0 < 0) || (y0 >= image.Height) ||
                (wx <= 0) || (hy <= 0))
                return;

            /* 计算实际高宽和区域起始点 */
            x_begin = x0 - wx;
            y_begin = y0 - hy;
            if (x_begin < 0) x_begin = 0;
            if (y_begin < 0) y_begin = 0;
            x_end = x0 + wx;
            y_end = y0 + hy;
            if (x_end >= image.Width) x_end = image.Width - 1;
            if (y_end >= image.Height) y_end = image.Height - 1;

            a2 = wx * wx + hy * hy;     /* 计算核函数半径平方a^2 */
            f = 0.0;                    /* 归一化系数 */

            for (int y = y_begin; y <= y_end; y++)
            {
                for (int x = x_begin; x <= x_end; x++)
                {
                    // image[height, width]
                    r = ((byte)image[y, x].Red) >> R_SHIFT;     /* 计算直方图 */
                    g = ((byte)image[y, x].Green) >> G_SHIFT;   /*移位位数根据R、G、B条数 */
                    b = ((byte)image[y, x].Blue) >> B_SHIFT;
                    index = r * G_BINS * B_BINS + g * B_BINS + b;
                    r2 = ((y - y0) * (y - y0) + (x - x0) * (x - x0)) * 1.0 / a2; /* 计算半径平方r^2 */
                    kr = 1 - r2;   /* 核函数k(r) = 1-r^2, |r| < 1; 其他值 k(r) = 0 */
                    f += kr;
                    colorHist[index] += kr;  /* 计算核密度加权彩色直方图 */
                }
            }

            /* 归一化直方图 */
            for (int i = 0; i < BIN_NUM; i++)
                colorHist[i] /= f;
        }
    }
}
