﻿using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    struct SpaceState                   /* 状态空间变量 */
    {
        public int xt;                  /* x坐标位置 */
        public int yt;                  /* x坐标位置 */
        public double v_xt;             /* x方向运动速度 */
        public double v_yt;             /* y方向运动速度 */
        public int Hxt;                 /* x方向半窗宽 */
        public int Hyt;                 /* y方向半窗宽 */
        public double at_dot;           /* 尺度变换速度 */
    }

    class ParticleFilter
    {
        /* 常量 */
        const double DELTA_T = 0.05;                /* 帧频，可以为30，25，15，10等 */
        const int POSITION_DISTURB = 15;            /* 位置扰动幅度   */
        const double VELOCITY_DISTURB = 40.0;       /* 速度扰动幅值   */
        const double SCALE_DISTURB = 0.0;           /* 窗宽高扰动幅度 */
        const double SCALE_CHANGE_D = 0.001;        /* 尺度变换速度扰动幅度 */
        const double SIGMA = 0.1;                   /* sigma = 0.1 */
        const double SIGMA2 = SIGMA * SIGMA * 2;    /* 2*sigma^2 */
        const double Pi_Thres = 0.9;                /* 权重阈值 */
        const double ALPHA_COEFFICIENT = 0.2;       /* 目标模型更新权重取0.1-0.3 */
        const double RECIP_SIGMA = 3.9894228040;    /* 1/(sqrt(2*pi)*sigma), 这里sigma = 0.1*/

        /* 粒子滤波相关变量 */
        int particleNum;            /* 粒子数 */
        SpaceState[] states;        /* 粒子状态 */
        double[] weights;           /* 例子权重 */
        Util randUtil;              /* 系统随机函数 */
        Feature feature;

        public ParticleFilter(int particleNumbers)
        {
            particleNum = particleNumbers;
            randUtil = new Util(0L);
            /* 方法1.使用color destribution */
            feature = new ColorDistribution();
            /* 方法2.使用EOH */
            //feature = new EOH();
        }

        public int Initialize(int x0, int y0, int wx, int hy, Image<Bgr, byte> image)
        {
            double []rn = new double[7];

            states = new SpaceState[particleNum];   /* 申请状态数组的空间 */
            weights = new double[particleNum];       /* 申请粒子权重数组的空间 */

            feature.Initialize(image, x0, y0, wx, hy);
            
            /* 初始化粒子状态(以(x0,y0,1,1,Wx,Hy,0.1)为中心呈N(0,0.6)正态分布) */
            states[0].xt = x0;  
            states[0].yt = y0;  
            states[0].v_xt = 0.0; // 1.0  
            states[0].v_yt = 0.0; // 1.0  
            states[0].Hxt = wx;  
            states[0].Hyt = hy;  
            states[0].at_dot = 0.0; // 0.1
            weights[0] = (1.0 / particleNum); /* 0.9; */
            for (int i = 1; i < particleNum; i++)  
            {  
                for (int j = 0; j < 7; j++ )
                    rn[j] = randUtil.RandGaussian(0, 0.6); /* 产生7个随机高斯分布的数 */  
                states[i].xt = (int)(states[0].xt + rn[0] * wx);  
                states[i].yt = (int)(states[0].yt + rn[1] * hy);  
                states[i].v_xt = (states[0].v_xt + rn[2] * VELOCITY_DISTURB);  
                states[i].v_yt = (states[0].v_yt + rn[3] * VELOCITY_DISTURB);  
                states[i].Hxt = (int)(states[0].Hxt + rn[4] * SCALE_DISTURB);  
                states[i].Hyt = (int)(states[0].Hyt + rn[5] * SCALE_DISTURB);  
                states[i].at_dot = (states[0].at_dot + rn[6] * SCALE_CHANGE_D);  
                /* 权重统一为1/N，让每个粒子有相等的机会 */
                weights[i] = weights[0];
            }

            return 1;  
        }

        public int ColorParticleTracking(Image<Bgr, byte> image, out int xc, out int yc,
            out int wx_h, out int hy_h, out double max_weight)
        {
            SpaceState estState; /* 估计状态 */
            /* 选择：选择样本，并进行重采样 */
            Select();
            /* 传播：采样状态方程，对状态变量进行预测 */
            Propagate(image);
            /* 观测：对状态量进行更新 */
            Observe(image);
            /* 估计：对状态量进行估计，提取位置量 */
            Estimation(out estState);
            /* 设置各项返回值 */
            xc = estState.xt;
            yc = estState.yt;
            wx_h = estState.Hxt;
            hy_h = estState.Hyt;
            /* 模型更新 */
            ModelUpdate(ref estState, Pi_Thres, image);
            /* 计算最大权重值 */
            max_weight = weights.Max();
            DrawPriticles(image);
            /* 进行合法性检验，不合法返回-1 */
            if (xc < 0 || yc < 0 || xc - wx_h < 0 || yc - hy_h < 0 
                ||xc + wx_h >= image.Width || yc + hy_h >= image.Height  
                || wx_h <= 0 || hy_h <= 0) 
                return (-1);
            else
                return (1);     
        }

        private void DrawPriticles(Image<Bgr, byte> image)
        {
            foreach (var particle in states)
            {
                image.Draw(new CircleF(new PointF(particle.xt, particle.yt), 1), new Bgr(Color.Green), -1);
            }
        } 

        private int ModelUpdate(ref SpaceState estState, double PiT, Image<Bgr, byte> image)  
        {
            return feature.ModelUpdate(image, ref estState);
        }

        private void Estimation(out SpaceState estimateState)  
        {  
            double at_dot = 0.0, 
                   Hxt = 0.0,
                   Hyt = 0.0,
                   v_xt = 0.0,
                   v_yt = 0.0,
                   xt = 0.0,
                   yt = 0.0,
                   weight_sum = 0.0;  
  
            for (int i = 0; i < particleNum; i++) /* 加权平均 */  
            {  
                at_dot += states[i].at_dot * weights[i];  
                Hxt += states[i].Hxt * weights[i];  
                Hyt += states[i].Hyt * weights[i];  
                v_xt += states[i].v_xt * weights[i];  
                v_yt += states[i].v_yt * weights[i];  
                xt += states[i].xt * weights[i];  
                yt += states[i].yt * weights[i];
                weight_sum += weights[i];  
            }

            /* 求平均 */  
            estimateState.at_dot = at_dot / weight_sum;
            estimateState.Hxt = (int)(Hxt / weight_sum + 0.5);
            estimateState.Hyt = (int)(Hyt / weight_sum + 0.5);
            estimateState.v_xt = v_xt / weight_sum;
            estimateState.v_yt = v_yt / weight_sum;
            estimateState.xt = (int)(xt / weight_sum + 0.5);
            estimateState.yt = (int)(yt / weight_sum + 0.5);  
  
            return;  
        } 

        private void Observe(Image<Bgr, byte> image)
        {
            feature.CalcWeights(image, states, weights);
            /* normalize weights */
            var sum = weights.Sum();
            for (int i = 0; i < particleNum; ++i)
                weights[i] /= sum;

            return;
        }

        private void Propagate(Image<Bgr, byte> image)
        {
            double[] rn = new double[7];  
  
            /* 对每一个状态向量state[i](共N个)进行更新 */  
            for (int i = 0; i < particleNum; i++ )  /* 加入均值为0的随机高斯噪声 */  
            {
                for (int j = 0; j < 7; j++)
                    rn[j] = randUtil.RandGaussian(0, 0.6); /* 产生7个随机高斯分布的数 */
                states[i].xt = (int)(states[i].xt + states[i].v_xt * DELTA_T + rn[0] * states[i].Hxt + 0.5);
                states[i].yt = (int)(states[i].yt + states[i].v_yt * DELTA_T + rn[1] * states[i].Hyt + 0.5);
                states[i].v_xt = (states[i].v_xt + rn[2] * VELOCITY_DISTURB);
                states[i].v_yt = (states[i].v_yt + rn[3] * VELOCITY_DISTURB);
                states[i].Hxt = (int)(states[i].Hxt + states[i].Hxt * states[i].at_dot + rn[4] * SCALE_DISTURB + 0.5);
                states[i].Hyt = (int)(states[i].Hyt + states[i].Hyt * states[i].at_dot + rn[5] * SCALE_DISTURB + 0.5);
                states[i].at_dot = (states[i].at_dot + rn[6] * SCALE_CHANGE_D);
            }
            return;  
        }

        private void Select()
        {
            SpaceState[] tmpState = new SpaceState[particleNum];
            int[] reSamlpeIndex = new int[particleNum];
            ImportanceSampling(reSamlpeIndex);
            for (int i = 0; i < particleNum; ++i)
                tmpState[i] = states[reSamlpeIndex[i]]; //temState为临时变量,其中state[i]用state[reSamlpeIndex[i]]来代替
            for (int i = 0; i < particleNum; ++i)
                states[i] = tmpState[i];

            return;
        }

        private void ImportanceSampling(int[] reSampleIndex)
        {
            // double rnum;
            double[] cumulateWeight;
            int index;

            cumulateWeight = new double[particleNum + 1];    /* 申请累计权重数组内存，大小为N+1 */  
            NormalizeCumulatedWeight(cumulateWeight);       /* 计算累计权重 */

            // 改为并行的版本，提高效率
            Parallel.For(0, reSampleIndex.Length, (i) => 
            {
                double rnum;
                rnum = randUtil.Rand0_1();       /* 随机产生一个[0,1]间均匀分布的数 */
                index = BinearySearch(rnum, cumulateWeight); /* 搜索<=rnum的最小索引j */
                if (index == particleNum)
                    index--;
                reSampleIndex[i] = index;   /* 放入重采样索引数组 */    
            });

            return;   
        }

        private int BinearySearch(double v, double[] cumulateWeight)
        {
            int l = 0, 
                r = particleNum - 1, 
                m = 0;

            while (r >= l)
            {
                m = (l + r) / 2;
                if (v >= cumulateWeight[m] && v < cumulateWeight[m + 1]) 
                    return (m);
                if (v < cumulateWeight[m]) 
                    r = m - 1;
                else 
                    l = m + 1;
            }

            return (0);  
        }

        private void NormalizeCumulatedWeight(double[] cumulateWeight)
        {
            for (int i = 0; i < particleNum + 1; i++)
                cumulateWeight[i] = 0;
            for (int i = 0; i < particleNum; i++)
                cumulateWeight[i + 1] = cumulateWeight[i] + weights[i];
            for (int i = 0; i < particleNum + 1; i++)
                cumulateWeight[i] = cumulateWeight[i] / cumulateWeight[particleNum];

            return;  
        }
    }
}
