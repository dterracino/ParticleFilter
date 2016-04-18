using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    class Util
    {
        /* 
         * 采用Park and Miller方法产生[0,1]之间均匀分布的伪随机数 
         * 算法详细描述见： 
         * [1] NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING. 
         * Cambridge University Press. 1992. pp.278-279. 
         * [2] Park, S.K., and Miller, K.W. 1988, Communications of the ACM,  
         * vol. 31, pp. 1192–1201. */  
        const long IA = 16807;
        const long IM = 2147483647;
        const double AM = (1.0/IM);
        const long IQ = 127773;
        const long IR = 2836;
        const long MASK = 123459876;

        long ran_seed = 0L;         /* 随机数种子，为全局变量，设置缺省值 */

        public Util(long seed)
        {
            this.SetSeed(seed);
        }

        public long SetSeed(long seed)
        {
            if (seed != 0) /* 如果传入的参数setvalue!=0，设置该数为种子 */
                ran_seed = seed;
            else                 /* 否则，利用系统时间为种子数 */
            {
                ran_seed = Environment.TickCount;
            }
            return (ran_seed);
        }

        public double RandGaussian(double u, double sigma)
        {
            double x1, x2, v1 = 0.0, v2 = 0.0;
            double s = 100.0;
            double y;

            /* 
            使用筛选法产生正态分布N(0,1)的随机数(Box-Mulles方法) 
            1. 产生[0,1]上均匀随机变量X1,X2 
            2. 计算V1=2*X1-1,V2=2*X2-1,s=V1^2+V2^2 
            3. 若s<=1,转向步骤4，否则转1 
            4. 计算A=(-2ln(s)/s)^(1/2),y1=V1*A, y2=V2*A 
            y1,y2为N(0,1)随机变量 
            */
            while (s > 1.0)
            {
                x1 = Rand0_1();
                x2 = Rand0_1();
                v1 = 2 * x1 - 1;
                v2 = 2 * x2 - 1;
                s = v1 * v1 + v2 * v2;
            }
            y = (double)(Math.Sqrt(-2.0 * Math.Log(s) / s) * v1);
            /* 
            根据公式 
            z = sigma * y + u 
            将y变量转换成N(u,sigma)分布 
            */
            return (sigma * y + u);  
        }

        //public double[] Rand0_1Array(int count)
        //{
        //    double[] retArray = new double[count];
        //    for (int i = 0; i < count; ++i)
        //    {
        //        retArray[i] = Rand0_1();
        //    }
        //    return retArray;
        //}

        /* 采用Park and Miller方法产生[0,1]之间均匀分布的伪随机数 
         * 算法详细描述见： 
         * [1] NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING. 
         * Cambridge University Press. 1992. pp.278-279. 
         * [2] Park, S.K., and Miller, K.W. 1988, Communications of the ACM,  
         * vol. 31, pp. 1192–1201. */
        public double Rand0_1()
        {
            long k;
            double ans;

            /* *idum ^= MASK;*/
            /* XORing with MASK allows use of zero and other */
            k = ran_seed / IQ;                              /* simple bit patterns for idum.                 */
            ran_seed = IA * (ran_seed - k * IQ) - IR * k;   /* Compute idum=(IA*idum) % IM without over- */
            if (ran_seed < 0) ran_seed += IM;               /* flows by Schrage’s method.               */
            ans = (double)(AM * ran_seed);                   /* Convert idum to a doubleing result.            */
            /* *idum ^= MASK; */
            /* Unmask before return. */
            return ans;
        }  
    }
}
