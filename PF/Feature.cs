using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    abstract class Feature
    {
        const double SIGMA = 0.1;                   /* sigma = 0.1 */
        const double RECIP_SIGMA = 3.9894228040;    /* 1/(sqrt(2*pi)*sigma), 这里sigma = 0.1*/
        const double SIGMA2 = SIGMA * SIGMA * 2;    /* 2*sigma^2 */

        public abstract void Initialize(Image<Bgr, byte> image, int x0, int y0, int wx, int hy);

        public abstract void CalcWeights(Image<Bgr, byte> image, SpaceState[] states, double[] weights);

        public abstract int ModelUpdate(Image<Bgr, byte> image, ref SpaceState estState);

        protected double CalcBhattacharyya(double[] hist1, double[] hist2)
        {
            double rho = 0.0;
            if (hist1.Length != hist2.Length)
                return rho;

            for (int i = 0; i < hist1.Length; i++)
                rho += Math.Sqrt(hist1[i] * hist2[i]);

            return (rho);
        }

        protected double CalcWeightedPi(double rho)
        {
            double pi_n, d2;

            d2 = 1 - rho;
            pi_n = RECIP_SIGMA * Math.Exp(-d2 / SIGMA2);
            //pi_n = Math.Exp(-d2 / SIGMA2);

            return (pi_n);
        }
    }
}
