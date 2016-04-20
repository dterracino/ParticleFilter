using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PF
{
    class Haar:Feature
    {
        /* 目标模型 */
        double[] modelHist;

        public override void Initialize(Image<Bgr, byte> image, int x0, int y0, int wx, int hy)
        {
            throw new NotImplementedException();
        }

        public override void CalcWeights(Image<Bgr, byte> image, SpaceState[] particles, double[] weights)
        {
            throw new NotImplementedException();
        }

        public override int ModelUpdate(Image<Bgr, byte> image, SpaceState estState)
        {
            throw new NotImplementedException();
        }
    }
}
