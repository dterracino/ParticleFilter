using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace PF
{
    public partial class MainForm : Form
    {
        const int MAX_COUNT = 30;
        Capture cap;
        ParticleFilter pf;
        bool track = false;

        // 做图相关
        bool selectRegion = false;
        bool isLeftDown = false;
        Point leftTop = new Point();
        Point rightBottom = new Point();
        int width, height;
        Image<Bgr, byte> currentFrame;
        Thread thread;
        bool pause = false;
        int lostCount = 0;

        public MainForm()
        {
            InitializeComponent();
        }

        private void btnSelVideo_Click(object sender, EventArgs e)
        {
            OpenFileDialog fDialog = new OpenFileDialog();
            fDialog.CheckFileExists = true;
            fDialog.InitialDirectory = @"F:\车辆";
            fDialog.Filter = "MP4文件|*.mp4|avi文件|*.avi";
            fDialog.Multiselect = false;
            if (fDialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                // 关闭上一个视频
                if (cap != null)
                {
                    thread.Abort();
                    thread = null;
                    track = false;
                }
                
                cap = new Capture(fDialog.FileName);
                width = (int)CvInvoke.cvGetCaptureProperty(cap, Emgu.CV.CvEnum.CAP_PROP.CV_CAP_PROP_FRAME_WIDTH);
                height = (int)CvInvoke.cvGetCaptureProperty(cap, Emgu.CV.CvEnum.CAP_PROP.CV_CAP_PROP_FRAME_HEIGHT);
                var fps = (int)CvInvoke.cvGetCaptureProperty(cap, Emgu.CV.CvEnum.CAP_PROP.CV_CAP_PROP_FPS);
                thread = new Thread(() =>
                {
                    while (true)
                    {
                        Thread.Sleep(1);
                        if (pause)
                            continue;

                        // 循环读入视频并使用粒子滤波跟踪
                        if ((currentFrame = cap.QueryFrame()) == null) break;


                        int rho_v;
                        int xc, yc, wh_x, hy_h;
                        double max_weight;
                        if (track)
                        {
                            // do track
                            rho_v = pf.ColorParticleTracking(currentFrame, out xc, out yc, out wh_x, out hy_h, out max_weight);
                            if (rho_v > 0 && max_weight > 0.0001)  /* 判断是否目标丢失 */
                            {
                                // draw blue
                                currentFrame.Draw(new Rectangle(xc - wh_x, yc - hy_h, wh_x * 2, hy_h * 2), new Bgr(Color.Blue), 1);
                                lostCount = 0;
                            }
                            else
                            {
                                // draw blue
                                currentFrame.Draw(new Rectangle(xc - wh_x, yc - hy_h, wh_x * 2, hy_h * 2), new Bgr(Color.Red), 1);
                                lostCount++;
                                if (lostCount >= MAX_COUNT)
                                {
                                    track = false;
                                    pf = null;
                                }
                            }
                        }
                        this.Invoke(new Action(() =>
                        {
                            pbxVideo.Image = currentFrame.ToBitmap();
                        }));
                    }
                });
                thread.Start();
            }
        }

        private void StartParticleFilter()
        {
            pf = new ParticleFilter(60);
            // 首先计算框在图中的坐标   
            int xc = (leftTop.X + rightBottom.X) / 2;
            int yc = (leftTop.Y + rightBottom.Y) / 2;
            int width_h = (rightBottom.X - leftTop.X)/ 2;
            int height_h = (rightBottom.Y - leftTop.Y) / 2;

            int x0 = xc * width / pbxVideo.Width;
            int y0 = yc * height / pbxVideo.Height;
            int wx = width_h * width / pbxVideo.Width;
            int hy = height_h * height / pbxVideo.Height;
            pf.Initialize(x0, y0, wx, hy, currentFrame);
            track = true;
        }

        void cap_ImageGrabbed(object sender, EventArgs e)
        {
            // 循环读入视频并使用粒子滤波跟踪
            currentFrame = (sender as Capture).RetrieveBgrFrame();

            int rho_v;
            int xc, yc, wh_x, hy_h;
            double max_weight;
            if (track)
            {
                // do track
                rho_v = pf.ColorParticleTracking(currentFrame, out xc, out yc, out wh_x, out hy_h, out max_weight);
                // draw
                if (rho_v > 0 && max_weight > 0.0001)  /* 判断是否目标丢失 */
                {
                    currentFrame.Draw(new Rectangle(xc - wh_x, yc - hy_h, wh_x * 2, hy_h * 2), new Bgr(Color.Red), 1);
                }

            }

            pbxVideo.Image = currentFrame.ToBitmap();
            Thread.Sleep(100);
        }

        private void btnSelRegion_Click(object sender, EventArgs e)
        {
            selectRegion = true;
        }

        private void pbxVideo_MouseDown(object sender, MouseEventArgs e)
        {
            if (selectRegion)
            {
                leftTop.X = e.X;
                leftTop.Y = e.Y;
                isLeftDown = true;
            }
        }

        private void pbxVideo_MouseMove(object sender, MouseEventArgs e)
        {
            if (selectRegion && isLeftDown)
            {
                rightBottom.X = e.X;
                rightBottom.Y = e.Y;
                pbxVideo.Refresh();
            }
        }

        private void pbxVideo_Paint(object sender, PaintEventArgs e)
        {
            if (selectRegion && isLeftDown)
            {
                var g = e.Graphics;
                g.DrawRectangle(Pens.Blue,
                    new Rectangle(leftTop.X, leftTop.Y, rightBottom.X - leftTop.X, rightBottom.Y - leftTop.Y));
            }
        }

        private void pbxVideo_MouseUp(object sender, MouseEventArgs e)
        {
            if (selectRegion && isLeftDown)
            {
                selectRegion = isLeftDown = false;
                // 选择了一个区域，开始粒子滤波
                StartParticleFilter();
            }
        }

        private void btnPause_Click(object sender, EventArgs e)
        {
            if ((sender as Button).Text == "暂停")
            {
                (sender as Button).Text = "开始";
            }
            else 
            {
                (sender as Button).Text = "暂停";
            }
            pause = !pause;
        }

        private void MainForm_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (thread != null && thread.IsAlive)
            {
                thread.Abort();
            }
        }
    }
}
