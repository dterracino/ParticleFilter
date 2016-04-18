namespace PF
{
    partial class MainForm
    {
        /// <summary>
        /// 必需的设计器变量。
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 清理所有正在使用的资源。
        /// </summary>
        /// <param name="disposing">如果应释放托管资源，为 true；否则为 false。</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows 窗体设计器生成的代码

        /// <summary>
        /// 设计器支持所需的方法 - 不要
        /// 使用代码编辑器修改此方法的内容。
        /// </summary>
        private void InitializeComponent()
        {
            this.pbxVideo = new System.Windows.Forms.PictureBox();
            this.btnSelVideo = new System.Windows.Forms.Button();
            this.btnPause = new System.Windows.Forms.Button();
            this.btnSelRegion = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.pbxVideo)).BeginInit();
            this.SuspendLayout();
            // 
            // pbxVideo
            // 
            this.pbxVideo.Location = new System.Drawing.Point(12, 12);
            this.pbxVideo.Name = "pbxVideo";
            this.pbxVideo.Size = new System.Drawing.Size(678, 522);
            this.pbxVideo.SizeMode = System.Windows.Forms.PictureBoxSizeMode.StretchImage;
            this.pbxVideo.TabIndex = 0;
            this.pbxVideo.TabStop = false;
            this.pbxVideo.Paint += new System.Windows.Forms.PaintEventHandler(this.pbxVideo_Paint);
            this.pbxVideo.MouseDown += new System.Windows.Forms.MouseEventHandler(this.pbxVideo_MouseDown);
            this.pbxVideo.MouseMove += new System.Windows.Forms.MouseEventHandler(this.pbxVideo_MouseMove);
            this.pbxVideo.MouseUp += new System.Windows.Forms.MouseEventHandler(this.pbxVideo_MouseUp);
            // 
            // btnSelVideo
            // 
            this.btnSelVideo.Location = new System.Drawing.Point(696, 21);
            this.btnSelVideo.Name = "btnSelVideo";
            this.btnSelVideo.Size = new System.Drawing.Size(75, 23);
            this.btnSelVideo.TabIndex = 1;
            this.btnSelVideo.Text = "选择视频";
            this.btnSelVideo.UseVisualStyleBackColor = true;
            this.btnSelVideo.Click += new System.EventHandler(this.btnSelVideo_Click);
            // 
            // btnPause
            // 
            this.btnPause.Location = new System.Drawing.Point(696, 63);
            this.btnPause.Name = "btnPause";
            this.btnPause.Size = new System.Drawing.Size(75, 23);
            this.btnPause.TabIndex = 2;
            this.btnPause.Text = "暂停";
            this.btnPause.UseVisualStyleBackColor = true;
            this.btnPause.Click += new System.EventHandler(this.btnPause_Click);
            // 
            // btnSelRegion
            // 
            this.btnSelRegion.Location = new System.Drawing.Point(696, 107);
            this.btnSelRegion.Name = "btnSelRegion";
            this.btnSelRegion.Size = new System.Drawing.Size(75, 23);
            this.btnSelRegion.TabIndex = 3;
            this.btnSelRegion.Text = "选择区域";
            this.btnSelRegion.UseVisualStyleBackColor = true;
            this.btnSelRegion.Click += new System.EventHandler(this.btnSelRegion_Click);
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(783, 546);
            this.Controls.Add(this.btnSelRegion);
            this.Controls.Add(this.btnPause);
            this.Controls.Add(this.btnSelVideo);
            this.Controls.Add(this.pbxVideo);
            this.Name = "MainForm";
            this.Text = "粒子滤波测试";
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.MainForm_FormClosing);
            ((System.ComponentModel.ISupportInitialize)(this.pbxVideo)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.PictureBox pbxVideo;
        private System.Windows.Forms.Button btnSelVideo;
        private System.Windows.Forms.Button btnPause;
        private System.Windows.Forms.Button btnSelRegion;
    }
}

