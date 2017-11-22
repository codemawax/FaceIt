using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace EigenfacesForRecognition
{
    public partial class DrawForm : Form
    {
        private int N1, N2;
        private PGMImage avgImage;
        private PGMImage[] image;

        public DrawForm(int N1, int N2, PGMImage avgImage, PGMImage[] image)
        {
            InitializeComponent();
            panel1.Paint += Panel1_Paint;
            panel1.Size = new Size(ClientSize.Width, ClientSize.Height);
            this.N1 = N1;
            this.N2 = N2;
            this.avgImage = avgImage;
            this.image = image;
        }

        private void Panel1_Paint(object sender, PaintEventArgs e)
        {
            int h = image[0].H, M = image.Length, w = image[0].W;
            int count = 0, nrow = N1 / 8;
            Graphics g = e.Graphics;

            for (int i = 0; i <= nrow; i++)
            {
                int left = 8;

                if (i == nrow)
                    left = N1 % 8;

                if (left < 8)
                {
                    for (int j = 0; j < left; j++)
                    {
                        int index = i * 8 * N2 + j * N2;

                        g.DrawImage(image[index].BMap, new Point(j * w, i * h));
                    }
                }

                else
                {
                    for (int j = 0; j < 8; j++)
                    {
                        int index = i * 8 * N2 + j * N2;

                        g.DrawImage(image[index].BMap, new Point(j * w, i * h));
                    }
                }

                count++;
            }

            int x = 0, y = (nrow + 1) * h;

            g.DrawImage(avgImage.BMap, new Point(x, y));
        }
    }
}