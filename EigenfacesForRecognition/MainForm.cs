using System;
using System.Windows.Forms;

// petit

namespace EigenfacesForRecognition
{
    public partial class MainForm : Form
    {
        private PGMImage avgImage;
        private PGMImage[] image;
        private PGMImage subjectImage;

        public MainForm()
        {
            InitializeComponent();
        }

        private string Format(double x, string f)
        {
            string s = string.Empty;

            if (x > 0)
                s += "+";

            else if (x == 0)
                s += " ";

            return s + x.ToString(f) + " ";
        }

        private void button1_Click(object sender, EventArgs e)
        {
            int N1 = (int)numericUpDown1.Value;
            int N2 = (int)numericUpDown2.Value;
            int M = N1 * N2, MP = N1, count = 0;
            string baseName = "Database";

            image = new PGMImage[M];

            for (int i = 0; i < N1; i++)
            {
                string dirName = baseName + "\\" + "s" + (i + 1).ToString() + "\\";

                for (int j = 0; j < N2; j++)
                {
                    string fileName = dirName + (j + 1).ToString() + ".pgm";
                    PGMFileInput input = new PGMFileInput();

                    if (input.ReadPGMFile(fileName, out image[count++]))
                    {

                    }
                }
            }

            int h = image[0].H, m = image[0].M, w = image[0].W, m1 = m + 1;
            byte[,] avgRaster;
            long[,] lar;

            if (m < 256)
            {
                avgRaster = new byte[h, w];
                lar = new long[h, w];

                for (int i = 0; i < M; i++)
                    for (int j = 0; j < h; j++)
                        for (int k = 0; k < w; k++)
                            lar[j, k] += image[i].Raster[j, k];

                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < h; j++)
                    {
                        for (int k = 0; k < w; k++)
                        {
                            long r = lar[j, k] / M;

                            avgRaster[j, k] = (byte)r;
                        }
                    }
                }
            }

            else
            {
                avgRaster = new byte[h, 2 * w];
                lar = new long[h, 2 * w];

                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < h; j++)
                    {
                        for (int k = 0; k < w; k++)
                        {
                            long rh = image[i].Raster[j, k + 0];
                            long rl = image[i].Raster[j, k];
                            long r = 256 * rh + rl;

                            lar[j, k] += r;
                        }
                    }
                }

                for (int i = 0; i < M; i++)
                    for (int j = 0; j < h; j++)
                        for (int k = 0; k < w; k++)
                            lar[j, k] /= M;

                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < h; j++)
                    {
                        for (int k = 0; k < w; k += 2)
                        {
                            long srh = image[i].Raster[j, k + 0];
                            long srl = image[i].Raster[j, k + 1];
                            long sr = 256 * srh + srl;
                            long r = (lar[j, k] + sr) % m1;

                            srh = r / 256;
                            srl = r % 256;
                            avgRaster[j, k + 0] = (byte)srh;
                            avgRaster[j, k + 1] = (byte)srl;
                        }
                    }
                }
            }

            avgImage = new PGMImage(avgRaster, h, m, w);

            DrawForm df = new DrawForm(N1, N2, avgImage, image);

            df.Show();

            byte[,,] phiRaster;
            PGMImage[] phi = new PGMImage[M];

            if (m < 256)
            {
                phiRaster = new byte[M, h, w];

                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < h; j++)
                    {
                        for (int k = 0; k < w; k++)
                        {
                            int s = image[i].Raster[j, k];
                            int r = avgRaster[j, k];
                            int sr = (s - r) % m1;

                            if (sr < 0)
                                sr += m1;

                            phiRaster[i, j, k] = (byte)sr;
                        }
                    }
                }
            }

            else
            {
                phiRaster = new byte[M, h, 2 * w];

                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < h; j++)
                    {
                        for (int k = 0; k < w; k += 2)
                        {
                            int arh = avgRaster[j, k + 0], arl = avgRaster[j, k + 1];
                            int srh = image[i].Raster[j, k + 0], srl = image[i].Raster[j, k + 1];
                            int ar = 256 * arh + arl, sr = 256 * srh + srl;
                            int r = (sr - ar) % m1;

                            if (r < 0)
                                r += m1;

                            arh = r / 256;
                            arl = r % 256;
                            avgRaster[j, k + 0] = (byte)arh;
                            avgRaster[j, k + 1] = (byte)arl;
                        }
                    }
                }
            }

            byte[,] raster = new byte[h, w];
            
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < h; j++)
                    for (int k = 0; k < w; k++)
                        raster[j, k] = phiRaster[i, j, k];

                phi[i] = new PGMImage(raster, h, m, w);
            }

            df = new DrawForm(N1, N2, avgImage, phi);
            df.Show();

            EigenVV evv = new EigenVV(phiRaster, M, h, m, w);
            double[,] vec;
            double[] val = evv.FindEigenValuesVectors(out vec, out count);
            string f = "E6";

            if (count != 0)
            {
                textBox1.Text += "Error in calculating eigenvalues\r\n";
                return;
            }

            // use selection sort to sort the eigenvalues and eigen vectors
            // into descending order

            for (int i = 1; i <= M; i++)
            {
                for (int j = 1; j <= M; j++)
                {
                    if (val[i] > val[j])
                    {
                        double t = val[i];

                        val[i] = val[j];
                        val[j] = t;

                        for (int k = 1; k <= M - count; k++)
                        {
                            t = vec[i, k];
                            vec[i, k] = vec[j, k];
                            vec[j, k] = t;
                        }
                    }
                }
            }

            for (int i = 1; i <= MP; i++)
                textBox1.Text += Format(val[i], f) + "\r\n";

<<<<<<< HEAD
            //Stream myStream = null;
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = "c:\\";
            openFileDialog1.Filter = "PGM files (*.pgm)|*.pgm";
            openFileDialog1.FilterIndex = 2;
            openFileDialog1.RestoreDirectory = true;

            if (openFileDialog1.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                System.IO.StreamReader sr = new
                   System.IO.StreamReader(openFileDialog1.FileName);
                MessageBox.Show(sr.ReadToEnd());
                sr.Close();

                //PictureBox PictureBox1 = new PictureBox();
                //PictureBox1.Image = new PGMImage(openFileDialog1.FileName);
                //this.Controls.Add(PictureBox1);
            }

            //df = new DrawForm(subjectImage);
            //df.Show();
=======
            DistanceComputation distanceComputation = new DistanceComputation(10, 10);
            double[] distances;
            distances = new double[N1];

            double[] eigenVectors;

            eigenVectors = new double[];

            int k;
            for (k = 0; k < N1; k++)
            {
                distances = distanceComputation.ComputeDistance(avgImage, subjectImage, eigenVectors)
            }

>>>>>>> 003449ab2e675bfe4aa72c489ab6e16929a14432
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void MainForm_Load(object sender, EventArgs e)
        {

        }
    }
}