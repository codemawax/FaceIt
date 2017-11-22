using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EigenfacesForRecognition
{
    public class PGMImage
    {
        private byte[,] raster;
        private int h, m, w;
        private Bitmap bitmap;

        public byte[,] Raster
        {
            get
            {
                return raster;
            }

            set
            {
                raster = value;
            }
        }

        public int H
        {
            get
            {
                return h;
            }
        }

        public int M
        {
            get
            {
                return m;
            }
        }

        public int W
        {
            get
            {
                return w;
            }
        }

        public Bitmap BMap
        {
            get
            {
                return bitmap;
            }
        }

        public PGMImage(byte[,] raster, int h, int m, int w)
        {
            this.raster = raster;
            this.h = h;
            this.m = m;
            this.w = w;
            bitmap = new Bitmap(w, h);

            if (m < 256)
            {
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        int r = raster[i, j];
                        Color c = Color.FromArgb(r, r, r);

                        bitmap.SetPixel(j, i, c);
                    }
                }
            }

            else
            {
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        int r = 256 * raster[i, j] + raster[i, j + 1];
                        Color c = Color.FromArgb(r, r, r);

                        bitmap.SetPixel(j, i, c);
                    }
                }
            }
        }
    }
}
