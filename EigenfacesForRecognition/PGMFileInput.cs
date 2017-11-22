using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace EigenfacesForRecognition
{
    class PGMFileInput
    {
        private char FirstNonWhitespace(BinaryReader br)
        {
            char c = (char)br.ReadByte();

            while (c == ' ' || c == '\t' || c == '\r' || c == '\n')
                c = (char)br.ReadByte();

            return c;
        }

        private int ReadPositiveNumber(BinaryReader br)
        {
            char c = FirstNonWhitespace(br);
            StringBuilder sb = new StringBuilder();

            if (c >= '0' && c <= '9')
            {
                while (c >= '0' && c <= '9')
                {
                    sb.Append(c);
                    c = (char)br.ReadByte();
                }

                return int.Parse(sb.ToString());
            }

            return -1;
        }

        public bool ReadPGMFile(string filename, out PGMImage image)
        {
            image = null;

            try
            {
                StreamReader sr = new StreamReader(filename);
                BinaryReader br = new BinaryReader(sr.BaseStream);

                // read header beginning with magic 'P', '5'

                byte b = br.ReadByte();

                if ((char)b != 'P')
                    return false;

                b = br.ReadByte();

                if ((char)b != '5')
                    return false;

                int w = ReadPositiveNumber(br);
                int h = ReadPositiveNumber(br);
                int m = ReadPositiveNumber(br);

                if (m <= 0 || m >= 65536)
                    return false;

                byte[,] raster;

                if (m < 256)
                {
                    raster = new byte[h, w];

                    for (int i = 0; i < h; i++)
                        for (int j = 0; j < w; j++)
                            raster[i, j] = br.ReadByte();
                }

                else
                {
                    raster = new byte[h, 2 * w];

                    for (int i = 0; i < h; i++)
                    {
                        for (int j = 0; j < w; j += 2)
                        {
                            raster[i, j + 0] = br.ReadByte();
                            raster[i, j + 1] = br.ReadByte();
                        }
                    }
                }

                image = new PGMImage(raster, h, m, w);
                return true;
            }
            catch (Exception exception)
            {
                MessageBox.Show(exception.ToString(), "Warning",
                    MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return false;
            }
        }
    }
}
