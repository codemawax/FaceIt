using System;
using EigenfacesForRecognition;

namespace EigenfacesForRecognition
{
    public class DistanceComputation
    {
        private int m_height, m_width, m_pixelCount, m_imageCount, m_rank, m_eigenFaceCount;
        private double[] m_weights;

        // Allocate ?
        private PGMImage[] m_eigenFaces;
        private PGMImage m_avgImage, m_diffImage, m_reconstructedEigenImage, m_normalizedReconstructedEigenImage;

        public DistanceComputation(int height, int width, int m)
        {
            this.m_height = height;
            this.m_width = width;
            this.m_pixelCount = height * width;

            byte[,] raster = new byte[height, width];
            m_avgImage = new PGMImage(raster, height, m, width);
            m_diffImage = new PGMImage(raster, height, m, width);
        }

        /**
 * Computes the EigenFaces matrix using a pixel matrix of multiple images.
 * @param pixelMatrix
 * @param meanColumn
 */
        //      def computeEigenFaces(pixelMatrix: Array[Array[Double]], meanColumn: Array[Double]): DoubleMatrix2D = {
        //  val diffMatrix = MatrixHelpers.computeDifferenceMatrixPixels(pixelMatrix, meanColumn)
        //  val eigenVectors = MatrixHelpers.computeEigenVectors(diffMatrix)
        //  computeEigenFaces(eigenVectors, diffMatrix)
        //}
     
        /**
       * Computes the EigenFaces matrix for a dataset of Eigen vectors and a diff matrix.
       */
        private void ComputeEigenFaces(PGMImage[] ai_diffImages, double[,] ai_eigenVectors)
        {
            m_imageCount = ai_eigenVectors.Length;
            m_rank = ai_eigenVectors.Length; // get height
            m_eigenFaceCount = ai_eigenVectors.Length; 

            int n, k, i, j;

            for (n = 0; n < m_rank; n++)
            {
                double sumSquare = 0.0F;

                for (i = 0; i < m_height; i++)
                {
                    for (j = 0; j < m_width; j++)
                    {
                        double eigenFacePixel = 0.0F;
                        for (k = 0; k < m_eigenFaceCount; k++)
                        {
                            Console.WriteLine("face " + k);
                            eigenFacePixel += (byte)(ai_diffImages[k].Raster[i, j] * ai_eigenVectors[n, k]);
                        }
                        m_eigenFaces[n].Raster[i, j] = (byte)eigenFacePixel;
                        sumSquare += (byte)(Math.Pow(eigenFacePixel, 2));
                    }
                }
                double norm = Math.Sqrt(sumSquare);
                for (i = 0; i < m_height; i++)
                {
                    for (j = 0; j < m_width; j++)
                    {
                        m_eigenFaces[n].Raster[i, j] = (byte)(m_eigenFaces[n].Raster[i, j] / norm);
                    }
                }
            }
        }

        /**
         * Calculates a distance score between a mean Pixels/EigenFaces model in comparison to an image subject.
         */
        public double ComputeDistance(PGMImage[] ai_diffImages, PGMImage ai_avgImage, PGMImage ai_subjectImage, double[,] ai_eigenVectors)
        {
            ComputeDiffImage(ref ai_subjectImage, ref ai_avgImage, ref m_diffImage);
            ComputeEigenFaces(ai_diffImages, ai_eigenVectors);
            ComputeWeights();
            ReconstructImageWithEigenFaces(ai_avgImage);
            return ComputeImageDistance(ai_subjectImage, m_reconstructedEigenImage);
        }

        /**
         * Computes the distance between two images.
         */
        private double ComputeImageDistance(PGMImage ai_image1, PGMImage ai_image2)
        {
            double distance = 0.0F;
            int i, j;
            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    distance += Math.Pow(ai_image1.Raster[i, j] - ai_image2.Raster[i, j], 2);
                }
            }
            return Math.Sqrt(distance / m_pixelCount);
        }

        /**
         * Computes the weights of faces vs. EigenFaces.
         */
        private void ComputeWeights()
        {
            int k, i, j;
            for (k = 0; k < m_eigenFaceCount; k++)
            {
                m_weights[k] = 0.0F;
                for (i = 0; i < m_height; i++)
                {
                    for (j = 0; j < m_width; j++)
                    {
                        m_weights[i] += m_diffImage.Raster[i, j] * m_eigenFaces[k].Raster[i, j];
                    }
                }
            }
        }

        /**
         * Computes the difference pixels between a subject image and a mean image.
         */
        private void ComputeDiffImage(ref PGMImage ai_subjectImage, ref PGMImage ai_avgImage, ref PGMImage ao_diffImage)
        {
            int i, j;
            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    ao_diffImage.Raster[i, j] = (byte)(ai_subjectImage.Raster[i, j] - ai_avgImage.Raster[i, j]);
                }
            }
        }

        /**
         * Reconstructs an image using Eigen Faces and weights.
         */
        private void ReconstructImageWithEigenFaces(PGMImage ai_avgImage)
        {
            int k, i, j;
            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    m_reconstructedEigenImage.Raster[i, j] = 0;
                }
            }

            for (k = 0; k < m_eigenFaceCount; k++)
            {
                m_weights[k] = 0.0F;
                for (i = 0; i < m_height; i++)
                {
                    for (j = 0; j < m_width; j++)
                    {
                        m_reconstructedEigenImage.Raster[i, j] += (byte)(m_diffImage.Raster[i, j] * m_eigenFaces[k].Raster[i, j]);
                    }
                }
            }

            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    m_reconstructedEigenImage.Raster[i, j] += ai_avgImage.Raster[i, j];
                }
            }

            double min = double.MaxValue;
            double max = double.MinValue;

            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    min = Math.Min(min, m_reconstructedEigenImage.Raster[i, j]);
                    max = Math.Max(max, m_reconstructedEigenImage.Raster[i, j]);
                }
            }

            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    m_normalizedReconstructedEigenImage.Raster[i, j] = (byte)(255.0 * (m_reconstructedEigenImage.Raster[i, j] - min) / (max - min));
                }
            }
        }
    }
}

