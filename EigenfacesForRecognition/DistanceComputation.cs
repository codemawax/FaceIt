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
        private PGMImage m_diffImage, m_reconstructedEigenImage;

        public DistanceComputation(int height, int width)
        {
            this.m_height = height;
            this.m_width = width;
            this.m_pixelCount = height * width;
        }

        /**
       * Computes the EigenFaces matrix for a dataset of Eigen vectors and a diff matrix.
       * @param eigenVectors
       * @param diffMatrix
       */
        private void ComputeEigenFaces(double[] ai_eigenVectors)
        {
            m_imageCount = ai_eigenVectors.Length;
            m_rank = ai_eigenVectors.Length; // Get height

            int n, k, i, j;

            for (n = 0; n < m_rank; n++)
            {
                for (k = 0; k < m_eigenFaceCount; k++)
                {
                    double sumSquare = 0.0F;
                    for (i = 0; i < m_height; i++)
                    {
                        for (j = 0; j < m_width; j++)
                        {
                            sumSquare += (byte)(m_diffImage.Raster[i, j] * ai_eigenVectors[k]);
                        }
                    }
                }
            }


//        (0 to(rank - 1)).foreach {
//                i =>
//var sumSquare = 0.0
//(0 to(pixelCount - 1)).foreach {
//                    j =>
//(0 to(imageCount - 1)).foreach {
//                        k =>
//eigenFaces(j)(i) += diffMatrix(j)(k) * eigenVectors.get(i, k)
//        }
//                    sumSquare += eigenFaces(j)(i) * eigenFaces(j)(i)
//      }
//                var norm = Math.sqrt(sumSquare)
//                (0 to(pixelCount - 1)).foreach {
//                    j =>
//eigenFaces(j)(i) /= norm
//                }
//            }
//            val eigenFacesMatrix = new DenseDoubleMatrix2D(pixelCount, rank)
//        eigenFacesMatrix.assign(eigenFaces)
      }

        /**
         * Calculates a distance score between a mean Pixels/EigenFaces model in comparison to an image subject.
         * @param meanPixels
         * @param eigenFaces
         * @param subjectPixels
         */
        public double ComputeDistance(PGMImage ai_avgImage, PGMImage ai_subjectImage, double[] ai_eigenVectors)
        {
            ComputeDiffImage(ai_subjectImage, ai_avgImage);
            ComputeEigenFaces(ai_eigenVectors);
            ComputeWeights();
            ReconstructImageWithEigenFaces(ai_avgImage);
            return ComputeImageDistance(ai_subjectImage, m_reconstructedEigenImage);
        }

        /**
         * Computes the distance between two images.
         * @param pixels1
         * @param pixels2
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
         * @param diffImagePixels
         * @param eigenFaces
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
         * @param subjectPixels
         * @param meanPixels
         */
        private void ComputeDiffImage(PGMImage ai_subjectImage, PGMImage ai_avgImage)
        {
            int i, j;
            for (i = 0; i < m_height; i++)
            {
                for (j = 0; j < m_width; j++)
                {
                    m_diffImage.Raster[i, j] = (byte)(ai_subjectImage.Raster[i, j] - ai_avgImage.Raster[i, j]);
                }
            }
        }

        /**
         * Reconstructs an image using Eigen Faces and weights.
         * @param weights
         * @param eigenFaces
         * @param meanPixels
         */
        private void ReconstructImageWithEigenFaces(PGMImage ai_avgImage)
        {
            int k, i, j;
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

            // TODO normalize


//            var min = Double.MaxValue
//    var max = -Double.MaxValue
//    (0 to(reconstructedPixels.length - 1)).foreach {
//                i =>
//min = Math.min(min, reconstructedPixels(i))
//      max = Math.max(max, reconstructedPixels(i))
//    }

//            val normalizedReconstructedPixels = new Array[Double](pixelCount)
//            (0 to(reconstructedPixels.length - 1)).foreach {
//                i =>
//normalizedReconstructedPixels(i) = (255.0 * (reconstructedPixels(i) - min)) / (max - min)
//            }
//            normalizedReconstructedPixels


  }



}

}

