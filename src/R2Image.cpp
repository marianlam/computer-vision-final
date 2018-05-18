// Source file for image class

// Include files
#include <algorithm>
#include <cmath>
#include <math.h>
#include <iostream>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include <stdlib.h>
#include "svd.h"
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
: pixels(NULL),
npixels(0),
width(0),
height(0)
{
}



R2Image::
R2Image(const char *filename)
: pixels(NULL),
npixels(0),
width(0),
height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
: pixels(NULL),
npixels(width * height),
width(width),
height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
: pixels(NULL),
npixels(width * height),
width(width),
height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
  pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
: pixels(NULL),
npixels(image.npixels),
width(image.width),
height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
  pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
  pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


double** R2Image::
svdTest(vector<int> x, vector<int> y, vector<int> xd, vector<int> yd, int size)
{
  // returns the H homography matrix computed using the DLT algorithm
  int n = size;

  double** matrixA = dmatrix(1,2*n,1,9);

  for (int i = 1; i <= 2*n; i++) {
    if (i % 2 != 0) {
      int index = (i-1)/2;
      matrixA[i][1] = 0;
      matrixA[i][2] = 0;
      matrixA[i][3] = 0;
      matrixA[i][4] = -x[index];
      matrixA[i][5] = -y[index];
      matrixA[i][6] = -1;
      matrixA[i][7] = yd[index]*x[index];
      matrixA[i][8] = yd[index]*y[index];
      matrixA[i][9] = yd[index];
    } else {
      int index = (i-2)/2;
      matrixA[i][1] = x[index];
      matrixA[i][2] = y[index];
      matrixA[i][3] = 1;
      matrixA[i][4] = 0;
      matrixA[i][5] = 0;
      matrixA[i][6] = 0;
      matrixA[i][7] = -xd[index]*x[index];
      matrixA[i][8] = -xd[index]*y[index];
      matrixA[i][9] = -xd[index];
    }
  }

  // // fit a 2D conic to five points
  // R2Point p1(1.2,3.5);
  // R2Point p2(2.1,2.2);
  // R2Point p3(0.2,1.6);
  // R2Point p4(0.0,0.5);
  // R2Point p5(-0.2,4.2);

  // // build the 5x6 matrix of equations
  // double** linEquations = dmatrix(1,5,1,6);
  //
  // linEquations[1][1] = p1[0]*p1[0];
  // linEquations[1][2] = p1[0]*p1[1];
  // linEquations[1][3] = p1[1]*p1[1];
  // linEquations[1][4] = p1[0];
  // linEquations[1][5] = p1[1];
  // linEquations[1][6] = 1.0;
  //
  // linEquations[2][1] = p2[0]*p2[0];
  // linEquations[2][2] = p2[0]*p2[1];
  // linEquations[2][3] = p2[1]*p2[1];
  // linEquations[2][4] = p2[0];
  // linEquations[2][5] = p2[1];
  // linEquations[2][6] = 1.0;
  //
  // linEquations[3][1] = p3[0]*p3[0];
  // linEquations[3][2] = p3[0]*p3[1];
  // linEquations[3][3] = p3[1]*p3[1];
  // linEquations[3][4] = p3[0];
  // linEquations[3][5] = p3[1];
  // linEquations[3][6] = 1.0;
  //
  // linEquations[4][1] = p4[0]*p4[0];
  // linEquations[4][2] = p4[0]*p4[1];
  // linEquations[4][3] = p4[1]*p4[1];
  // linEquations[4][4] = p4[0];
  // linEquations[4][5] = p4[1];
  // linEquations[4][6] = 1.0;
  //
  // linEquations[5][1] = p5[0]*p5[0];
  // linEquations[5][2] = p5[0]*p5[1];
  // linEquations[5][3] = p5[1]*p5[1];
  // linEquations[5][4] = p5[0];
  // linEquations[5][5] = p5[1];
  // linEquations[5][6] = 1.0;
  //
  // printf("\n Fitting a conic to five points:\n");
  // printf("Point #1: %f,%f\n",p1[0],p1[1]);
  // printf("Point #2: %f,%f\n",p2[0],p2[1]);
  // printf("Point #3: %f,%f\n",p3[0],p3[1]);
  // printf("Point #4: %f,%f\n",p4[0],p4[1]);
  // printf("Point #5: %f,%f\n",p5[0],p5[1]);
  //

  // compute the SVD
  double singularValues[10]; // 1..9
  double** nullspaceMatrix = dmatrix(1,9,1,9);

  svdcmp(matrixA, 2*n, 9, singularValues, nullspaceMatrix);

  // get the result
  // printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for (int i = 2; i < 10; i++) if (singularValues[i] < singularValues[smallestIndex]) smallestIndex = i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  // printf("\n Transformation Matrix H: \n %f, %f, %f,\n %f, %f, %f,\n %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);

  double** matrixH = dmatrix(1,3,1,3);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      switch(i) {
        case 1:
        matrixH[i][j] = nullspaceMatrix[j][smallestIndex];
        break;
        case 2:
        matrixH[i][j] = nullspaceMatrix[j+3][smallestIndex];
        break;
        case 3:
        matrixH[i][j] = nullspaceMatrix[j+6][smallestIndex];
        break;
      }
    }
  }

  // // make sure the solution is correct:
  // printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] +
  // p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] +
  // p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] +
  // p1[0]*nullspaceMatrix[4][smallestIndex] +
  // p1[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);
  //
  // printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] +
  // p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] +
  // p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] +
  // p2[0]*nullspaceMatrix[4][smallestIndex] +
  // p2[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);
  //
  // printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] +
  // p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] +
  // p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] +
  // p3[0]*nullspaceMatrix[4][smallestIndex] +
  // p3[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);
  //
  // printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] +
  // p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] +
  // p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] +
  // p4[0]*nullspaceMatrix[4][smallestIndex] +
  // p4[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);
  //
  // printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] +
  // p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] +
  // p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] +
  // p5[0]*nullspaceMatrix[4][smallestIndex] +
  // p5[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);
  //
  // R2Point test_point(0.34,-2.8);
  //
  // printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
  // test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
  // test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
  // test_point[0]*nullspaceMatrix[4][smallestIndex] +
  // test_point[1]*nullspaceMatrix[5][smallestIndex] +
  // nullspaceMatrix[6][smallestIndex]);

  return matrixH;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}


void R2Image::
SobelX(void)
{
  // convert image to grayscale
  (*this).ChangeSaturation(0);
  // make copy of image
  R2Image temp(*this);
  // initialize kernel
  const int SOBEL_X_KERNEL[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
  // Apply the Sobel oprator to the image in X direction
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j) = SOBEL_X_KERNEL[0][0] * temp.Pixel(i-1,j-1)
      + SOBEL_X_KERNEL[0][2] * temp.Pixel(i-1,j+1)
      + SOBEL_X_KERNEL[1][0] * temp.Pixel(i,j-1)
      + SOBEL_X_KERNEL[1][2] * temp.Pixel(i,j+1)
      + SOBEL_X_KERNEL[2][0] * temp.Pixel(i+1,j-1)
      + SOBEL_X_KERNEL[2][2] * temp.Pixel(i+1,j+1);
      // Pixel(i,j).Clamp();
    }
  }
}


void R2Image::
SobelY(void)
{
  // convert image to grayscale
  (*this).ChangeSaturation(0);
  // make copy of image
  R2Image temp(*this);
  // initialize kernel
  const int SOBEL_Y_KERNEL[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
  // Apply the Sobel oprator to the image in Y direction
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j) = SOBEL_Y_KERNEL[0][0] * temp.Pixel(i-1,j-1)
      + SOBEL_Y_KERNEL[0][1] * temp.Pixel(i-1,j)
      + SOBEL_Y_KERNEL[0][2] * temp.Pixel(i-1,j+1)
      + SOBEL_Y_KERNEL[2][0] * temp.Pixel(i+1,j-1)
      + SOBEL_Y_KERNEL[2][1] * temp.Pixel(i+1,j)
      + SOBEL_Y_KERNEL[2][2] * temp.Pixel(i+1,j+1);
      // Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}

void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {

      double red = Pixel(i,j).Red();
      double green = Pixel(i,j).Green();
      double blue = Pixel(i,j).Blue();
      double gray = (red+green+blue)/3;

      // extrapolate grayness
      double newRed = -gray*(factor-1)+red*factor;
      double newGreen = -gray*(factor-1)+green*factor;
      double newBlue = -gray*(factor-1)+blue*factor;

      // normalize
      if (newRed > 1) newRed = 1;
      if (newGreen > 1) newGreen = 1;
      if (newBlue > 1) newBlue = 1;

      if (newRed < 0) newRed = 0;
      if (newGreen < 0) newGreen = 0;
      if (newBlue < 0) newBlue = 0;

      // set new rgb
      Pixel(i,j).SetRed(newRed);
      Pixel(i,j).SetGreen(newGreen);
      Pixel(i,j).SetBlue(newBlue);
    }
  }
}

// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // make copy of image
  R2Image temp(*this);
  // calculate kernel width
  int kernel_width = 6 * sigma + 1;
  // initialize vector
  vector<float> kernel;
  // compute kernel (set of weights)
  int midpoint = kernel_width / 2;
  float weight = 0;
  for (int i = 0; i < kernel_width; i++) {
    weight = exp(-(pow(i - midpoint, 2) / (2 * sigma * sigma)))
    / sqrt(2 * M_PI * sigma * sigma);
    kernel.push_back(weight);
  }
  // sum weights
  float sum = 0;
  for (int i = 0; i < kernel_width; i++) {
    // print weights
    // cout << kernel[i] << endl;
    sum += kernel[i];
  }
  // print sum
  // cout << "sum = " << sum << endl;

  float normalized_sum = 0;
  // normalize kernel
  for (int i = 0; i < kernel_width; i++) {
    kernel[i] = kernel[i] / sum;
    normalized_sum += kernel[i];
  }
  // print normalized sum
  // cout << "normalized sum = " << normalized_sum << endl; // should equal 1

  // Gaussian blur of the image. Separable solution is preferred
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      // horizontal pass
      R2Pixel *pixel_copy = new R2Pixel();
      for (int k = 0; k < kernel_width; k++) {
        int index = j-midpoint+k;
        // prevent looping out of bounds
        if (index < 0) {
          index = 0;
        }
        if (index >= height) {
          index = height - 1;
        }
        *pixel_copy += kernel[k] * temp.Pixel(i,index);
      }
      Pixel(i,j) = *pixel_copy;
      // Pixel(i,j).Clamp();
    }
  }

  // make another copy of image
  R2Image temp2(*this);

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      // vertical pass
      R2Pixel *pixel_copy = new R2Pixel();
      for (int k = 0; k < kernel_width; k++) {
        int index = i-midpoint+k;
        // prevent looping out of bounds
        if (index < 0) {
          index = 0;
        }
        if (index >= width) {
          index = width - 1;
        }
        *pixel_copy += kernel[k] * temp2.Pixel(index,j);
      }
      Pixel(i,j) = *pixel_copy;
      // Pixel(i,j).Clamp();
    }
  }
}

// EXTRA CREDIT #1
void R2Image::
MedianFilter(void)
{
  R2Pixel medianPixel;
  vector<Tuple> sampling_window;
  for (int i = 0; i < 9; i++) {
    sampling_window.push_back(Tuple(0, R2Pixel(0.0,0.0,0.0,0.0)));
  }
  // Implementation of a median filter
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      sampling_window[0] = Tuple(Pixel(i-1,j-1).Red(), Pixel(i-1,j-1));
      sampling_window[1] = Tuple(Pixel(i-1,j).Red(), Pixel(i-1,j));
      sampling_window[2] = Tuple(Pixel(i-1,j+1).Red(), Pixel(i-1,j+1));

      sampling_window[3] = Tuple(Pixel(i,j-1).Red(), Pixel(i,j-1));
      sampling_window[4] = Tuple(Pixel(i,j).Red(), Pixel(i,j));
      sampling_window[5] = Tuple(Pixel(i,j+1).Red(), Pixel(i,j+1));

      sampling_window[6] = Tuple(Pixel(i+1,j-1).Red(), Pixel(i+1,j-1));
      sampling_window[7] = Tuple(Pixel(i+1,j).Red(), Pixel(i+1,j));
      sampling_window[8] = Tuple(Pixel(i+1,j+1).Red(), Pixel(i+1,j+1));

      // sort vector in descending order
      sort(sampling_window.rbegin(), sampling_window.rend());
      medianPixel = sampling_window[4].pixel;
      SetPixel(i,j,medianPixel);
    }
  }
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      if (i == 0 || (i == width - 1) || j == 0 || (j == height - 1)) {
        SetPixel(i,j,R2Pixel(1.0,1.0,1.0,1.0));
      }
    }
  }
}

void R2Image::
BilateralFilter(void) {
  // Implementation of a bilateral filter
}

void R2Image::
Harris(double sigma)
{
  // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
  // Output should be 50% grey at flat regions, white at corners and black/dark near edges
  recentLocations.clear();

  R2Image img1(*this);
  R2Image img2(*this);
  R2Image img3(*this);
  R2Image temp(*this);

  img1.SobelX();
  img2.SobelY();
  img3.SobelX();

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      img3.Pixel(i,j) = img3.Pixel(i,j) * img2.Pixel(i,j); // Ixy
      img1.Pixel(i,j) = img1.Pixel(i,j) * img1.Pixel(i,j); // Ixx
      img2.Pixel(i,j) = img2.Pixel(i,j) * img2.Pixel(i,j); // Iyy
    }
  }

  img1.Blur(sigma);
  img2.Blur(sigma);
  img3.Blur(sigma);

  vector<Feature> featureVec;

  for (int i = 0.1 * width; i < (width - (0.1 * width)); i++) {
    for (int j = 0.1 * height; j < (height - (0.1 * height)); j++) {
      temp.Pixel(i,j) = img1.Pixel(i,j)*img2.Pixel(i,j)-img3.Pixel(i,j)*img3.Pixel(i,j)
      -0.04*((img1.Pixel(i,j)+img2.Pixel(i,j))*(img1.Pixel(i,j)+img2.Pixel(i,j)));

      featureVec.push_back(Feature(i, j, temp.Pixel(i,j)));

      temp.Pixel(i,j).SetRed(temp.Pixel(i,j).Red()+0.5);
      temp.Pixel(i,j).SetGreen(temp.Pixel(i,j).Green()+0.5);
      temp.Pixel(i,j).SetBlue(temp.Pixel(i,j).Blue()+0.5);
      temp.Pixel(i,j).SetAlpha(temp.Pixel(i,j).Alpha()+0.5);
      temp.Pixel(i,j).Clamp();
    }
  }

  // sort vector in descending order
  sort(featureVec.rbegin(), featureVec.rend());

  // add 150 elements with highest harris scores to recentLocations vector
  R2Pixel *redPixel = new R2Pixel(1,0,0,0);
  Feature currentFeature;
  int counter = 0;
  int index = 0;
  bool window_isEmpty;
  
    while (counter < 150) {
    currentFeature = featureVec[index];
    window_isEmpty = true;
    for (int i = -10; i < 10; i++) {
      for (int j = -10; j <= 10; j++) {
        if ((currentFeature.centerX + i > 0)
        && (currentFeature.centerX + i < width)
        && (currentFeature.centerY + j > 0)
        && (currentFeature.centerY + j < height)) {
          if (temp.Pixel(currentFeature.centerX + i, currentFeature.centerY + j) == *redPixel) {
            window_isEmpty = false;
          }
        }
      }
    }
    if (window_isEmpty) {
      temp.SetPixel(currentFeature.centerX, currentFeature.centerY, *redPixel);
      recentLocations.push_back(currentFeature);
      counter++;
    }
    index++;
  }
}

double** R2Image::
frameProcessing(R2Image * otherImage)
{
  int px;
  int py;
  double sum;
  int min;
  int minX;
  int minY;
  vector<Feature> min_ssd;
  vector<Feature> vec_copy;

  for(int a = 0; a < recentLocations.size(); a++) {
    px = recentLocations[a].centerX;
    py = recentLocations[a].centerY;
    min = 100000;
    minX = 0;
    minY = 0;
      
 if ((px + (-0.01 * width)) >= 0 && (px + (-0.01 * width)) < width && (py + (-0.01 * height)) >= 0 && (py + (-0.01 * height)) < height) {
      for (int b = px + (-0.1 * width); b < px + (0.1 * width); b++) {
        for (int c = py + (-0.1 * height); c < py + (0.1 * height); c++) {
          R2Pixel *ssd = new R2Pixel();
          for (int u = -6; u < 6; u++) {
            for (int v = -6; v < 6; v++) {
              // calculate SSD for this region
              *ssd += (otherImage->Pixel(b+u, c+v) - Pixel(px+u,py+v)) * (otherImage->Pixel(b+u, c+v) - Pixel(px+u,py+v));
              sum = ssd->Red() + ssd->Green() + ssd->Blue();
            }
          }
          if (sum < min) {
            min = sum;
            minX = b;
            minY = c;
          }
        }
      }
      // add pixel with smallest SSD to vector
      vec_copy.push_back(Feature(px, py, 0));
      min_ssd.push_back(Feature(minX, minY, min));
    }
  }

  //mark tracked features
  for (int i = 0; i < min_ssd.size(); i++) {
      R2Pixel redPixel(1.0,0.0,0.0,1.0);
      if ((min_ssd[i].centerX - 5 > 0)
      && (min_ssd[i].centerX + 5 < width)
      && (min_ssd[i].centerY - 5 > 0)
      && (min_ssd[i].centerY + 5 < height)) {
        for(int u = -5; u < 5; u++){
          for(int v = -5; v < 5; v++){
            otherImage->SetPixel(min_ssd[i].centerX+u, min_ssd[i].centerY+v, redPixel);
          }
        }
      }
  }

  // run RANSAC and get optimal H transformation matrix
  int rand_index;

  int inliers;
  int N = 0;
  int max_inliers = 0;
  vector<Feature> vecInliers; 
  vector<Feature> vecMaxInliers; 

  int x1, x2, x3, x4;
  int y1, y2, y3, y4;
  int xd1, xd2, xd3, xd4;
  int yd1, yd2, yd3, yd4;

  double** matrixH;
  double** matrixH_optimal;
  vector<int> vectorA;
  vectorA.push_back(0);
  vectorA.push_back(0);
  vectorA.push_back(0);

  vector<double> vectorHA;
  vectorHA.push_back(0);
  vectorHA.push_back(0);

  double HA_x;
  double HA_y;
  double HA_z;
  double distance;

  // loop N (1000) number of times
  while (N <= 1000) {
    inliers = 0;

    // randomly choose 4 points
    rand_index = rand() % min_ssd.size();
    x1 = vec_copy[rand_index].centerX;
    y1 = vec_copy[rand_index].centerY;
    xd1 = min_ssd[rand_index].centerX;
    yd1 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x2 = vec_copy[rand_index].centerX;
    y2 = vec_copy[rand_index].centerY;
    xd2 = min_ssd[rand_index].centerX;
    yd2 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x3 = vec_copy[rand_index].centerX;
    y3 = vec_copy[rand_index].centerY;
    xd3 = min_ssd[rand_index].centerX;
    yd3 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x4 = vec_copy[rand_index].centerX;
    y4 = vec_copy[rand_index].centerY;
    xd4 = min_ssd[rand_index].centerX;
    yd4 = min_ssd[rand_index].centerY;

    vector<int> x;
    vector<int> y;
    vector<int> xd;
    vector<int> yd;

    x.push_back(x1);
    x.push_back(x2);
    x.push_back(x3);
    x.push_back(x4);

    y.push_back(y1);
    y.push_back(y2);
    y.push_back(y3);
    y.push_back(y4);

    xd.push_back(xd1);
    xd.push_back(xd2);
    xd.push_back(xd3);
    xd.push_back(xd4);

    yd.push_back(yd1);
    yd.push_back(yd2);
    yd.push_back(yd3);
    yd.push_back(yd4);

    // estimate and get H matrix
    matrixH = svdTest(x, y, xd, yd, 4);

    // loop over all A->B matches
    for (int i = 0; i < min_ssd.size(); i++) {
      // matrix multiplication to get HA
      vectorA[0] = vec_copy[i].centerX;
      vectorA[1] = vec_copy[i].centerY;
      vectorA[2] = 1;

      HA_x = matrixH[1][1]*vectorA[0] + matrixH[1][2]*vectorA[1] + matrixH[1][3]*vectorA[2];
      HA_y = matrixH[2][1]*vectorA[0] + matrixH[2][2]*vectorA[1] + matrixH[2][3]*vectorA[2];
      HA_z = matrixH[3][1]*vectorA[0] + matrixH[3][2]*vectorA[1] + matrixH[3][3]*vectorA[2];

      HA_x = HA_x/HA_z;
      HA_y = HA_y/HA_z;

      vectorHA[0] = HA_x;
      vectorHA[1] = HA_y;

      // find the distance between HA and B
      distance = sqrt(pow(vectorHA[0] - min_ssd[i].centerX, 2) + pow(vectorHA[1] - min_ssd[i].centerY, 2));

      // if distance is beneath threshold, increment counter
      if (distance <= 5.0) {
        inliers++;
        vecInliers.push_back(Feature(HA_x,HA_y,0));

      }
    }
    // keep track of H matrix with the most supporters
    if (inliers > max_inliers) {
      max_inliers = inliers;
      matrixH_optimal = matrixH;
      vecMaxInliers = vecInliers;
    }
    N++;
  }

  vector<int> optimal_x;
  vector<int> optimal_y;
  vector<int> optimal_xd;
  vector<int> optimal_yd;

  // re-calculate H matrix with ALL good points
  for (int i = 0; i < vecMaxInliers.size(); i++) {
    vectorA[0] = vec_copy[i].centerX;
    vectorA[1] = vec_copy[i].centerY;
    vectorA[2] = 1;

    HA_x = matrixH_optimal[1][1]*vectorA[0] + matrixH_optimal[1][2]*vectorA[1] + matrixH_optimal[1][3]*vectorA[2];
    HA_y = matrixH_optimal[2][1]*vectorA[0] + matrixH_optimal[2][2]*vectorA[1] + matrixH_optimal[2][3]*vectorA[2];
    HA_z = matrixH_optimal[3][1]*vectorA[0] + matrixH_optimal[3][2]*vectorA[1] + matrixH_optimal[3][3]*vectorA[2];

    HA_x = HA_x/HA_z;
    HA_y = HA_y/HA_z;

    vectorHA[0] = HA_x;
    vectorHA[1] = HA_y;

    distance = sqrt(pow(vectorHA[0] - min_ssd[i].centerX, 2) + pow(vectorHA[1] - min_ssd[i].centerY, 2));

    if (distance <= 5.0) {
      optimal_x.push_back(vec_copy[i].centerX);
      optimal_y.push_back(vec_copy[i].centerY);
      optimal_xd.push_back(min_ssd[i].centerX);
      optimal_yd.push_back(min_ssd[i].centerY);
    }
  }

  // draw lines
  for (int i = 0; i < optimal_x.size(); i++) {
    otherImage -> line(optimal_x[i], optimal_xd[i], optimal_y[i], optimal_yd[i], 0, 1, 0);
  }

  bestHMatrix = svdTest(optimal_x, optimal_y, optimal_xd, optimal_yd, max_inliers);

  // cout << "Optimal Number of Inliers: " << max_inliers << endl;

  
  return bestHMatrix;
}

void R2Image::
Sharpen()
{
  // make copy of image
  R2Image temp(*this);
  // initialize kernel
  const int SHARPEN_KERNEL[3][3]= {{0, -1, 0}, {-1, 5, -1}, {0, -1, 0}};
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j) = SHARPEN_KERNEL[0][1] * temp.Pixel(i-1,j)
      + SHARPEN_KERNEL[1][0] * temp.Pixel(i,j-1)
      + SHARPEN_KERNEL[1][1] * temp.Pixel(i,j)
      + SHARPEN_KERNEL[1][2] * temp.Pixel(i,j+1)
      + SHARPEN_KERNEL[2][1] * temp.Pixel(i+1,j);
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if (x0 > 3 && x0 < width - 3 && y0 > 3 && y0 < height - 3) {
    for (int x = x0 - 3; x <= x0 + 3; x++) {
      for (int y = y0 - 3; y <= y0 + 3; y++) {
        Pixel(x,y).Reset(r,g,b,1.0);
      }
    }
  }
  if (x0 > x1) {
    int x = y1;
    y1 = y0;
    y0 = x;

    x = x1;
    x1 = x0;
    x0 = x;
  }
  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;
  if (deltax != 0) deltaerr = fabs(float(float(deltay)/deltax));
  // Assume deltax != 0 (line is not vertical),
  // note that this division needs to be done in a way that preserves the fractional part
  int y = y0;
  for (int x = x0; x <= x1; x++) {
    Pixel(x,y).Reset(r,g,b,1.0);
    error = error + deltaerr;
    if (error >= 0.5) {
      if (deltay > 0) y = y + 1;
      else y = y - 1;
      error = error - 1.0;
    }
  }
}

double** R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
  // into this image with a 50% opacity.

  int px;
  int py;
  double sum;
  double min;
  int minX;
  int minY;
  vector<Feature> min_ssd;
  vector<Feature> vec_copy;

  for (int a = 0; a < recentLocations.size(); a++) {
    px = recentLocations[a].centerX;
    py = recentLocations[a].centerY;
    min = 1000000;
    minX = 0;
    minY = 0;

    if (px + (-0.1 * width) >= 0 && px + (-0.1 * width) < width && py + (-0.1 * height) >= 0 && py + (-0.1 * height) < height) {
      for (int b = px + (-0.1 * width); b < px + (0.1 * width); b++) {
        for (int c = py + (-0.1 * height); c < py + (0.1 * height); c++) {
          R2Pixel *ssd = new R2Pixel();
          for (int u = -3; u <= 3; u++) {
            for (int v = -3; v <= 3; v++) {
              // calculate SSD for this region
              *ssd += (otherImage->Pixel(b+u, c+v) - Pixel(px+u,py+v)) * (otherImage->Pixel(b+u, c+v) - Pixel(px+u,py+v));
              sum = ssd->Red() + ssd->Green() + ssd->Blue();
            }
          }

          if (sum < min) {
            min = sum;
            minX = b;
            minY = c;
          }
        }
      }
      // add pixel with smallest SSD to vector
      vec_copy.push_back(Feature(px, py, 0));
      min_ssd.push_back(Feature(minX, minY, min));
    }
  }

  int rand_index;

  int inliers;
  int N = 0;
  int max_inliers = 0;

  int x1, x2, x3, x4;
  int y1, y2, y3, y4;
  int xd1, xd2, xd3, xd4;
  int yd1, yd2, yd3, yd4;

  double** matrixH;
  double** matrixH_optimal;
  vector<int> vectorA;
  vectorA.push_back(0);
  vectorA.push_back(0);
  vectorA.push_back(0);

  vector<double> vectorHA;
  vectorHA.push_back(0);
  vectorHA.push_back(0);

  double HA_x;
  double HA_y;
  double HA_z;
  double distance;

  // loop N (1000) number of times
  while (N <= 1000) {
    inliers = 0;

    // randomly choose 4 points
    rand_index = rand() % min_ssd.size();
    x1 = vec_copy[rand_index].centerX;
    y1 = vec_copy[rand_index].centerY;
    xd1 = min_ssd[rand_index].centerX;
    yd1 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x2 = vec_copy[rand_index].centerX;
    y2 = vec_copy[rand_index].centerY;
    xd2 = min_ssd[rand_index].centerX;
    yd2 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x3 = vec_copy[rand_index].centerX;
    y3 = vec_copy[rand_index].centerY;
    xd3 = min_ssd[rand_index].centerX;
    yd3 = min_ssd[rand_index].centerY;

    rand_index = rand() % min_ssd.size();
    x4 = vec_copy[rand_index].centerX;
    y4 = vec_copy[rand_index].centerY;
    xd4 = min_ssd[rand_index].centerX;
    yd4 = min_ssd[rand_index].centerY;

    vector<int> x;
    vector<int> y;
    vector<int> xd;
    vector<int> yd;

    x.push_back(x1);
    x.push_back(x2);
    x.push_back(x3);
    x.push_back(x4);

    y.push_back(y1);
    y.push_back(y2);
    y.push_back(y3);
    y.push_back(y4);

    xd.push_back(xd1);
    xd.push_back(xd2);
    xd.push_back(xd3);
    xd.push_back(xd4);

    yd.push_back(yd1);
    yd.push_back(yd2);
    yd.push_back(yd3);
    yd.push_back(yd4);

    // estimate and get H matrix
    matrixH = svdTest(x, y, xd, yd, 4);

    // loop over all A->B matches
    for (int i = 0; i < min_ssd.size(); i++) {
      // matrix multiplication to get HA
      vectorA[0] = vec_copy[i].centerX;
      vectorA[1] = vec_copy[i].centerY;
      vectorA[2] = 1;

      HA_x = matrixH[1][1]*vectorA[0] + matrixH[1][2]*vectorA[1] + matrixH[1][3]*vectorA[2];
      HA_y = matrixH[2][1]*vectorA[0] + matrixH[2][2]*vectorA[1] + matrixH[2][3]*vectorA[2];
      HA_z = matrixH[3][1]*vectorA[0] + matrixH[3][2]*vectorA[1] + matrixH[3][3]*vectorA[2];

      HA_x = HA_x/HA_z;
      HA_y = HA_y/HA_z;

      vectorHA[0] = HA_x;
      vectorHA[1] = HA_y;

      // find the distance between HA and B
      distance = sqrt(pow(vectorHA[0] - min_ssd[i].centerX, 2) + pow(vectorHA[1] - min_ssd[i].centerY, 2));

      // if distance is beneath threshold, increment counter
      if (distance <= 4.0) {
        inliers++;
      }
    }
    // keep track of H matrix with the most supporters
    if (inliers > max_inliers) {
      max_inliers = inliers;
      matrixH_optimal = matrixH;
    }
    N++;
  }

  vector<int> optimal_x;
  vector<int> optimal_y;
  vector<int> optimal_xd;
  vector<int> optimal_yd;

  // re-calculate H matrix with ALL good points
  for (int i = 0; i < min_ssd.size(); i++) {
    vectorA[0] = vec_copy[i].centerX;
    vectorA[1] = vec_copy[i].centerY;
    vectorA[2] = 1;

    HA_x = matrixH_optimal[1][1]*vectorA[0] + matrixH_optimal[1][2]*vectorA[1] + matrixH_optimal[1][3]*vectorA[2];
    HA_y = matrixH_optimal[2][1]*vectorA[0] + matrixH_optimal[2][2]*vectorA[1] + matrixH_optimal[2][3]*vectorA[2];
    HA_z = matrixH_optimal[3][1]*vectorA[0] + matrixH_optimal[3][2]*vectorA[1] + matrixH_optimal[3][3]*vectorA[2];

    HA_x = HA_x/HA_z;
    HA_y = HA_y/HA_z;

    vectorHA[0] = HA_x;
    vectorHA[1] = HA_y;

    distance = sqrt(pow(vectorHA[0] - min_ssd[i].centerX, 2) + pow(vectorHA[1] - min_ssd[i].centerY, 2));

    if (distance <= 4.0) {
      optimal_x.push_back(vec_copy[i].centerX);
      optimal_y.push_back(vec_copy[i].centerY);
      optimal_xd.push_back(min_ssd[i].centerX);
      optimal_yd.push_back(min_ssd[i].centerY);
    }
  }

  // draw lines
  for (int i = 0; i < optimal_x.size(); i++) {
    line(optimal_x[i], optimal_xd[i], optimal_y[i], optimal_yd[i], 0, 1, 0);
  }

  matrixH_optimal = svdTest(optimal_x, optimal_y, optimal_xd, optimal_yd, max_inliers);

  cout << "Optimal Number of Inliers: " << max_inliers << endl;

  cout << "\nMatrix H: " << endl;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      cout << matrixH_optimal[i][j] << "\t";
    }
    cout << "\n";
  }

  return matrixH_optimal;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage, double** bestHMatrix)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.


  // make copy of image
  R2Image temp(*this);

  double x;
  double y;
  double z;
  vector<int> xd;
  vector<int> yd;
  R2Pixel imageB_pixel;
  R2Pixel blended_pixel;
  R2Pixel *white = new R2Pixel(1.0,1.0,1.0,1.0);

  double determinant = 0; 
  double** minorMatrix = dmatrix(1,3,1,3);
  

  cout << "\nMatrix H: " << endl;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      cout << bestHMatrix[i][j] << "\t";
    }
    cout << "\n";
  }

  minorMatrix[1][1]=(bestHMatrix[2][2]*bestHMatrix[3][3])-(bestHMatrix[2][3]*bestHMatrix[3][2]);                //***********************************
  minorMatrix[1][2]=-1*((bestHMatrix[2][1]*bestHMatrix[3][3])-(bestHMatrix[2][3]*bestHMatrix[3][1]));        //                                                                 *
  minorMatrix[1][3]=(bestHMatrix[2][1]*bestHMatrix[3][2])-(bestHMatrix[2][2]*bestHMatrix[3][1]);                //                                                                 *
  minorMatrix[2][1]=-1*((bestHMatrix[1][2]*bestHMatrix[3][3])-(bestHMatrix[3][2]*bestHMatrix[1][3]));        //                                                                 *
  minorMatrix[2][2]=(bestHMatrix[1][1]*bestHMatrix[3][3])-(bestHMatrix[1][3]*bestHMatrix[3][1]);                // CALCULATION OF COFACTOR          *
  minorMatrix[2][3]=-1*((bestHMatrix[1][1]*bestHMatrix[3][2])-(bestHMatrix[1][2]*bestHMatrix[3][1]));        //                                                                 *
  minorMatrix[3][1]=bestHMatrix[1][2]*bestHMatrix[2][3]-bestHMatrix[1][3]*bestHMatrix[2][2];                //                                                                 *
  minorMatrix[3][2]=-1*((bestHMatrix[1][1]*bestHMatrix[2][3])-(bestHMatrix[1][3]*bestHMatrix[2][1]));        //                                                                 *
  minorMatrix[3][3]=(bestHMatrix[1][1]*bestHMatrix[2][2])-(bestHMatrix[1][2]*bestHMatrix[2][1]);                //***********************************


  determinant = (bestHMatrix[1][1]*minorMatrix[1][1])
               -(bestHMatrix[1][2]*minorMatrix[1][2])
               +(bestHMatrix[1][3]*minorMatrix[1][3]);
  
  cout<< "Determinant: " << determinant << endl;

  double** adjointMatrix = dmatrix(1,3,1,3);
  for(int i = 1 ;  i <= 3; i++){
    for(int j = 1; j <= 3; j++){
      adjointMatrix[i][j] = minorMatrix[j][i];
    }
  }

  double** inverseMatrix = dmatrix(1,3,1,3);
   for(int i = 1 ;  i <= 3; i++){
    for(int j = 1; j <= 3; j++){
      inverseMatrix[i][j] = adjointMatrix[i][j] * (1/determinant);
    }
  }

  cout << "\nInverse Matrix H: " << endl;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      cout << inverseMatrix[i][j] << "\t";
    }
    cout << "\n";
  }

  for(double i = 0; i < otherImage -> width; i++){
    for(double j=0; j < otherImage -> height; j++){
       x = inverseMatrix[1][1]*i + inverseMatrix[1][2]*j + inverseMatrix[1][3];
       y = inverseMatrix[2][1]*i + inverseMatrix[2][2]*j + inverseMatrix[2][3];
       z = inverseMatrix[3][1]*i + inverseMatrix[3][2]*j + inverseMatrix[3][3];

       x = x/z; 
       y = y/z; 

       if ((x > 0 && x < width) && (y > 0 && y < height)){
         imageB_pixel = otherImage -> Pixel(i,j);
         blended_pixel = (imageB_pixel + temp.Pixel(x,y))/2;

         SetPixel(x,y,blended_pixel);

       }
    }
  }

  

  // apply H matrix
  // for (double i = 0; i < width; i++) {
  //   for (double j = 0; j < height; j++) {


  //     x = bestHMatrix[1][1]*i + bestHMatrix[1][2]*i + bestHMatrix[1][3];
  //     y = bestHMatrix[2][1]*i + bestHMatrix[2][2]*j + bestHMatrix[2][3];
  //     z = bestHMatrix[3][1]*i + bestHMatrix[3][2]*j + bestHMatrix[3][3];


  //     x = x/z;
  //     y = y/z;
      
  //     if ((x > 0 && x < width) && (y > 0 && y < height)) {

  //       imageB_pixel = otherImage->Pixel(i,j);
  //       blended_pixel = (imageB_pixel + temp.Pixel(x,y)) / 2;
  //       SetPixel(i,j,blended_pixel);

  //     }
  //     // } else {
  //     //   SetPixel(i,j,*white);
  //     // }
  //   }
  // }


}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        // Read pixel values
        int red, green, blue;
        if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
          fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
          fclose(fp);
          return 0;
        }

        // Assign pixel values
        double r = (double) red / max_value;
        double g = (double) green / max_value;
        double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
extern "C" {
  #   define XMD_H // Otherwise, a conflict with INT32
  #   undef FAR // Otherwise, a conflict with windows.h
  #   include "jpeg/jpeglib.h"
};
#endif



int R2Image::
ReadJPEG(const char *filename)
{
  #ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
  #else
  fprintf(stderr, "JPEG not supported");
  return 0;
  #endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
  #ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
  #else
  fprintf(stderr, "JPEG not supported");
  return 0;
  #endif
}
