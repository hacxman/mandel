#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

typedef struct {
  double r;
  double g;
  double b;
} color_dt;

void writeimage(char * filename, color_dt *color, int depth, int w, int h) {
  FILE * fout = fopen(filename, "w+");
  double maxval = (2 << (depth - 1)) - 1;
  fprintf(fout, "P3\n%d %d\n%d\n", w, h, (int)maxval);
  int i,j;
  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      int r = color[i*h + j].r * maxval;
      int g = color[i*h + j].g * maxval;
      int b = color[i*h + j].b * maxval;
      r = r > maxval ? maxval : (r < 0 ? 0 : r);
      g = g > maxval ? maxval : (g < 0 ? 0 : g);
      b = b > maxval ? maxval : (b < 0 ? 0 : b);
      fprintf(fout, "%d %d %d ", r, g, b);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}

inline int iterate(complex double z, complex double c, int mi, complex double *_z) {
  int it = 0;
  for (it = 0 ; it < mi; it++) {
    if (creal(z)*creal(z) + cimag(z)*cimag(z) < 2.0f*2.0f) break;
    z = z*z + c;
//    it++;
//    if (it >= mi) {
//      *_z = z;
//      return it;
//    }
  }
  *_z = z;
  return it;
}

inline complex double xytoC(int x, int y, int w, int h, int cx, int cy, 
                     double mx, double my, double cdx, double cdy, double zoom) {

  double a = cdx - mx*2*((x - cx)/zoom) / (double)w;
  double b = cdy - my*2*((y - cy)/zoom) / (double)h;
  return a + b*_Complex_I;

}

int main(void) {
  color_dt * img;
  int w = 10240, h = 10240;
//  int w = 2048, h = 2048;
  img = malloc(sizeof(color_dt) * w * h);
  int i,j;
  int maxit = 400;
  for (i = 0; i < w; i++) {
//    printf("%d/%d      \r", i, w);
    //fflush(stdout);
  #pragma omp parallel for
    for (j = 0; j < h; j++) {
      complex double _z;
      int r = iterate(0+0i, xytoC(i,j,w,h,w/2,h/2,2.0,2.0,0.00f,0.0f,1.0f), maxit, &_z);
      double lg = log2(log10(cabs(_z))/log10(2.0));
      double vz = (r - lg) / (double)maxit;
      img[i*h+j].r = sin(maxit/10*vz)*sin(maxit/10*vz);
      img[i*h+j].g = cos(1.5*maxit/10*vz)*cos(1.5*maxit/10*vz); //1.0f - r/(double)maxit;
      img[i*h+j].b = sin(1.7*maxit/10*vz)*sin(1.7*maxit/10*vz);
    }
  }
  //writeimage("test.pnm", img, 16, w, h);
}

