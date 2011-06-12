#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <xmmintrin.h>

#include <assert.h>

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

const int count = 16;
complex double * kokot(complex double z[count], complex double c[count]) {
  int i;
  for (i = 0; i < count; i++) {
    z[0] = z[0]*z[0] + c[0];
  }
  return z;
}

void test_fill() {
  int w = 1024, h = 1024;
  color_dt * img = malloc(sizeof(color_dt) * w * h);
  int i,j;
  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      img[i*h+j].r = (double)i / (double)w; //0.3;
      img[i*h+j].g = 0;
      img[i*h+j].b = 0;
    }
  }

  writeimage("test.pnm", img, 16, w, h);
}

typedef double d2 __attribute__ ((vector_size (16)));
typedef union {double v[2]; d2 reg;} d2_u;
typedef struct {int v[2];} i2;

inline i2 iterate(d2 zre, d2 zim, d2 cre, d2 cim, int mi, d2_u * _zre, d2_u * _zim) { //(complex double z, complex double c, int mi, complex double *_z) {
  int it = 0;
  d2 nzre, nzim;
  d2 end_c;
  i2 r; r.v[0] = -1; r.v[1] = -1;
  for (it = 0 ; it < mi; it++) {
//    z = z*z + c;
    nzre = zre * zre - zim * zim;
    nzim = zim * zre + zre * zim;
    nzre += cre;
    nzim += cim;
    zre = nzre;
    zim = nzim;
    
    end_c = zre * zre + zim * zim;
  
    int tmp; 
    tmp = 0; 
    if ((((d2_u)end_c).v[tmp] > 4.0) && (r.v[tmp] < 0)) {
      _zre->v[tmp] = ((d2_u)zre).v[tmp];
      _zim->v[tmp] = ((d2_u)zim).v[tmp];
      r.v[tmp] = it;
    }// else {
      //_zre->v[tmp] = _zre->v[tmp];
      //_zim->v[tmp] = _zim->v[tmp];
      //r.v[tmp] = r.v[tmp];
    //};

    tmp = 1;
    if ((((d2_u)end_c).v[tmp] > 4.0) && (r.v[tmp] < 0)) {
      _zre->v[tmp] = ((d2_u)zre).v[tmp];
      _zim->v[tmp] = ((d2_u)zim).v[tmp];
      r.v[tmp] = it;
    }// else {
    //  _zre->v[tmp] = _zre->v[tmp];
    //  _zim->v[tmp] = _zim->v[tmp];
    //  r.v[tmp] = r.v[tmp];
    //};

/*
    _zre->v[1] = 
        (((d2_u)end_c).v[1] > 4.0) && (r.v[1] < 0) ? 
          ((d2_u)zre).v[1] : 
          _zre->v[1];
    _zim->v[1] = 
        (((d2_u)end_c).v[1] > 4.0) && (r.v[1] < 0) ? 
          ((d2_u)zim).v[1] : 
          _zim->v[1];
    r.v[1] = (((d2_u)end_c).v[1] > 4.0) && (r.v[1] < 0) ? it : r.v[1];*/

    if ( (r.v[0] != -1) && (r.v[1] != -1) ) {
//      printf("%lf %lf\n", ((d2_u)end_c).v[0], ((d2_u)end_c).v[1]);
      break;
    }
  }
//  *_zre = zre;
//  *_zim = zim;
  return r;
}


d2 xytoC (d2 xy, d2 wh, d2 cxy, d2 mxy, d2 cdxy, double zoom) {  /*(int x, int y, int w, int h, int cx, int cy, 
                     double mx, double my, double cdx, double cdy, double zoom) {*/

//  double a = cdx - mx*2*((x - cx)/zoom) / (double)w;
//  double b = cdy - my*2*((y - cy)/zoom) / (double)h;
//  return a + b*_Complex_I;
  d2_u zoom_v; 
  zoom_v.v[0] = zoom;
  zoom_v.v[1] = zoom;
  
  d2_u two;
  two.v[0] = 2.0f;
  two.v[1] = 2.0f;

  return cdxy - mxy*two.reg*((xy - cxy)/zoom_v.reg) / wh;
}

int main(void) {
  color_dt * img;
  int w = 10240, h = 10240;
//  int w = 2048, h = 2048;
  img = malloc(sizeof(color_dt) * w * h);
  assert(img);
  int i,j;
  int maxit = 400;
  int * its = malloc(sizeof(int) * w * h);
  d2_u * out_zs = malloc(sizeof(d2) * w * h);
  assert(its);
  assert(out_zs);
        d2_u zero; zero.v[0] = 0.0f; zero.v[1] = 0.0f;
        d2_u two;  two.v[0] = 2.0f; two.v[1] = 2.0f;

  for (i = 0; i < w; i++) {
//    printf("%d/%d      \r", i, w);
    //fflush(stdout);
    #pragma omp parallel for
    for (j = 0; j < h; j+=2) {
        d2_u _zre;
        d2_u _zim;
        d2_u xy; xy.v[0] = i; xy.v[1] = j;
        d2_u wh; wh.v[0] = w; wh.v[1] = h;
        d2 coord1 = xytoC(xy.reg, wh.reg, wh.reg/two.reg, two.reg, zero.reg, 1.0f);
        //printf("c1.x=%lf c1.y=%lf\n",((d2_u)coord1).v[0], ((d2_u)coord1).v[1]);

        xy.v[0] = i; xy.v[1] = j+1;
        d2 coord2 = xytoC(xy.reg, wh.reg, wh.reg/two.reg, two.reg, zero.reg, 1.0f);

        d2_u re; 
          re.v[0] = ((d2_u)coord1).v[0];
          re.v[1] = ((d2_u)coord2).v[0];
        d2_u im; 
          im.v[0] = ((d2_u)coord1).v[1];
          im.v[1] = ((d2_u)coord2).v[1];


        i2 r = iterate(zero.reg, zero.reg, re.reg, im.reg , maxit, &_zre, &_zim);
        // * ((i2*)(its + i*h+j)) = r;
        img[i*h+j].r = _zre.v[0];
        img[i*h+j].g = _zim.v[0];
        img[i*h+j+1].r = _zre.v[1];
        img[i*h+j+1].g = _zim.v[1];

        img[i*h+j].b = r.v[0];
        img[i*h+j+1].b = r.v[1];
/*
        out_zs[i*h+j].v[0] = _zre.v[0];
        out_zs[i*h+j].v[1] = _zim.v[0];
        out_zs[i*h+j+1].v[0] = _zre.v[1];
        out_zs[i*h+j+1].v[1] = _zim.v[1];*/
//      complex double _z;
//      int r = iterate(0+0i, xytoC(i,j,w,h,w/2,h/2,2.0,2.0,0.00f,0.0f,1.0f), maxit, &_z);
    }
  };
  for (i = 0; i < w; i++) {
    #pragma omp parallel for
    for (j = 0; j < h; j++) {
        //double lg = log2(log10(cabs(_z))/log2(2.0));
        complex double d = img[i*h+j].r + img[i*h+j].g*_Complex_I;
        double lg = log2(log10(cabs(d))/log10(2.0));
        int r = img[i*h+j].b; //its[i*h+j];
//        printf("%lf\n", lg);
        double vz = (r - lg) / (double)maxit;
        img[i*h+j].r = sin(maxit/10*vz)*sin(maxit/10*vz);
        img[i*h+j].g = cos(1.5*maxit/10*vz)*cos(1.5*maxit/10*vz); //1.0f - r/(double)maxit;
        img[i*h+j].b = sin(1.7*maxit/10*vz)*sin(1.7*maxit/10*vz);
    }
  }
/*      img[i*h+j].r = sin(maxit/10*vz)*sin(maxit/10*vz);
      img[i*h+j].g = cos(1.5*maxit/10*vz)*cos(1.5*maxit/10*vz); //1.0f - r/(double)maxit;
      img[i*h+j].b = sin(1.7*maxit/10*vz)*sin(1.7*maxit/10*vz);
    }
  }*/
//  writeimage("test.pnm", img, 16, w, h);
}

