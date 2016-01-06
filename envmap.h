// Cass Everitt - Dec 30, 2014

#ifndef ENVMAP_H_
#define ENVMAP_H_

#include <cmath>
#include <stdio.h>
#include <algorithm>

using namespace std;

#define PI 3.1415927

struct CubeIndex {
  CubeIndex(){}
  CubeIndex( int _face, double _x, double _y ) : face(_face), x(_x), y(_y) {}
  int face;
  double x, y;
};

struct Vec2 {
  union {
    double v[2];
    struct { double x, y; };
  };
};

Vec2 make_Vec2( double x, double y ) {
  Vec2 v;
  v.x = x; v.y = y;
  return v;
}

struct Vec3 {
  union {
    double v[3];
    struct { double x, y, z; };
  };
};

Vec3 make_Vec3( double x, double y, double z ) {
  Vec3 v;
  v.x = x; v.y = y; v.z = z;
  return v;
}

Vec3 make_Vec3( const double * vec ) {
  Vec3 v;
  v.x = vec[0]; v.y = vec[1]; v.z = vec[2];
  return v;
}



template <typename T> struct vector_traits { static const int num_elements = 0; };
template <> struct vector_traits<Vec3> { static const int num_elements = 3; };
template <> struct vector_traits<Vec2> { static const int num_elements = 2; };

template <typename T> T operator + ( const T & v, const T & w ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] + w.v[n];
  }
  return o;
}

template <typename T> T operator - ( const T & v, const T & w ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] - w.v[n];
  }
  return o;
}

template <typename T> T operator + ( const T & v, double f ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] - f;
  }
  return o;
}

template <typename T> T operator - ( const T & v, double f ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] - f;
  }
  return o;
}

template <typename T> T operator * ( const T & v, double f ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] * f;
  }
  return o;
}

template <typename T> T operator * ( const T & v, const T & u ) {
  Vec3 o;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] * u.v[n];
  }
  return o;
}

template <typename T> T operator / ( const T & v, double f ) {
  Vec3 o = v;
  for(int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    o.v[n] = v.v[n] / f;
  }
  return o;
}

template <typename T> double dot( const T & a, const T & b ) {
  double sum = 0.0;
  for( int n = 0; n < vector_traits<T>::num_elements; n++ ) {
    sum += a.v[n] * b.v[n];
  }
  return sum;
}

template <typename T> double length( const T & v ) {
  return sqrt( dot( v, v ) );
}

template <typename T> T normalize( const T & v ) {
  double len = length( v );
  T n = v / len;
  return n;
}

Vec3 cross( const Vec3 & a, const Vec3 & b ) {
  Vec3 c;
  c.x = a.y * b.z - a.z * b.y;
  c.y = a.z * b.x - a.x * b.z;
  c.z = a.x * b.y - a.y * b.x;
  return c;
}

double deg2rad( double deg ) {
  return (PI / 180.0) * deg;
}

double rad2deg( double rad ) {
  return (180.0 / PI) * rad;
}

Vec3 equirect_to_direction( const Vec2 & e ) { 
  double lon = PI * e.x;
  double lat = (PI / 2.0) * e.y;
  Vec3 v;
  v.x = cos( lat ) * sin( lon );  
  v.y = sin( lat );
  v.z = - cos( lat ) * cos( lon );
  return v;
}

Vec2 direction_to_equirect( const Vec3 & v ) { 
  Vec2 p;
  double xz = sqrt( v.x * v.x + v.z * v.z );
  p.x = atan2( v.x, v.z ) / ( 2.0 * PI );
  p.y = ( atan2( v.y, xz ) / PI ) + 0.5;
  if( p.x < 0.0 ) {
    p.x += 1.0;
  }
  p.x = 1.0 - p.x;
  return p;
}

// remaps [-1,1] to itself with a different sample distribution
#if 0
double unpinch( double x ) {
  double ax = abs(x);
  double y = 1.25 * x - 0.25 * x * ax;
  return y;
}

double pinch( double x ) {
  double y = 2.5 - 0.5 * sqrt( 25.0 - 16.0 * abs(x) );
  return x > 0 ? y : -y;
}
#endif
#if 0
double unpinch( double x ) {
  double ax = abs(x);
  double y = 1.5 * x - 0.5 * x * ax;
  return y;
}

double pinch( double x ) {
  double y = 1.5 - 0.5 * sqrt( 9.0 - 8.0 * abs(x) );
  return x > 0 ? y : -y;
}
#endif
#if 1
double unpinch( double x ) {
  double ax = abs(x);
  double y = 1.375 * x - 0.375 * x * ax;
  return y;
}

double pinch( double x ) {
  double y = 1.833333 - 0.16666667 * sqrt( 121.0 - 96.0 * abs(x) );
  return x > 0 ? y : -y;
}
#endif

Vec3 cube_index_to_direction( const CubeIndex & ci ) {
  Vec3 v;
  v.x = ci.x; v.y = ci.y; v.z = 1.0;
  v = normalize( v );
  Vec3 swz;
  switch( ci.face ) {
    case 0: // +Y
      swz.z =  v.x;
      swz.x = -v.y;
      swz.y =  v.z;
      break;
    case 1: // +Z
      swz.y = -v.x;
      swz.x = -v.y;
      swz.z =  v.z;
      break;
    case 2: // -Y
      swz.z = -v.x;
      swz.x = -v.y;
      swz.y = -v.z;
      break;
    case 3: // -X
      swz.z = -v.x;
      swz.y =  v.y;
      swz.x = -v.z;
      break;
    case 4: // -Z
      swz.x =  v.x;
      swz.y =  v.y;
      swz.z = -v.z;
      break;
    case 5: // +X
      swz.z =  v.x;
      swz.y =  v.y;
      swz.x =  v.z;
      break;
    default:
      break;
  }
  return swz;
}

Vec3 unpinched_cube_index_to_direction( CubeIndex ci ) {
  ci.x = pinch( ci.x );
  ci.y = pinch( ci.y );
  return cube_index_to_direction( ci );
}

Vec3 pinched_cube_index_to_direction( CubeIndex ci ) {
  ci.x = unpinch( ci.x );
  ci.y = unpinch( ci.y );
  return cube_index_to_direction( ci );
}

Vec2 cube_index_to_baseball_cover_index( const CubeIndex & ci ) {
  Vec2 v;
  v.x = (ci.x + 1.0) / 6.0;
  v.y = (ci.y + 1.0) / 4.0;
  if( ci.face > 2 ) {
    v.y += 0.5f;
  }
  switch( ci.face ) {
    case 1:
    case 4:
      v.x += 1.0/3.0;
      break;
    case 2:
    case 5:
      v.x += 2.0/3.0;
      break;
    default:
      break;
  }
  return v;
}

CubeIndex baseball_cover_index_to_cube_index( const Vec2 & bbci ) {
  CubeIndex ci;

  int col = floor(bbci.x * 3.0);
  int row = floor(bbci.y * 2.0);
  ci.face = 3 * row + col;
  ci.x = (bbci.x * 3.0 - col) * 2.0 - 1.0;
  ci.y = (bbci.y * 2.0 - row) * 2.0 - 1.0;
  return ci;
}

Vec3 unpinched_baseball_cover_index_to_direction( const Vec2 & bbci ) {
  CubeIndex ci = baseball_cover_index_to_cube_index( bbci );
  Vec3 dir = unpinched_cube_index_to_direction( ci );
  return dir;
}

Vec3 pinched_baseball_cover_index_to_direction( const Vec2 & bbci ) {
  CubeIndex ci = baseball_cover_index_to_cube_index( bbci );
  Vec3 dir = pinched_cube_index_to_direction( ci );
  return dir;
}

Vec3 baseball_cover_index_to_direction( const Vec2 & bbci ) {
  CubeIndex ci = baseball_cover_index_to_cube_index( bbci );
  Vec3 dir = cube_index_to_direction( ci );
  return dir;
}

CubeIndex unpinch_cube_index( const CubeIndex & ci ) {
  CubeIndex rci;
  rci.face = ci.face;
  rci.x = unpinch(ci.x);
  rci.y = unpinch(ci.y);
  return rci;
}

CubeIndex pinch_cube_index( const CubeIndex & ci ) {
  CubeIndex rci;
  rci.face = ci.face;
  rci.x = pinch(ci.x);
  rci.y = pinch(ci.y);
  return rci;
}

CubeIndex direction_to_cube_index( const Vec3 &v ) {
  CubeIndex ci;

  Vec3 av;
  av.x = abs(v.x);
  av.y = abs(v.y);
  av.z = abs(v.z);

  int fc =  0;
  fc |= (av.x > av.y ? 1 : 0 ) << 0;
  fc |= (av.x > av.z ? 1 : 0 ) << 1;
  fc |= (av.y > av.z ? 1 : 0 ) << 2;
  
  switch(fc) {
    case 0: // Z
    case 1:
      if( v.z > 0 ) {
        ci.face = 1;
        ci.x = -v.y/av.z;
        ci.y = -v.x/av.z;
      } else {
        ci.face = 4;
        ci.x = v.x/av.z;
        ci.y = v.y/av.z;
      }
      break;

    case 3: // X
    case 7: 
      if( v.x > 0 ) {
        ci.face = 5;
        ci.x = v.z/av.x;
        ci.y = v.y/av.x;
      } else {
        ci.face = 3;
        ci.x = -v.z/av.x;
        ci.y =  v.y/av.x;
      }
      break;
    
    case 4: // Y
    case 6:
      if( v.y > 0 ) {
        ci.face = 0;
        ci.x =  v.z/av.y;
        ci.y = -v.x/av.y;
      } else {
        ci.face = 2;
        ci.x = -v.z/av.y;
        ci.y = -v.x/av.y;
      }

    case 2: // impossible: x <= y <= z incompatible with x > z
    case 5: // impossible: x  > y  > z incompatible with x <= z
    default:
      break;
  }
  return ci;
}

float lerp( float a, float b, float f ) {
  return a * (1.0f - f) + b * f;
}

float bilerp( float p00, float p10, float p01, float p11, float fx, float fy ) {
  float y0 = lerp( p00, p10, fx );
  float y1 = lerp( p01, p11, fx );
  return lerp( y0, y1, fy );
}

void sample_bilerp( const unsigned char * img, int width, int height, int comp, float s, float t, unsigned char *dst ) {
   float fi = s * width + 0.5f;
   float fj = t * height + 0.5f;
   int i = fi - 1.0f;
   int j = fj - 1.0f;
   int i0 = max( i, 0 );
   int i1 = min( i + 1, width - 1 );
   int j0 = max( j, 0 );
   int j1 = min( j + 1, height - 1 );
   float ifrac = fi - floor( fi );
   float jfrac = fj - floor( fj );
   const unsigned char * p00 = img + ( ( j0 * width + i0 ) * comp );
   const unsigned char * p10 = img + ( ( j0 * width + i1 ) * comp );
   const unsigned char * p01 = img + ( ( j1 * width + i0 ) * comp );
   const unsigned char * p11 = img + ( ( j1 * width + i1 ) * comp );
   for( int c = 0; c < comp; c++ ) {
      dst[c] = (unsigned char)bilerp( p00[c], p10[c], p01[c], p11[c], ifrac, jfrac );
   }
}

double raysphere( Vec3 c, Vec3 ray ) {
  Vec3 o = make_Vec3( 0.0, 0.0, 0.0 );

  ray = normalize( ray );

  Vec3 ominusc = o - c;

  double loc = dot( ray, ominusc );

  // radius == 1.0
  double root = loc * loc - dot( ominusc, ominusc ) + 1.0;

  if( root <= 0.0 ) {
    return 0.0;
  }

  root = sqrt( root );

  if( root < loc ) {
    return 0.0;
  }

  return root - loc;
}

#endif // ENVMAP_H_


