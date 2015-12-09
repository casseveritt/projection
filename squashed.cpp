// squashed
// Cass Everitt
// Dec 7, 2015


#include <stdio.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "envmap.h"
#include "linear.h"

void resample_fwd( Vec3 disp, double squash, unsigned char * src, int srcEdgeLength, int comp, unsigned char * dst, int dstEdgeLength ) {

  int w = srcEdgeLength * 3;
  int h = srcEdgeLength * 2;

  int dw = dstEdgeLength * 3;
  int dh = dstEdgeLength * 2;

  disp = normalize( disp );

  r3::Vec3d d( &disp.x );
  r3::Vec3d z( 0.0f, 1.0f, 0.0f );
  r3::Rotationd fwd( d, z );
  r3::Rotationd inv = fwd.Inverse();

  for( int y = 0; y < dh; y++ ) {
    unsigned char * line = dst + (dh - 1 - y) * dw * 4;
    Vec2 p;
    p.y = ( y + 0.5f ) / dh;
    for( int x = 0; x < dw; x++ ) {
      p.x = ( x + 0.5f ) / dw;
      Vec3 dir = baseball_cover_index_to_direction( p );
      r3::Vec3d d2( & dir.x );
      fwd.MultVec( d2 );
      d2 *= r3::Vec3d( squash, 1.0, squash );
      inv.MultVec( d2 );
      Vec3 v = make_Vec3( d2.Ptr() );
      CubeIndex ci = direction_to_cube_index( v );
      Vec2 bbci = cube_index_to_baseball_cover_index( ci );
      sample_bilerp( src, w, h, comp, bbci.x, 1.0f - bbci.y, line );
      line[3] = 255;
      line += 4;
    }
  }

}

void resample_inv( Vec3 disp, double squash, unsigned char * src, int srcEdgeLength, int comp, unsigned char * dst, int dstEdgeLength ) {

  int w = srcEdgeLength * 3;
  int h = srcEdgeLength * 2;

  int dw = dstEdgeLength * 3;
  int dh = dstEdgeLength * 2;

  disp = normalize( disp );

  r3::Vec3d d( &disp.x );
  r3::Vec3d z( 0.0f, 1.0f, 0.0f );
  r3::Rotationd fwd( d, z );
  r3::Rotationd inv = fwd.Inverse();

  for( int y = 0; y < dh; y++ ) {
    unsigned char * line = dst + (dh - 1 - y) * dw * 4;
    Vec2 p;
    p.y = ( y + 0.5f ) / dh;
    for( int x = 0; x < dw; x++ ) {
      p.x = ( x + 0.5f ) / dw;
      Vec3 dir = baseball_cover_index_to_direction( p );
      r3::Vec3d d2( & dir.x );
      fwd.MultVec( d2 );
      d2 *= r3::Vec3d( 1.0/squash, 1.0, 1.0/squash );
      inv.MultVec( d2 );
      Vec3 v = make_Vec3( d2.Ptr() );
      CubeIndex ci = direction_to_cube_index( v );
      Vec2 bbci = cube_index_to_baseball_cover_index( ci );
      sample_bilerp( src, w, h, comp, bbci.x, 1.0f - bbci.y, line );
      line[3] = 255;
      line += 4;
    }
  }


}


int main(int argc, char **argv) {

  printf( "recenter projection\n" );

  int w = 0;
  int h = 0;
  int comp = 0;

  unsigned char * bbc = stbi_load( argv[1], &w, &h, &comp, 0 );

  printf( "width = %d, height = %d, num_components = %d\n", w, h, comp );

  int dw = w;
  int dh = h;


  Vec3 disp;
  double lat, lon, r, squash;

  if( argc == 6 ) {
	lat = atof( argv[2] );
  	lon = atof( argv[3] );
  	r = atof( argv[4] );
  	squash = atof( argv[5] );
  } else { // pick Miami
  	lat = 25;
  	lon = -80.21;
  	r = .8;
  	squash = 0.25;
  }

  disp.x = sin( deg2rad( lon ) ) * cos( deg2rad( -lat ) );
  disp.z = cos( deg2rad( lon ) ) * cos( deg2rad( -lat ) );
  disp.y = sin( deg2rad( -lat ) );

  disp = disp * r;


  printf( "Displacement: ( %lf, %lf, %lf )\n", disp.x, disp.y, disp.z );

  unsigned char * rfwd = new unsigned char[ dw * dh * 4 ];
  unsigned char * rinv = new unsigned char[ dw * dh * 4 ];

  resample_inv( disp, squash, bbc, h / 2, comp, rinv, dh / 2 );  
  resample_fwd( disp, squash, rinv, dh / 2, 4, rfwd, dh / 2 );

  stbi_image_free( bbc );

  printf( "writing squash_inv.png\n" );
  stbi_write_png( "squash_inv.png", dw, dh, 4, rinv, 0);
  printf( "writing squash_fwd.png\n" );
  stbi_write_png( "squash_fwd.png", dw, dh, 4, rfwd, 0);

  delete [] rfwd;
  delete [] rinv;

  return 0;
}
