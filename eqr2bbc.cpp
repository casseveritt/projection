// equirect to baseball-cover
// Cass Everitt
// Dec 14, 2015


#include <stdio.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "envmap.h"

int main(int argc, char **argv) {

  printf( "equirect to baseball-cover conversion\n" );

  int w = 0;
  int h = 0;
  int comp = 0;

  if ( argc != 3 ) {
  	printf( "usage: eqr2bbc <eqr-image> <bbc-edge-length>\n");
  	return 1;
  }
  int edgeLength = atof( argv[2] );

  unsigned char * eqr = stbi_load( argv[1], &w, &h, &comp, 0 );

  printf( "width = %d, height = %d, num_components = %d\n", w, h, comp );

  int dw = edgeLength * 3;
  int dh = edgeLength * 2;
  unsigned char * r = new unsigned char[ dw * dh * 4 ];

  for( int y = 0; y < dh; y++ ) {
    unsigned char * line = r + (dh - 1 - y) * dw * 4;
    Vec2 p;
    p.y = ( y + 0.5f ) / dh;
    for( int x = 0; x < dw; x++ ) {
      p.x = ( x + 0.5f ) / dw;
      Vec3 dir = baseball_cover_index_to_direction( p );
      Vec2 eqrtc = direction_to_equirect( dir );
      sample_bilerp( eqr, w, h, comp, eqrtc.x, 1.0f - eqrtc.y, line );
      line[3] = 255;
      line += 4;
    }
  }

  stbi_write_png( "bbc.png", dw, dh, 4, r, 0);
  delete [] r;

  stbi_image_free( eqr );

  return 0;
}
