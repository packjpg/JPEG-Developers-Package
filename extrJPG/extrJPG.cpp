#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VER	"0.1"
// #define LFS
#define PAUSE

void dump_JPEGs( char* filename )
{
	FILE* fp_in = NULL; // file pointer
	FILE* fp_out = NULL; // file pointer
	fpos_t jpg0; // start adress of JPG file
	fpos_t jpg1; // end adress of JPG file
	fpos_t jpgc; // current adress of JPG file
	unsigned int pos = 0;
	#if defined LFS
	__int64 len = 0;
	#else
	unsigned int len;
	#endif
	
	char fn[ 256 ]; // file name of output
	char dm[ 256 ]; // dimensions of JPEG file string
	int mode = 0; // 0 -> no JPEG detected; 1 -> JPEG header; 2 -> JPEG image data
	
	int wdth = 0; // width of JPEG file
	int hght = 0; // height of JPEG file
	int bpp  = 0; // bpp of JPEG file
	
	int inum = 0; // internal # of JPG
	int onum = 0; // external # of JPG
	
	int c;
	
	
	// write status message to screen
	fprintf( stdout, "\nProcessing file \"%s\"...\n", filename );	
	
	// open input file
	fp_in = fopen( filename, "rb" );
	if ( fp_in == NULL ) {
		fprintf( stdout, "error: couldn't open input file\n" );
		return;
	}
	
	// proceed through file
	while ( ( c = getc( fp_in ) ) != EOF ) {
		switch ( mode ) {
			case 0:
				if ( c == 0xFF ) {
					c = getc( fp_in );
					if ( c == 0xD8 ) {
						fgetpos( fp_in, &jpg0 );
						jpg0 -= 2;
						mode = 1;
					}
				}
				break;
			
			case 1:
				if ( c != 0xFF ) mode = 0;
				else {
					c = getc( fp_in ); // check marker id
					if ( c == 0xD9 ) { // end of image
						pos = (unsigned int) jpg0;
						// find free filename
						while( true ) {
							sprintf( fn, "%s_%.8X_%s_%i.jpg", filename, pos, dm, onum );
							fp_out = fopen( fn, "rb" );
							if ( fp_out == NULL ) break;
							else {
								fclose( fp_out );
								onum++;
							}
						}
						while( fp_out == NULL ) fp_out = fopen( fn, "wb" );
						fgetpos( fp_in, &jpg1 );
						len = jpg1 - jpg0;						
						fprintf( stdout, "%.3i: %s JPEG @ 0x%.8X -> %14s; %6ikb\n", ++inum, ( ( jpg0 == 0 ) && ( getc( fp_in ) == EOF ) ) ? "original" : "embedded", pos, dm, len >> 10 );
						fsetpos( fp_in, &jpg0 );
						while ( len-- > 0 ) {
							c = getc( fp_in );
							putc( c, fp_out );
						}
						fclose( fp_out );
						fsetpos( fp_in, &jpg0 );
						getc( fp_in );
						getc( fp_in );
						mode = 0;
					}
					else {
						if ( c == 0xDA ) mode = 2; // start-of-scan
						len = getc( fp_in );
						len = ( (len << 8) + (getc( fp_in )&255) ) - 2;
						if ( ( c == 0xC0 ) || ( c == 0xC1 ) || ( c == 0xC2 ) ) { // sof marker
							bpp = getc( fp_in );
							hght = getc( fp_in );
							hght = ( hght << 8 ) + ( getc( fp_in ) % 255 );
							wdth = getc( fp_in );
							wdth = ( wdth << 8 ) + ( getc( fp_in ) % 255 );
							bpp = bpp * getc( fp_in );
							len = len - 6;
							sprintf( dm, "%ix%ix%i", wdth, hght, bpp );
						}
						while( len-- > 0 )
							getc( fp_in );
					}
				}
				break;
				
			case 2:
				if ( c == 0xFF ) {
					c = getc( fp_in ); // check id
					if ( ( c != 0x00 ) &&
						 ( c != 0xD0 ) &&
						 ( c != 0xD1 ) &&
						 ( c != 0xD2 ) &&
						 ( c != 0xD3 ) &&
						 ( c != 0xD4 ) &&
						 ( c != 0xD5 ) &&
						 ( c != 0xD6 ) &&
						 ( c != 0xD7 ) ) {
						fgetpos( fp_in, &jpgc );
						jpgc -= 2;
						fsetpos( fp_in, &jpgc );
						mode = 1;
					}
				}
				break;
		}
	}
	fclose( fp_in );
	fprintf( stdout, "-> done, %i JPEG file(s) found\n", inum );
}

int main(int argc, char** argv)
{
	#if defined LFS
	fprintf( stdout, "extrJPG v%s (LFS) by Matthias Stirner\n", VER );
	#else
	fprintf( stdout, "extrJPEG v%s by Matthias Stirner\n", VER );
	#endif
	
	if ( argc < 2 )
		fprintf( stdout, "Usage: %s [name of file(s)]\n", argv[ 0 ] );
	
	for ( argv++; --argc > 0; argv++ )
		dump_JPEGs( (*argv) );	
	
	#if defined PAUSE
		// pause before exit
		fprintf( stdout, "\n\n" );
		fprintf( stdout, "< press ENTER >\n" );
		fgetc( stdin );
	#endif
}
