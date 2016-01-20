#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys\stat.h>

#define ABS(v1)				( (v1 < 0) ? -v1 : v1 )
#define ABSDIFF( v1, v2 )	( (v1 > v2) ? (v1 - v2) : (v2 - v1) )
#define VPOS( p, w, h )		( (p / h) + ((p % h) * w) )
#define MEM_ERRMSG			"out of memory error"


/* -----------------------------------------------
	declarations - functions (main)
	----------------------------------------------- */

void initialize_options( int argc, char** argv );
void output_config( void );
void showhelp( void );

bool makepgm( void );
bool conv_raw( void );
bool calc_dim( void );
bool proc_img( void );
bool linesmatch( bool smart );

signed int convert_uval( unsigned int uval );
unsigned int get_uval( unsigned char* cval );
void create_img_fn( void );


/* -----------------------------------------------
	declarations - enumerations
	----------------------------------------------- */
	
enum METHOD { fixd = 1, fixa = 2, lma_s = 3, lma_f = 4, autofind = 5 };
enum VCONV  { flatten = 1, scale = 2, logarithmic = 3 };
enum PPROC  { none = 1, signum = 2, invert = 3, emphasize = 4 };
enum PTYPE  { pbm = 4, pgm = 5, ppm = 6 };


/* -----------------------------------------------
	declarations - program info
	----------------------------------------------- */
	
static char* headline = "--- makepgm v2.0a by Matthias Stirner ---";
static char* comment = "created by makepgm";


/* -----------------------------------------------
	declarations - storage
	----------------------------------------------- */

char*          raw_fn; // raw filename
char*          img_fn; // image filename
char           statusmessage[ 512 ];

int            rawsize;
int            imgsize;

unsigned char* rawdata;
unsigned char* imgdata;

char*          pchar;  // identifier char for image (= extension)
int            bpp;    // bit per pixel

unsigned int   range;  // range of values
unsigned int   middle; // middle value


/* -----------------------------------------------
	declarations - settings
	----------------------------------------------- */
	
METHOD method   = autofind;  // method for image conversion, default: auto
VCONV  vconv    = flatten;   // value conversion method, default: flatten
PPROC  pproc    = none;      // image post processing, default: none
PTYPE  ptype    = pgm;       // output image type, default: pgm

int    bpv      = 1;         // byte (!) per value
bool   vsigned  = false;     // value signed or not

int    width    = 768;       // width of picture
int    height   = 512;       // height of picture
bool   widthfixed  = false;  // if true, width is fixed
bool   heightfixed = false;  // if true, height is fixed

double aspx     = 4.0;       // aspect ratio x
double aspy     = 3.0;       // aspect ratio y
double aspect = aspx / aspy; // aspect for picture

char** filelist = NULL;	     // list of files to process 
int    hdrsize  = 0;         // headersize, subtracted from raw size



/* -----------------------------------------------
	reads in arguments, calls other functions
	----------------------------------------------- */
		
int main( int argc, char** argv )
{	
	int errorcount = 0;
	int convcount = 0;
	
	
	// print headline
	printf( "\n%s\n\n", headline );
	
	// initialize options
	initialize_options( argc, argv );
	// ...then check if something went wrong
	if ( ( (widthfixed)  && (width < 1)  ) ||
		 ( (heightfixed) && (height < 1) ) ||
		 ( aspect <= 0 ) || ( *filelist == NULL ) ||
		 ( bpv > 4 ) || ( bpv < 1 ) ) {
		showhelp();
		return 0;
	}	
	// output config
	output_config();	
	
	// main function loop
	while ( *filelist != NULL ) {	
		raw_fn = *(filelist++);
		create_img_fn();
		
		printf( "...processing \"%s\" -> ", raw_fn );
		if ( makepgm() ) {
			printf( "%s\n",	statusmessage );
			convcount++;
		}
		else {
			printf( "error: %s\n", statusmessage );
			errorcount++;
		}
	}
	
	printf( "\n\n-> %i file(s) converted, %i errors\n", convcount, errorcount );
	
	
	return 1;
}


/* -----------------------------------------------
	reads in commandline arguments
	----------------------------------------------- */
	
void initialize_options( int argc, char** argv )
{
	char** tmp_flp;
	int tmp_val;
	int i;
	
	
	// get memory for filelist & preset with NULL
	filelist = new char*[ argc ];
	for ( i = 0; i < argc; i++ )
		filelist[ i ] = NULL;
	
	// preset temporary filelist pointer
	tmp_flp = filelist;
	
	
	// read in arguments
	while ( --argc > 0 ) {
		argv++;		
		// switches begin with '-'
		if ( sscanf( (*argv), "-hd%i", &tmp_val ) == 1 ) {
			hdrsize = tmp_val;
		}
		else if ( sscanf( (*argv), "-a%lf:%lf", &aspx, &aspy ) == 2 ) {
	   		aspect = ( (aspx > 0) && (aspy > 0) ) ? ( aspx / aspy ) : 0;
			method = fixa;
	  		widthfixed = false;
	  		heightfixed = false;	  			
   		}
		else if ( sscanf( (*argv), "-w%i", &tmp_val ) == 1 ) {
			width = tmp_val;
			method = fixd;
			widthfixed = true;
		}
		else if ( sscanf( (*argv), "-h%i", &tmp_val ) == 1 ) {
			height = tmp_val;
			method = fixd;
			heightfixed = true;
		}
		else if ( sscanf( (*argv), "-s%i", &tmp_val ) == 1 ) {
			bpv = tmp_val;
			vsigned = true;
		}
		else if ( sscanf( (*argv), "-u%i", &tmp_val ) == 1 ) {
			bpv = tmp_val;
			vsigned = false;
		}
		else if ( strcmp( (*argv), "-lms" ) == 0 ) {
			method = lma_s;
	   		widthfixed = false;
	   		heightfixed = false;
	   	}
		else if ( strcmp( (*argv), "-lmf" ) == 0 ) {
			method = lma_f;
	   		widthfixed = false;
	   		heightfixed = false;
	   	}
   		else if ( strcmp( (*argv), "-flt" ) == 0) {
			vconv = flatten;
   		}
   		else if ( strcmp( (*argv), "-scl" ) == 0) {
			vconv = scale;
   		}
   		else if ( strcmp( (*argv), "-log" ) == 0) {
			vconv = logarithmic;
   		}
		else if ( strcmp( (*argv), "-sgn" ) == 0) {
			pproc = signum;
   		}
		else if ( strcmp( (*argv), "-neg" ) == 0) {
			pproc = invert;
   		}
		else if ( strcmp( (*argv), "-emz" ) == 0) {
			pproc = emphasize;
   		}
   		else if ( strcmp( (*argv), "-ppm" ) == 0) {
			ptype = ppm;
   		}
   		else if ( strcmp( (*argv), "-pgm" ) == 0) {
			ptype = pgm;
   		}
   		else if ( strcmp( (*argv), "-pbm" ) == 0) {
			ptype = pbm;
   		}
		else {
			// if argument is not switch, it's a filename
			*(tmp_flp++) = *argv;
		}		
	}
	
	// set bpp & pchar
	switch ( ptype ) {
	
		case pbm:
			pchar = "pbm";
			bpp = 1;
			break;
			
		case pgm:
			pchar = "pgm";
			bpp = 8;
			break;
			
		case ppm:
			pchar = "ppm";
			bpp = 24;
			break;
	}
	
	// set range & middle
	middle = 1 << ( (8 * bpv) - 1 );
	range  = ( middle - 1 ) + middle;
}


/* -----------------------------------------------
	output configuration to screen
	----------------------------------------------- */
	
void output_config( void )
{	
	// method info
	switch ( method ) {	
		case fixd:
			if ( widthfixed )  printf( "-> width fixed at %i\n", width );
			if ( heightfixed ) printf( "-> height fixed at %i\n", height );
			break;			
		case fixa:
			printf( "-> aspect ratio fixed at %i:%i\n", (int) aspx, (int) aspy );
			break;
		case lma_s:
			printf( "-> width/height is matched (smart)\n" );
			break;
		case lma_f:
			printf( "-> width/height is matched (full)\n" );
			break;
		case autofind:
			printf( "-> width/height is found automatic\n" );
			break;
	}
	
	// value conversion info
	printf( "-> value conversion method: " );
	switch ( vconv ) {	
		case flatten:
			if ( bpv > 1 ) printf( "flat" );
			else printf( "none" );
			break;			
		case scale:
			if ( bpv > 1 ) printf( "scale" );
			else printf( "none" );
			break;
		case logarithmic:
			printf( "log" );
			break;	
	}
	printf( "\n" );
	
	// post processing info
	if ( pproc != none ) {
		printf( "-> image post processing: " );
		switch ( pproc ) {	
			case none:
				printf( "none" );
				break;			
			case signum:
				printf( "signum" );
				break;
			case invert:
				printf( "invert" );
				break;
			case emphasize:
				printf( "emphasize" );
				break;	
		}
		printf( "\n" );
	}
	
	// input / output info
	printf( "-> input set to %s %ibit\n", (vsigned) ? "signed" : "unsigned", bpv * 8 );
	printf( "-> output set to %ibpp %s\n", bpp, pchar );

	printf( "\n\n" );
}


/* -----------------------------------------------
	shows help in case of wrong input
	----------------------------------------------- */
	
void showhelp( void )
{
	printf( "Usage: makepgm [switches] [filename]\n" );
	printf( "\n" );
	printf( "   [-hd??] -> subtract ?? byte from begin of data\n" );
	printf( "   [-w???] -> fixed width for image, must be > 0\n" );
	printf( "   [-h???] -> fixed height for image, must be > 0\n" );
	printf( "   [-a?:?] -> fixed aspect ratio, must be > 0\n" );
	printf( "    [-lms] -> do smart linematch to find width\n" );
	printf( "    [-lmf] -> do full linematch to find width\n" );
	printf( "     [-s?] -> set input value signed / ? byte, max 4\n" );	
	printf( "     [-u?] -> set input value unsigned / ? byte, max 4\n" );
	printf( "\n" );
	printf( "    [-flt] -> flatten values > 1 byte\n" );
	printf( "    [-scl] -> scales values to 1 byte (linear)\n" );
	printf( "    [-log] -> scales values to 1 byte (logarithmic)\n" );
	printf( "\n" );
	printf( "    [-sgn] -> convert values to signum\n" );
	printf( "    [-neg] -> invert image colors\n" );
	printf( "    [-emz] -> emphasize zeroes\n" );
	printf( "\n" );
	printf( "    [-ppm] -> create 24bpp colored ppm-image\n" );
	printf( "    [-pgm] -> create 8bpp grayscale pgm-image\n" );
	printf( "    [-pbm] -> create 1bpp monochrome pbm-image\n" );	
	printf( "\n" );
	printf( "Examples: \"makepgm -hd36 lena.pgm\"\n"  );
	printf( "          \"makepgm -w512 -h512 -ppm lena.raw\"\n" );
	printf( "          \"makepgm -a16:9 -ppm kodim01.raw\"\n" );
	printf( "          \"makepgm -lms -u2 data.hex\"\n" );
}


/* -----------------------------------------------
	converts one file
	----------------------------------------------- */
	
bool makepgm( void )
{
	FILE* fp;
	
	
	// open data file
	fp = fopen( raw_fn, "rb" );
	if ( fp == NULL ) {
		sprintf( statusmessage, "couldn't read data" );
		return false;
	}
	
	// find out sizes, subtract headerinfo
	fseek( fp, 0, SEEK_END );	
	rawsize = ftell( fp ) - hdrsize;
	imgsize = rawsize / bpv;
	if ( rawsize <= 0 ) {
		if ( hdrsize > 0 ) sprintf( statusmessage, "too small for header subtraction" );
		else sprintf( statusmessage, "file is empty" );
		return false;
	}
	
	// alloc memory for raw data
	rawdata = ( unsigned char* ) calloc ( rawsize, 1 );
	if ( rawdata == NULL ) {
		sprintf( statusmessage, MEM_ERRMSG );
		return false;
	}
	// alloc memory for image data
	imgdata = ( unsigned char* ) calloc ( imgsize, 1 );
	if ( imgdata == NULL ) {
		sprintf( statusmessage, MEM_ERRMSG );
		return false;
	}
	
	// read in raw data and close file
	fseek( fp, hdrsize, SEEK_SET );
	fread( rawdata, 1, rawsize, fp );
	fclose( fp );	
	
	// convert raw data to image data 
	if ( !conv_raw() )
		return false;
	
	// find out correct width & height
	if ( !calc_dim() )
		return false;
	
	// post process image data
	if ( !proc_img() )
		return false;
	
	// open image file
	fp = fopen( img_fn, "wb" );
	if ( fp == NULL ) {
		sprintf( statusmessage, "couldn't write %s", pchar );
		return false;
	}
	
	// write PBM/PGM/PPM header
	fprintf( fp, "P%i\n", ptype );
	fprintf( fp, "# %s\n", comment );
	fprintf( fp, "%i %i\n", width, height );
	if ( ptype != pbm )
		fprintf( fp, "255\n" );
		
	// write image data & close file
	fwrite( imgdata, 1, imgsize, fp );
	fclose( fp );
	
	// free allocated memory
	free( rawdata );
	free( imgdata );	
	
	// success message
	sprintf( statusmessage, "converted to %ix%ix%i %s",
				width, height, bpp, pchar );
				
	
	return true;
}


/* -----------------------------------------------
	convert raw data to image data 
	----------------------------------------------- */
	
bool conv_raw( void )
{
	signed int sval;
	unsigned int uval;
	int scs; // linear scale factor
	double scl; // logarithmic scale factor
	int i;
	
	
	if ( vsigned ) { // if value is signed
	
		scs = ( 1 << ((bpv - 1) * 8) );
		scl = 128.0 / ( bpv * log2( 128.0 ) );
		
		for ( i = 0; i < imgsize; i++ )
		{
			uval = get_uval( &(rawdata[ i * bpv ]) );
			sval = ( uval >= middle ) ? ( uval - range - 1 ) : uval;
			
			switch ( vconv ) {
				case flatten: // "flatten" - method
					sval = ( sval < -128 ) ? -128 : sval;
					sval = ( sval >  127 ) ?  127 : sval;
					break;
				case scale: // "scale" - method
					sval /= scs;
					break;
				case logarithmic: // "logarithmic" - method
					if ( sval < 0 )
						sval = (signed int) -floor( log2( ABS(sval) ) * scl );
					if ( sval > 0 )
						sval = (signed int) floor( log2( ABS(sval) ) * scl );
					break;
			}
	   	
			imgdata[ i ] = (unsigned char) ( sval + 128 );
		}		
	}
	else { // if value is not signed
	
		scs = ( 1 << ((bpv - 1) * 8) );
		scl = 256.0 / ( bpv * log2( 256.0 ) );
		
		for ( i = 0; i < imgsize; i++ )
		{
			uval = get_uval( &(rawdata[ i * bpv ]) );			

			switch ( vconv ) {
				case flatten: // "flatten" - method
					uval = ( uval > 255 ) ? 255 : uval;
					break;
				case scale: // "scale" - method
					uval /= scs;
					break;
				case logarithmic: // "logarithmic" - method
					if ( uval > 0 )
						uval = (unsigned int) floor( log2( uval ) * scl );
					break;
			}
		
			imgdata[ i ] = (unsigned char) uval;
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	calculate correct dimensions for image
	----------------------------------------------- */
	
bool calc_dim( void )
{
	// find out correct width & height
	switch ( method )
	{	
		case fixd: // fixed dimension method
			if ( ( widthfixed ) && ( !heightfixed ) )
				height = (int) ceil( imgsize / ( 0.125 * bpp * width ) );
			if ( ( !widthfixed ) && ( heightfixed ) )
				width  = (int) ceil( imgsize / ( 0.125 * bpp * height) );
			break;
			
		case fixa: // fixed aspect method
			width  = (int) ceil( sqrt( (aspect * imgsize) / (0.125 * bpp) ) );
			height = (int) ceil( imgsize / ( 0.125 * bpp * width ) );
			break;
			
		case lma_s: // smart linesmatching
			if ( !linesmatch( true ) ) return false;
			height = (int) ceil( imgsize / ( 0.125 * bpp * width ) );
			break;
			
		case lma_f: // full linesmatching
			if ( !linesmatch( false ) ) return false;
			height = (int) ceil( imgsize / ( 0.125 * bpp * width ) );
			break;
			
		case autofind:
			if ( !linesmatch( true ) )
			if ( !linesmatch( false ) )
			width  = (int) ceil( sqrt( (aspect * imgsize) / (0.125 * bpp) ) );
			height = (int) ceil( imgsize / ( 0.125 * bpp * width ) );
			break;
	}
	
	
	return true;
}


/* -----------------------------------------------
	post-processes image data
	----------------------------------------------- */
	
bool proc_img( void )
{
	int i;
	
	switch ( pproc )
	{
		case none: // no postprocessing
			break;
		case signum: // signum postprocessing
			for ( i = 0; i < imgsize; i++ ) {
				imgdata[ i ] = ( imgdata[ i ] > 128 ) ? 255 : imgdata[ i ];
				imgdata[ i ] = ( imgdata[ i ] < 128 ) ?   0 : imgdata[ i ];
			}
			break;
		case invert: // invert image colors
			for ( i = 0; i < imgsize; i++ )
				imgdata[ i ] = 255 - imgdata[ i ];
			break;
		case emphasize: // emphasize zeroes
			if ( vsigned )
				for ( i = 0; i < imgsize; i++ )
					imgdata[ i ] = ( imgdata[ i ] == 128 ) ? 0 : 255;
			else
				for ( i = 0; i < imgsize; i++ )
					imgdata[ i ] = ( imgdata[ i ] == 0 ) ? 0 : 255;
	}
	
	
	return true;
}


/* -----------------------------------------------
	matches lines of image to find right width
	----------------------------------------------- */
	
bool linesmatch( bool smart )
{
	int     wmin     = 8 * bpp;            // minimum width => 64 pixels / 8 bit
	int     hmin     = 8 * bpp;            // minimum height => 64 pixels / 8 bit
	int     wmax     = ( imgsize / hmin ); // maximum width (calculated)
	int     maxsize  = 64 * 1024 * 1024;   // maximum size for image
	int     minsize  = wmin * hmin;        // minimum size for image
	
	int*    llen;    // lengths of lines to test
	double* lmvl;    // line match values - the lower the better
	int     plenc;   // count of possible line lengths
	int     alenc;   // count of active line lenghts
	
	double  avrg;
	int     oldpos, newpos;
	int     v1, v2;
	int     w, h;
	
	int     bmatch;
	bool    stuck;
	int     i, j;
	
	
	// security checks
	if ( imgsize < minsize ) {
		sprintf( statusmessage, "file is too small for matching" );
		return false;
	}	
	if ( imgsize > maxsize ) {
		sprintf( statusmessage, "file is too big for matching" );
		return false;
	}
	
	
	if ( smart ) // for smart lines matching
	{
		// count (smart) possible lengths of lines
		plenc = 0;
		for ( w = wmin; w <= wmax; w++ )
			if ( imgsize % w == 0 )	plenc++;
		
		alenc = plenc;
		
		// security check
		if ( plenc == 0 ) {
			sprintf( statusmessage, "possible width not found" );
			return false;
		}
		
		// alloc memory
		llen = ( int* ) calloc ( plenc, sizeof( int ) );
		if ( llen == NULL ) {
			sprintf( statusmessage, MEM_ERRMSG );
			return false;
		}
		lmvl = ( double* ) calloc ( plenc, sizeof( double ) );
		if ( lmvl == NULL ) {
			sprintf( statusmessage, MEM_ERRMSG );
			return false;
		}
		
		// fill linelength array
		i = 0;		
		for ( w = wmin; w <= wmax; w++ ) {
			if ( imgsize % w == 0 )	{
				llen[ i ] = w;
				lmvl[ i ] = 0;
				i++;
			}
		}
	}
	else // for full lines matching
	{
		plenc = wmax - wmin + 1;
		alenc = plenc;
		
		// alloc memory
		llen = ( int* ) calloc ( plenc, sizeof( int ) );
		if ( llen == NULL ) {
			sprintf( statusmessage, MEM_ERRMSG );
			return false;
		}
		lmvl = ( double* ) calloc ( plenc, sizeof( double ) );
		if ( lmvl == NULL ) {
			sprintf( statusmessage, MEM_ERRMSG );
			return false;
		}		
		
		//  use every possible & impossible line length
		for ( i = 0; i < plenc; i++ ) {
			llen[ i ] = i + wmin;
			lmvl[ i ] = 0;
		}
	}
	
	// exclude linelengths not matching for format
	if  ( ptype == ppm )
	for ( i = 0; i < plenc; i++ )
	if  ( llen[ i ] % 3 > 0 ) {
		llen[ i ] = 0;
		alenc--;
	}
	
	// the main part: matching lines
	newpos = 0;
	while ( alenc > 1 )
	{
		avrg = 0;
		oldpos = newpos;
		newpos = (int) ceil( (double) ( imgsize - 1 ) / ( alenc - 1 ) );
		newpos = ( newpos < hmin ) ? hmin : newpos;
				
		// calculate line match values
		for ( i = 0; i < plenc; i++ )
		if ( llen[ i ] > 0 )
		{			
			w = llen[ i ];
			h = imgsize / llen[ i ];
			
			// lines matching
			v1 = imgdata[ VPOS( oldpos, w, h ) ];
			for ( j = oldpos + 1; j <= newpos; j++ ) {
				v2 = v1;
				v1 = imgdata[ VPOS( j, w, h ) ];				
				lmvl[ i ] += ABSDIFF( v1, v2 );
			}
			
			avrg += lmvl[ i ];				
		}
		
		// compare match values, throw out if matchvalue > average
		avrg /= alenc;
		stuck = true;
		for ( i = 0; i < plenc; i++ )
		if ( ( llen[ i ] > 0 ) && ( lmvl[ i ] > avrg ) ) {
			stuck = false;
			llen[ i ] = 0;
			alenc--;
		}
		
		// to prevent endless loops with one-color images
		if ( stuck ) {
			sprintf( statusmessage, "could not find match" );
			return false;
		}
	}
	
	// get best match = last one standing :-)
	bmatch = plenc - 1;
	while ( llen[ bmatch ] == 0 )
		bmatch--;
	
	// ...and write it to width
	width = ( llen[ bmatch ] * 8 ) / bpp;
	
	// free memory
	free( llen );
	free( lmvl );
	
	
	return true;	
}


/* -----------------------------------------------
	gets uint from array of chars
	----------------------------------------------- */
	
unsigned int get_uval( unsigned char* cval )
{
	unsigned int uval = 0;
	int i;
	
	for ( i = 0; i < bpv; i++ )
		uval += ( ( (unsigned char) cval[i] ) << ( 8 * i ) );
	
	return uval;
}


/* -----------------------------------------------
	creates filename for image
	----------------------------------------------- */
	
void create_img_fn( void )
{
	int i;
	
	
	img_fn = new char[ strlen( raw_fn ) + 8 ];
	
	// copy filename + '.' + extension into one string -> imgfilename
	sprintf( img_fn, "%s.%s", raw_fn, pchar );
	
	// find position of dot in filename and replace it with '_'
	for ( i = strlen( raw_fn ) - 1; i >= 0; i-- )
	if ( raw_fn[ i ] == '.' ) {
		img_fn[ i ] = '_';
		break;
	}
}
