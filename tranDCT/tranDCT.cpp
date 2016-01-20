#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <ctime>

#if defined(UNIX) || defined (__LINUX__)
	#include <unistd.h>
#else
	#include <io.h>
#endif

#define ABS(v1)			( (v1 < 0) ? -v1 : v1 )
#define ABSDIFF(v1,v2)	( (v1 > v2) ? (v1 - v2) : (v2 - v1) )
#define ROUND_F(v1)		( (v1 < 0) ? (int) (v1 - 0.5) : (int) (v1 + 0.5) )
#define CLAMPED(l,h,v)	( ( v < l ) ? l : ( v > h ) ? h : v )

#define DCT_BCH(nx)		( ( imgwidth   + ( nx - 1 ) ) / nx )
#define DCT_BCV(ny)		( ( imgheight  + ( ny - 1 ) ) / ny )
#define DCT_BCS(nx,ny)	( DCT_BCH(nx) * DCT_BCV(ny) )

#define DCT_IMGW(nx)	( DCT_BCH( nx ) * nx )
#define DCT_IMGH(ny)	( DCT_BCV( ny ) * ny )
#define DCT_IMGS(nx,ny)	( DCT_IMGW( nx ) * DCT_IMGH( ny ) )

#define COS_DCT(l,s,n)  ( cos( ( ( 2 * l + 1 ) * s * M_PI ) / ( 2 * n ) ) )
#define C_DCT(n)		( ( n == 0 ) ? ( 1 ) : ( sqrt( 2 ) ) )
#define CFPOS(x,y,u,v)	( ( ( ( y * bx_c ) + x ) * ns_c ) + ( ( v * nx_c ) + u ) )

#define MEM_ERRMSG	"out of memory error"
#define FRD_ERRMSG	"could not read file / file not found: %s"
#define FWR_ERRMSG	"could not write file / file writeprotected: %s"
#define	WR_C(f,c,p)	( out_format == CIMG ) ? fwrite( &( p = CLAMPED( 0, 255, c + 128 ) ), sizeof( char ), 1, f ) : fwrite( &c, sizeof( CTYPE ), 1, fp )

#define CSV_SEP		';'
#define MISS_REP	'-'

#define FAST_DCT	true
#define CTYPE		int
#define DCT_SCALE	sqrt( 8 )
#define NX_MAX		64
#define NY_MAX		64
#define NX_DEF		8
#define NY_DEF		8


/* -----------------------------------------------
	struct & enum declarations
	----------------------------------------------- */

enum O_FORM
			{ 	IMG   =  1, COEF  =  2, CSV   =  3,
				RAWP  =  4,	CRAW  =  5,	CIMG  =  6 };
				
enum I_FORM
			{   PGM = 0, DCT = 1, TDS = 2, UNK = 3 };
			
struct TDS_COL
{
	int   t_nx; // blocksize x
	int   t_ny; // blocksize y
	
	int   x_stp; // step x
	int   y_stp; // step y
	
	int   x_pos; // relative position x
	int   y_pos; // relative position y
	int   u_pos; // absolute position u
	int   v_pos; // absolute position v
	
	int   x_til; // relative position x (last)
	int   y_til; // relative position y (last)
	int   u_til; // absolute position u (last)
	int   v_til; // absolute position v (last)
	
	char  name[ 256 ]; // name of column
	
	int   numrows; // number of elements in column
	int   nummiss; // number of missing elements
	signed CTYPE*  coldata; // column data
	unsigned char* chklist; // check list (incomplete yes/no)
	
	TDS_COL* prev; // previous column
	TDS_COL* next; // next column
};


/* -----------------------------------------------
	function declarations: main interface
	----------------------------------------------- */

void initialize_options( int argc, char** argv );
void process_file( void );
void execute( bool (*function)() );
void get_status( bool (*function)() );
void show_help( void );


/* -----------------------------------------------
	function declarations: main functions
	----------------------------------------------- */
	
bool check_file( void );
bool parse_tdasfile( FILE* fp );
bool trans_image( void );
bool trans_coeffs( void );
bool write_file( void );
bool read_column( void );
bool write_csv( void );
bool reset_buffers( void );


/* -----------------------------------------------
	function declarations: DCT
	----------------------------------------------- */

bool  prepare_dct( int nx, int ny );
float idct_2d_bnx( signed CTYPE* F, int ix, int iy );
float idct_2d_bny( signed CTYPE* F, int ix, int iy );
float fdct_2d_bnx( unsigned char* f, int iu, int iv );
float fdct_2d_bny( unsigned char* f, int iu, int iv );
float idct_2d_fst( signed CTYPE* F, int ix, int iy );
float fdct_2d_fst( unsigned char* f, int iu, int iv );
void  block_idct_2d( signed CTYPE* F, unsigned char* f );
void  block_fdct_2d( unsigned char* f, signed CTYPE* F );


/* -----------------------------------------------
	function declarations: miscelaneous helpers
	----------------------------------------------- */

char* fetch_command( char* text, char* cmd );
char* create_filename( const char* base, const char* extension );
void set_extension( char* destination, const char* origin, const char* extension );
void add_underscore( char* filename );


/* -----------------------------------------------
	global variables: data storage
	----------------------------------------------- */

signed CTYPE*   coefdata = NULL;	// DCT coefficients
unsigned char*  imgdata  = NULL;	// image data

int    imgsize   = 0;		// size of image
int    imgwidth  = 0;		// width of image
int    imgheight = 0;		// height of image

int    cofsize   = 0;		// size of coefficient image;
int    cofwidth  = 0;		// width of coefficient image
int    cofheight = 0;		// heigt of coefficient image


/* -----------------------------------------------
	global variables: DCT data
	----------------------------------------------- */

float (*idct_2d)( signed CTYPE*, int, int )  = NULL; // inverse dct  function
float (*fdct_2d)( unsigned char*, int, int ) = NULL; // forward dct function

float* icos_idct_nx = NULL;		// precalculated values for idct, nx table
float* icos_idct_ny = NULL;		// precalculated values for idct, ny table
float* icos_fdct_nx = NULL;		// precalculated values for fdct, nx table
float* icos_fdct_ny = NULL;		// precalculated values for fdct, ny table

float* icos_idct_fst = NULL;	// precalculated values for idct, fast table
float* icos_fdct_fst = NULL;	// precalculated values for fdct, fast table

int    bs_c      = 0;		// block count (current)
int    bx_c      = 0;		// block count x (current)
int    by_c      = 0;		// block count y (current)

int    ns_c      = 0;		// dct block size (current)
int    nx_c      = 0;		// dct block size x (current)
int    ny_c      = 0;		// dct block size y (current)


/* -----------------------------------------------
	global variables: dependency examination
	----------------------------------------------- */

TDS_COL* curr_col = NULL; // current column
bool inc_incp = false; // include incomplete lines yes/no
bool cmt_incp = false; // comment before incomplete lines yes/no
bool chk_ssiz = false; // check wether all columns have the same size
bool inc_name = false; // include names in header of csv file
bool inc_chkl = false; // include check list in csv file
bool inc_desc = false; // include descriptions (first col)
bool wrt_gplt = false; // write gnuplot .plt file


/* -----------------------------------------------
	global variables: info about files
	----------------------------------------------- */
	
char*  inpfilename = NULL;	// name of input file
char*  outfilename = NULL;	// name of output file
char*  tmpfilename = NULL;	// temporary file name

I_FORM in_format  = UNK;	// type of input file
O_FORM out_format = IMG;	// type of output file

char** filelist  = NULL;	// list of files to process 
int    file_cnt  = 0;		// count of files in list
int    file_no   = 0;		// number of current file

char   pgm_comment[ 256 ];	// pgm file comment (gets trashed, anyways)


/* -----------------------------------------------
	global variables: messages
	----------------------------------------------- */

char statusmessage[ 128 ];
char errormessage [ 128 ];
bool (*errorfunction)();
int  errorlevel;
// meaning of errorlevel:
// -1 -> wrong input
// 0 -> no error
// 1 -> warning
// 2 -> fatal error


/* -----------------------------------------------
	global variables: settings
	----------------------------------------------- */

bool   verbosity  = false;	// level of verbosity
bool   overwrite  = false;	// overwrite files yes / no
int    err_tresh  = 1;		// error threshold ( proceed on warnings yes (2) / no (1) )

int    nx_s		 = NX_DEF;	// dct block size x setting
int    ny_s      = NY_DEF;	// dct block size y setting

bool   def_actn  = true;	// use default action yes/no
int    outp_mode = 0;		// output mode for DCT coefficients
FILE*  msgout    = stdout;	// stream for output of messages


/* -----------------------------------------------
	global variables: info about program
	----------------------------------------------- */

const unsigned char prgversion   = 9;
static const char*  subversion   = "b";
static const char*  versiondate  = "05/11/2008";
static const char*  author       = "Matthias Stirner";


/* -----------------------------------------------
	main-function
	----------------------------------------------- */

int main( int argc, char** argv )
{	
	sprintf( statusmessage, "no statusmessage specified" );
	sprintf( errormessage, "no errormessage specified" );
	
	clock_t begin, end;
	
	int error_cnt = 0;
	int warn_cnt  = 0;
	int speed     = 0;
	
	
	// read options from command line
	initialize_options( argc, argv );
	
	// write program info to screen
	fprintf( msgout,  "\n--> tranDCT v%i.%i%s (%s) by %s <--\n\n",
			prgversion / 10, prgversion % 10, subversion, versiondate, author );
	
	// check if user input is wrong, show help screen if it is
	if ( file_cnt == 0 ) {
		show_help();
		return -1;
	}
	
	// (re)set program has to be done first
	reset_buffers();
	
	// process file(s) - this is the main function routine
	begin = clock();
	for ( file_no = 0; file_no < file_cnt; file_no++ ) {		
		process_file();
		if ( errorlevel >= err_tresh ) error_cnt++;
		else {
			if ( errorlevel == 1 ) warn_cnt++;
		}
	}
	end = clock();
	
	// calculate speed
	speed = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC ); 
	
	// show statistics
	fprintf( msgout,  "\n\n-> %i file(s) processed in %i ms, %i error(s), %i warning(s)\n", file_cnt, speed, error_cnt, warn_cnt );
	
	
	return file_cnt;
}


/* ----------------------- Begin of main interface functions -------------------------- */


/* -----------------------------------------------
	reads in commandline arguments
	----------------------------------------------- */
	
void initialize_options( int argc, char** argv )
{
	char** tmp_flp;
	int i;
	
	
	// get memory for filelist & preset with NULL
	filelist = ( char** ) calloc( argc, sizeof( char* ) );
	// filelist = new char*[ argc ];
	for ( i = 0; i < argc; i++ )
		filelist[ i ] = NULL;
	
	// preset temporary filelist pointer
	tmp_flp = filelist;
	
	
	// read in arguments
	while ( --argc > 0 ) {
		argv++;		
		// switches begin with '-'
		if ( strcmp((*argv), "-v" ) == 0 ) {
			verbosity = true;
		}
		else if ( strcmp((*argv), "-o" ) == 0 ) {
			overwrite = true;
		}
		else if ( strcmp((*argv), "-p" ) == 0 ) {
			err_tresh = 2;
		}
		else if ( sscanf( (*argv), "-b%ix%i", &nx_s, &ny_s ) == 2 ) {			
			nx_s = ( nx_s < 1 ) ? 1 : nx_s;
			ny_s = ( ny_s < 1 ) ? 1 : ny_s;
			nx_s = ( nx_s > NX_MAX ) ? NX_MAX : nx_s;
			ny_s = ( ny_s > NY_MAX ) ? NY_MAX : ny_s;
		}		
		else if ( strcmp( (*argv), "-img" ) == 0 ) {
			out_format = IMG;
			def_actn   = false;
		}
		else if ( strcmp( (*argv), "-dct" ) == 0 ) {
			out_format = COEF;
			def_actn   = false;
		}
		else if ( strcmp( (*argv), "-rawp" ) == 0 ) {
			out_format = RAWP;
			def_actn   = false;
		}
		else if ( sscanf( (*argv), "-craw%i", &outp_mode ) == 1 ) {
			outp_mode  = CLAMPED( 0, 3, outp_mode );
			out_format = CRAW;
			def_actn   = false;
		}
		else if ( sscanf( (*argv), "-cimg%i", &outp_mode ) == 1 ) {
			outp_mode  = CLAMPED( 0, 3, outp_mode );
			out_format = CIMG;
			def_actn   = false;
		}
		else if ( strcmp( (*argv), "-csv" ) == 0 ) {
			out_format = CSV;
			def_actn   = false;
		}
		else {
			// if argument is not switch, it's a filename
			*(tmp_flp++) = *argv;
		}		
	}
	
	// count number of files (or filenames) in filelist
	for ( file_cnt = 0; filelist[ file_cnt ] != NULL; file_cnt++ );
	
	// this is some kind of a 'hack', putting the first .tds file found in the file list right to the beginning
	// this is needed for comfortable drag & drop handling
	for ( i = 0; i < file_cnt; i++ ) {
		if ( strstr( filelist[ i ], ".tds" ) != NULL ) {
			filelist[ argc - 1 ] = filelist[ i ];
			for ( ; i > 0; i-- ) {
				filelist[ i ] = filelist[ i - 1 ];
				filelist[ i - 1 ] = filelist[ argc - 1 ];
			}
			break;
		}
	}
}


/* -----------------------------------------------
	processes one file
	----------------------------------------------- */

void process_file( void )
{
	clock_t begin, end;
	
	char* actionmsg  = NULL;
	char* errtypemsg = NULL;
	
	int speed;
	
	
	errorfunction = NULL;
	errorlevel = 0;
	
	
	// main function routine
	begin = clock();
	
	fprintf( msgout,  "\nProcessing file %i of %i \"%s\" -> ",
				file_no + 1, file_cnt, filelist[ file_no ] );
	if ( verbosity )
		fprintf( msgout,  "\n----------------------------------------" );
	
	// check input file and determine filetype
	actionmsg = (char*) "Parsing";
	execute( check_file );
	
	// get specific action message
	if ( in_format == UNK ) actionmsg = (char*) "unknown filetype";
	else if ( out_format == CSV ) actionmsg = (char*) "Analysing";
	else if ( in_format != TDS ) actionmsg = (char*) "Converting";
	
	if ( !verbosity ) fprintf( msgout, "%s -> ", actionmsg );
	
	if ( in_format != TDS ) switch ( out_format )
	{		
		case CSV:
			// skip if scripting file missing
			if ( curr_col == NULL )
				break;
			// go back to image if dct as input
			if ( in_format == DCT )	{
				execute( trans_coeffs );
			}
			// go back to first column	
			while ( curr_col->prev != NULL )
				curr_col = curr_col->prev;
			// go through all columns, collect data
			while ( true ) {
				if ( ( nx_c != curr_col->t_nx ) || ( ny_c != curr_col->t_ny ) ) {
					// set nx/ny
					nx_s = curr_col->t_nx;
					ny_s = curr_col->t_ny;						
					// transform image
					execute( trans_image );
				}
				//analyse dependencies
				execute( read_column );
				
				// switch to next column
				if ( curr_col->next != NULL )
					curr_col = curr_col->next;
				else break;
			}
			// write results to file
			execute( write_csv );
			break;
					
		case IMG:
		case RAWP:
			// transform image first
			if ( in_format == PGM )
				execute( trans_image );
			// transform coeffs / write file
			execute( trans_coeffs );
			execute( write_file );
			break;
			
		case COEF:
		case CRAW:
		case CIMG:
			// transform image first
			if ( in_format == PGM )
				execute( trans_image );
			// only transform if needed
			if ( ( nx_c != nx_s ) || ( ny_c != ny_s ) ) {
				execute( trans_coeffs );
				execute( trans_image );
			}
			// write output
			execute( write_file );
			break;
	}
	// reset buffers
	reset_buffers();
	
	end = clock();
	
	
	// speed calculation
	speed = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );

	
	if ( !verbosity ) {
        if ( errorlevel < err_tresh ) 
            fprintf( msgout,  "%i ms", speed );
        else fprintf( msgout,  "ERROR" );
        if ( errorlevel > 0 )
            fprintf( msgout,  "\n" );
    } else {
        fprintf( msgout,  "\n----------------------------------------\n" );
        if ( errorlevel < err_tresh ) fprintf( msgout,  "-> %s -> %i ms\n", actionmsg, speed );
	}
	
	switch ( errorlevel )
	{
		case 0:
			errtypemsg = (char*) "none";
			break;
			
		case 1:
			if ( errorlevel < err_tresh )
				errtypemsg = (char*) "warning (ignored)";
			else
				errtypemsg = (char*) "warning (skipped file)";
			break;
		
		case 2:
			errtypemsg = (char*) "fatal error";
			break;
	}
	
	if ( errorlevel > 0 )
	{
		get_status( errorfunction );
		fprintf( stderr, " %s -> %s:\n", statusmessage, errtypemsg  );
		fprintf( stderr, " %s\n", errormessage );
		if ( verbosity )
			fprintf( stderr, " (in file \"%s\")\n", filelist[ file_no ] );
	}
	
	if ( verbosity )
		fprintf( msgout,  "\n" );
}


/* -----------------------------------------------
	main-function execution routine
	----------------------------------------------- */

void execute( bool (*function)() )
{	
	clock_t begin, end;
	bool success;
	int i;
	
	
	if ( errorlevel < err_tresh )
	{
		// get statusmessage
		get_status( function );
		// write statusmessage
		if ( verbosity ) {
			fprintf( msgout,  "\n%s ", statusmessage );
			for ( i = strlen( statusmessage ); i <= 30; i++ )
				fprintf( msgout,  " " );			
		}
		
		// set starttime
		begin = clock();
		// call function
		success = ( *function )();
		// set endtime
		end = clock();
		
		if ( ( errorlevel > 0 ) && ( errorfunction == NULL ) )
			errorfunction = function;
		
		// write statusmessage
		if ( success ) {
			if ( verbosity ) fprintf( msgout,  "%6ims",
				(int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC ) );
		}
		else {
			errorfunction = function;
			if ( verbosity ) fprintf( msgout,  "%8s", "ERROR" );
		}
	}
}


/* -----------------------------------------------
	gets statusmessage for function
	----------------------------------------------- */
	
void get_status( bool (*function)() )
{	
	if ( function == NULL ) {
		sprintf( statusmessage, "unknown action" );
	}
	else if ( function == *check_file ) {
		sprintf( statusmessage, "Parsing input file" );
	}
	else if ( function == *reset_buffers ) {
		sprintf( statusmessage, "Resetting program" );
	}
	else if ( function == *trans_image ) {
		sprintf( statusmessage, "Performing forward DCT (%ix%i)", nx_s, ny_s );
	}
	else if ( function == *trans_coeffs ) {
		sprintf( statusmessage, "Performing inverse DCT (%ix%i)", nx_c, ny_c );
	}
	else if ( function == *write_file ) {
		sprintf( statusmessage, "Writing output file" );
	}
	else if ( function == *read_column ) {
		sprintf( statusmessage, "Reading one column" );
	}
	else if ( function == *write_csv ) {
		sprintf( statusmessage, "Writing CSV/PLT file" );
	}
}


/* -----------------------------------------------
	shows help in case of wrong input
	----------------------------------------------- */
	
void show_help( void )
{	
	fprintf( msgout, "\n" );
	fprintf( msgout, "Usage: tranDCT [switches] [filename(s)]");
	fprintf( msgout, "\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-v]     set verbose messages\n" );
	fprintf( msgout, " [-o]     overwrite existing files\n" );
	fprintf( msgout, " [-p]     proceed on warnings\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-b?x?]  specify blocksize for DCT (def: %ix%i) (max: %ix%i)\n", NX_DEF, NY_DEF, NX_MAX, NY_MAX );
	fprintf( msgout, " [-???]   specify output format (see below)\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-img]   pixel data with header (.pgm)\n" );
	fprintf( msgout, " [-dct]   DCT coefficients with header (.dct)\n" );
	fprintf( msgout, " [-rawp]  raw pixel data (.raw)\n" );
	fprintf( msgout, " [-craw?] raw DCT coefficents (0=std,1=dhf,2=squ,3=unc) (.rwc)\n" );
	fprintf( msgout, " [-cimg?] raw DCT coefficents as image (0=std,1=dhf,2=squ,3=unc) (.pgm)\n" );
	fprintf( msgout, " [-csv]   csv file, read in scripting file first (.csv)\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, "Examples: \"trandct -v -o -b4x8 baboon.pgm\"\n" );
	fprintf( msgout, "          \"trandct -b2x2 -cimg2 *.pgm\"\n" );	
}

/* ----------------------- End of main interface functions -------------------------- */

/* ----------------------- Begin of main functions -------------------------- */


/* -----------------------------------------------
	check file and determine filetype
	----------------------------------------------- */
	
bool check_file( void )
{
	FILE* fp;
	
	const char* out_ext = "out";
	const char* tmp_ext = "tmp";
	
	char inp[ 256 ];
	int namelen = strlen( filelist[ file_no ] ) + 5;
	int tmp;
	
	// free memory from filenames if needed
	if ( inpfilename != NULL ) free( inpfilename );
	if ( outfilename != NULL ) free( outfilename );
	if ( tmpfilename != NULL ) free( tmpfilename );
	// alloc memory for filenames
	inpfilename = (char*) calloc( namelen, sizeof( char ) );
	outfilename = (char*) calloc( namelen, sizeof( char ) );
	tmpfilename = (char*) calloc( namelen, sizeof( char ) );
	
	// check if file exists
	if ( access( filelist[ file_no ], 0 ) == -1 ) {
		sprintf( errormessage, "file \"%s\" doesn't exist", filelist[ file_no ]);
		errorlevel = 2;
		return false;
	}
	
	// check if file is readable
	if ( access( filelist[ file_no ], 4 ) == -1 ) {
		sprintf( errormessage, "file \"%s\" is not readable", filelist[ file_no ]);
		errorlevel = 2;
		return false;
	}
	
	// open file, check again
	fp = fopen( filelist[ file_no ], "rb" );
	if ( fp == NULL ) {
		sprintf( errormessage, FRD_ERRMSG, filelist[ file_no ] );
		errorlevel = 2;
		return false;
	}
	
	// check first line
	fgets( inp, 255, fp ); // image id
	
	// check file id, determine filetype
	if ( strncmp( inp, "P5", 2 ) == 0 ) {
		// file is PGM
		in_format = PGM;
		
		// read header
		fgets( inp, 255, fp ); // comment or width/height
		if ( sscanf( inp, "%i %i", &imgwidth, &imgheight ) != 2 ) { // check for comment
			strcpy( pgm_comment, inp ); // copy comment
			fgets( inp, 255, fp );
			if ( sscanf( inp, "%i %i", &imgwidth, &imgheight ) != 2 ) { // read in width/height
				sprintf( errormessage, "error in header structure" );
				errorlevel = 2;
				return false;
			}
		}
		else { // write empty comment
			strcpy( pgm_comment, "" );
		}
		fgets( inp, 255, fp ); // maximum pixel value
		if ( sscanf( inp, "%i", &tmp ) != 1 ) {
			sprintf( errormessage, "error in header structure" );
			errorlevel = 2;
			return false;
		}
		if ( ( tmp >= 256 ) || ( tmp <= 0 ) ) {
			sprintf( errormessage, "bitrate (%i) is not supported", tmp );
			errorlevel = 2;
			return false;
		}
		// set other info
		imgsize = imgwidth * imgheight;
		nx_c 	= 0;
		ny_c    = 0;
		// read in image data
		imgdata = (unsigned char*) calloc( imgsize, sizeof( char ) );
		if ( imgdata == NULL ) {
			sprintf( errormessage, MEM_ERRMSG );
			errorlevel = 2;
			return false;
		}
		if ( fread( imgdata, sizeof( char ), imgsize, fp ) != (unsigned) imgsize ) {
			sprintf( errormessage, "not enough information in file" );
			errorlevel = 2;
			return false;
		}
	}
	else if ( strncmp( inp, "C5", 2 ) == 0 ) {
		// file is DCT coefficients
		in_format = DCT;
		// read header
		fgets( inp, 255, fp ); // comment or width/height
		if ( sscanf( inp, "%i %i", &imgwidth, &imgheight ) != 2 ) { // check for comment
			strcpy( pgm_comment, inp ); // copy comment
			fgets( inp, 255, fp );
			if ( sscanf( inp, "%i %i", &imgwidth, &imgheight ) != 2 ) { // read in width/height
				sprintf( errormessage, "error in header structure" );
				errorlevel = 2;
				return false;
			}
		}
		fgets( inp, 255, fp ); // blocksizes x/y
		if ( sscanf( inp, "%ix%i", &nx_c, &ny_c ) != 2 ) {
			sprintf( errormessage, "error in header structure" );
			errorlevel = 2;
			return false;
		}
		
		// set other info
		imgsize = imgwidth * imgheight;
		cofsize = DCT_IMGS( nx_c, ny_c );
		
		// read in coefficient data
		coefdata = (signed CTYPE*) calloc( cofsize, sizeof( CTYPE ) );
		if ( coefdata == NULL ) {
			sprintf( errormessage, MEM_ERRMSG );
			errorlevel = 2;
			return false;
		}
		if ( fread( coefdata, sizeof( CTYPE ), cofsize, fp ) != (unsigned) cofsize ) {
			sprintf( errormessage, "not enough information in file" );
			errorlevel = 2;
			return false;
		}
	}
	else if ( strncmp( inp, "<TDAS>", 6 ) == 0 ) {
		// file is dependency analysis script
		in_format = TDS;
		if ( !parse_tdasfile( fp ) )
			return false;
	}
	else {
		// file is neither
		in_format = UNK;
		fclose( fp );
		sprintf( errormessage, "filetype of file \"%s\" is unknown", filelist[ file_no ] );
		errorlevel = 2;
		return false;		
	}
	
	// close file
	fclose ( fp );
	
	// create filenames if needed, do some error checks
	if ( in_format != TDS )
	{		
		// error check - can't write csv or plt without scripting file
		if ( ( curr_col == NULL ) && ( out_format == CSV ) ) {
			sprintf( errormessage, "can't write .csv without scripting file" );
			errorlevel = 2;
			return false;
		}
		
		// set default action if needed
		if ( def_actn ) {
			if ( curr_col != NULL )
				out_format = CSV;
			else if ( in_format == PGM )
				out_format = COEF;
			else
				out_format = IMG;
		}
		
		// set extensions
		switch ( out_format )
		{
			case IMG:
				out_ext = "pgm";
				tmp_ext = "tmp";
				break;
				
			case COEF:
				out_ext = "dct";
				tmp_ext = "tmp";
				break;
				
			case RAWP:
				out_ext = "raw";
				tmp_ext = "tmp";
				break;
				
			case CRAW:
				out_ext = "rwc";
				tmp_ext = "tmp";
				break;
				
			case CIMG:
				out_ext = "pgm";
				tmp_ext = "tmp";
				break;
				
			case CSV:
				out_ext = "csv";
				tmp_ext = "plt";
				break;
		}
		
		// copy input filename		
		strcpy( inpfilename, filelist[ file_no ] );
		
		// create main output filename
		set_extension( outfilename, filelist[ file_no ], (char*) out_ext );
		if ( !overwrite ) {
			while ( access( outfilename, 0 ) == 0 ) {
				namelen += sizeof( char );
				outfilename = (char*) realloc( outfilename, namelen );
				add_underscore( outfilename );
			}
		}
		
		// create temp output filename
		set_extension( tmpfilename, filelist[ file_no ], tmp_ext );
		if ( !overwrite ) {
			while ( access( tmpfilename, 0 ) == 0 ) {
				namelen += sizeof( char );
				tmpfilename = (char*) realloc( tmpfilename, namelen );
				add_underscore( tmpfilename );
			}
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	parse dependency analysis script
	----------------------------------------------- */

bool parse_tdasfile( FILE* fp )
{
	TDS_COL* new_col = NULL;
	char inp[ 1024 ];
	char cmd[ 1024 ];
	int ln = 1;
	char* cur;
	
	
	// check for existing scripting data
	if ( curr_col != NULL ) {
		sprintf( errormessage, "only one scripting file is allowed per session" );
		errorlevel = 2;
		return false;
	}
	
	// parse file
	while ( fgets( inp, 1023, fp ) != NULL )
	{
		// reset current pointer, increment line number
		cur = inp;
		ln++;
		// return false & clean up if some error happened
		if ( errorlevel == 2 ) {
			if ( new_col != NULL )
				free( new_col );
			while ( curr_col != NULL ) {
				new_col = curr_col;
				curr_col = curr_col->prev;
				free( new_col );
			}
			return false;
		}
		// parse current line
		while ( ( cur = fetch_command( cur, cmd ) ) != NULL )
		{			
			if ( cmd[ 0 ] == '#' ) {
				// ignore comments
				break;
			}
			else if ( strcmp( cmd, "[inc_incp]" ) == 0 ) {
				// include incomplete lines
				inc_incp = true;			
			}
			else if ( strcmp( cmd, "[cmt_incp]" ) == 0 ) {
				// comment before incomplete lines
				cmt_incp = true;			
			}
			else if ( strcmp( cmd, "[inc_name]" ) == 0 ) {
				// include names in csv
				inc_name = true;			
			}
			else if ( strcmp( cmd, "[inc_chkl]" ) == 0 ) {
				// include checklist
				inc_chkl = true;			
			}
			else if ( strcmp( cmd, "[inc_desc]" ) == 0 ) {
				// include descriptions (first collumn)
				inc_desc = true;			
			}
			else if ( strcmp( cmd, "[chk_ssiz]" ) == 0 ) {
				// check for same size
				chk_ssiz = true;			
			}
			else if ( strcmp( cmd, "[wrt_gplt]" ) == 0 ) {
				// write gnuplot .plt file
				wrt_gplt = true;			
			}
			else if ( strcmp( cmd, "[COL]" ) == 0 ) {
				// start new column
				// check if old column complete
				if ( new_col != NULL ) {
					sprintf( errormessage, "columns have to be closed with '[/COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// alloc new column
				new_col = ( TDS_COL* ) calloc( 1, sizeof( TDS_COL ) );
				if ( new_col == NULL ) {
					sprintf( errormessage, MEM_ERRMSG );
					errorlevel = 2;
					break;
				}
				// set the initial values
				new_col->t_nx    = -1;
				new_col->t_ny    = -1;
				new_col->x_stp   = -1;
				new_col->y_stp   = -1;
				new_col->x_pos   =  0;
				new_col->y_pos   =  0;
				new_col->u_pos   = -1;
				new_col->v_pos   = -1;
				new_col->x_til   =  0;
				new_col->y_til   =  0;
				new_col->u_til   = -1;
				new_col->v_til   = -1;
				new_col->name[0] = '\0';
				new_col->numrows = 0;
				new_col->nummiss = 0;
				new_col->coldata = NULL;
				new_col->chklist = NULL;
				new_col->prev    = curr_col;
				new_col->next    = NULL;
			}
			else if ( strcmp( cmd, "[/COL]" ) == 0 ) {
				// finish column
				// check for completeness
				if ( ( new_col->t_nx > NX_MAX ) || ( new_col->t_ny > NY_MAX ) ) {
					sprintf( errormessage, "dct blocksize %ix%i is not allowed (line #%i)", new_col->t_nx, new_col->t_ny, ln );
					errorlevel = 2;
					break;
				}
				if ( ( new_col->t_nx < 1 ) || ( new_col->t_ny < 1 ) ) {
					sprintf( errormessage, "missing dct blocksize specification, use 'bls ?x?' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				if ( ( new_col->x_stp < 1 ) || ( new_col->y_stp < 1 ) ) {
					sprintf( errormessage, "missing step size or step size not allowed, use 'stp ?/?' (%i/%i) (line #%i)", ln, new_col->x_stp, new_col->y_stp );
					errorlevel = 2;
					break;
				}
				if ( ( new_col->u_pos < 0 ) || ( new_col->v_pos < 0 ) ) {
					sprintf( errormessage, "missing get coordinates, use 'get ?/?/?/?' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				if ( ( new_col->u_pos >= new_col->t_nx ) || ( new_col->v_pos >= new_col->t_ny ) ) {
					sprintf( errormessage, "get coordinates ?/?/%i/%i are not allowed (line #%i)", new_col->u_pos, new_col->v_pos, ln );
					errorlevel = 2;
					break;
				}
				// if this test passed, add the column here
				if ( curr_col != NULL )
					curr_col->next = new_col;
				curr_col = new_col;
				new_col = NULL;
			}
			else if ( strcmp( cmd, "nam" ) == 0 ) { // declare name for column
				// is column declared?
				if ( new_col == NULL ) {
					sprintf( errormessage, "syntax error: declare a new column first with '[COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// get next command
				if ( ( cur = fetch_command( cur, cmd ) ) == NULL ) {
					sprintf( errormessage, "syntax error: missing name specification 'nam \"[NAME]\"' (line #%i)", ln );
					errorlevel = 2;
				}
				// check command
				if ( sscanf( cmd, "%s", new_col->name ) != 1 ) {
					sprintf( errormessage, "syntax error: error in name specification (line #%i)", ln );
					errorlevel = 2;
					break;
				}
			}
			else if ( strcmp( cmd, "bls" ) == 0 ) { // specifiy blocksize
				// is column declared?
				if ( new_col == NULL ) {
					sprintf( errormessage, "syntax error: declare a new column first with '[COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// get next command
				if ( ( cur = fetch_command( cur, cmd ) ) == NULL ) {
					sprintf( errormessage, "syntax error: missing block size specification 'bls ?/?' (line #%i)", ln );
					errorlevel = 2;
				}
				// check command
				if ( sscanf( cmd, "%ix%i", &(new_col->t_nx), &(new_col->t_ny) ) != 2 ) {
					sprintf( errormessage, "syntax error: error in block size specification (line #%i)", ln );
					errorlevel = 2;
					break;
				}
			}
			else if ( strcmp( cmd, "stp" ) == 0 ) { // specifiy stepsize
				// is column declared?
				if ( new_col == NULL ) {
					sprintf( errormessage, "syntax error: declare a new column first with '[COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// get next command
				if ( ( cur = fetch_command( cur, cmd ) ) == NULL ) {
					sprintf( errormessage, "syntax error: missing step size specification 'stp ?/?' (line #%i)", ln );
					errorlevel = 2;
				}
				// check command
				if ( sscanf( cmd, "%i/%i", &(new_col->x_stp), &(new_col->y_stp) ) != 2 ) {
					sprintf( errormessage, "syntax error: error in step size specification (line #%i)", ln );
					errorlevel = 2;
					break;
				}
			}
			else if ( strcmp( cmd, "get" ) == 0 ) { // specifiy get coordinates
				// is column declared?
				if ( new_col == NULL ) {
					sprintf( errormessage, "syntax error: declare a new column first with '[COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// get next command
				if ( ( cur = fetch_command( cur, cmd ) ) == NULL ) {
					sprintf( errormessage, "syntax error: missing get coordinates specification 'get ?/?/?/?' (line #%i)", ln );
					errorlevel = 2;
				}
				// check command
				if ( sscanf( cmd, "%i/%i/%i/%i", &(new_col->x_pos), &(new_col->y_pos), &(new_col->u_pos), &(new_col->v_pos) ) != 4 ) {
					sprintf( errormessage, "syntax error: error in get coordinates specification (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// set "till" parameters
				new_col->x_til = new_col->x_pos;
				new_col->y_til = new_col->y_pos;
				new_col->u_til = new_col->u_pos;
				new_col->v_til = new_col->v_pos;
			}
			else if ( strcmp( cmd, "til" ) == 0 ) { // specifiy get coordinates
				// is column declared?
				if ( new_col == NULL ) {
					sprintf( errormessage, "syntax error: declare a new column first with '[COL]' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				if ( ( new_col->u_til < 0 ) || ( new_col->v_til < 0 ) ) {
					sprintf( errormessage, "syntax error: 'til ?/?/?/?' before 'get ?/?/?/?' (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// get next command
				if ( ( cur = fetch_command( cur, cmd ) ) == NULL ) {
					sprintf( errormessage, "syntax error: missing get coordinates specification 'til ?/?/?/?' (line #%i)", ln );
					errorlevel = 2;
				}
				// check command
				if ( sscanf( cmd, "%i/%i/%i/%i", &(new_col->x_til), &(new_col->y_til), &(new_col->u_til), &(new_col->v_til) ) != 4 ) {
					sprintf( errormessage, "syntax error: error in get coordinates specification (line #%i)", ln );
					errorlevel = 2;
					break;
				}
				// check back with "pos" parameters
				if ( ( new_col->x_til < new_col->x_pos ) || ( new_col->y_til < new_col->y_pos ) ||
					 ( new_col->u_til < new_col->u_pos ) || ( new_col->v_til < new_col->v_pos ) ) {
					sprintf( errormessage, "inconsistency in 'get'/'til' parameters (line #%i)", ln );
					errorlevel = 2;
					break;
				}
			}
			else { // syntax error / unknown command
				sprintf( errormessage, "syntax error: unknown command: \"%s\" (line #%i)", cmd, ln );
				errorlevel = 2;
				break;
			}
		}
	}
	
	// check for current column
	if ( curr_col == NULL ) {
		sprintf( errormessage, "file contains no columns data" );
		errorlevel = 2;
		return false;
	}
	
	// go back to first column	
	while ( curr_col->prev != NULL )
		curr_col = curr_col->prev;
	
	
	return true;
}


/* -----------------------------------------------
	transform image to DCT coefficients
	----------------------------------------------- */
bool trans_image( void )
{
	signed CTYPE  F[ NX_MAX * NY_MAX ]; // block (coeffs)
	unsigned char f[ NX_MAX * NY_MAX ]; // block (image)
	int iu, iv;
	int x, y;
	int dp;
	int i;
	
	
	// prepare DCT trasnformation
	if ( !prepare_dct( nx_s, ny_s ) ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// alloc memory
	if ( coefdata != NULL ) free( coefdata );
	coefdata = (signed CTYPE*) calloc( cofsize, sizeof( CTYPE ) );
	if ( coefdata == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// convert to DCT coefficients	
	for ( y = 0; y < imgheight; y += ny_c ) {
		for ( x = 0; x < imgwidth; x += nx_c ) {
			// fill current block
			for ( iv = 0; iv < ny_c; iv++ ) {
				for ( iu = 0; iu < nx_c; iu++ ) {
					if ( ( ( x + iu ) < imgwidth ) && ( ( y + iv ) < imgheight ) ) {						
						f[ iu + ( nx_c * iv ) ] = imgdata[ ( x + iu ) + ( imgwidth * ( y + iv ) ) ];
					}
					else
						f[ iu + ( nx_c * iv ) ] = 128;
				}
			}
			// transform block
			block_fdct_2d( f, F );
			// copy to DCT coeffs
			dp = ( ( x / nx_c ) + ( bx_c * ( y / ny_c ) ) ) * ns_c;
			for ( i = 0; i < ns_c; i++ )
				coefdata[ dp + i ] = F[ i ];
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	transform DCT coefficients to image
	----------------------------------------------- */
bool trans_coeffs( void )
{
	signed CTYPE  F[ NX_MAX * NY_MAX ]; // block (coeffs)
	unsigned char f[ NX_MAX * NY_MAX ]; // block (image)
	int iu, iv;
	int x, y;
	int dp;
	int i;
	
	
	// prepare DCT trasnformation
	if ( !prepare_dct ( nx_c, ny_c ) ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// alloc memory
	if ( imgdata != NULL ) free( imgdata );
	imgdata = (unsigned char*) calloc( imgsize, sizeof( char ) );
	if ( imgdata == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// convert to pixel data
	for ( y = 0; y < imgheight; y += ny_c ) {
		for ( x = 0; x < imgwidth; x += nx_c ) {
			// fill current block
			dp = ( ( x / nx_c ) + ( bx_c * ( y / ny_c ) ) ) * ns_c;
			for ( i = 0; i < ns_c; i++ )
				F[ i ] = coefdata[ dp + i ];
			// transform block
			block_idct_2d( F, f );
			// copy to image data
			for ( iv = 0; iv < ny_c; iv++ ) {
				for ( iu = 0; iu < nx_c; iu++ ) {
					if ( ( ( x + iu ) < imgwidth ) && ( ( y + iv ) < imgheight ) )
						imgdata[ ( x + iu ) + ( imgwidth * ( y + iv ) ) ] = f[ iu + ( nx_c * iv ) ];
				}
			}
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	write output file
	----------------------------------------------- */
bool write_file( void )
{
	FILE* fp;
	
	signed CTYPE coef;
	unsigned char pix;
	int bp, dp;
	int hp, vp;
	
	// open file for output
	fp = fopen( outfilename, "wb" );
	if ( fp == NULL ){
		sprintf( errormessage, FWR_ERRMSG, outfilename );
		errorlevel = 2;
		return false;
	}
	
	// write output file
	switch ( out_format )
	{			
		case IMG:
			// write pgm header
			fprintf( fp, "P5\n%i %i\n255\n", imgwidth, imgheight );
			// write image data
			fwrite( imgdata, sizeof( char ), imgsize, fp );
			break;
			
		case COEF:
			// write dct file header
			fprintf( fp, "C5\n%i %i\n%ix%i\n", imgwidth, imgheight, nx_c, ny_c );
			// write coefficient data
			fwrite( coefdata, sizeof( CTYPE ), cofsize, fp );
			break;
			
		case RAWP:
			// write image data
			fwrite( imgdata, sizeof( char ), imgsize, fp );
			break;
			
		case CRAW:
		case CIMG:			
			switch ( outp_mode )
			{
				case 0:
					// simple block per block 'dhuf'
					if ( out_format == CIMG ) { // write pgm header
						fprintf( fp, "P5\n%i %i\n255\n", ns_c, cofsize / ns_c );
					}
					// actual conversion starts here
					for ( dp = 0; dp < cofsize; dp++ ) {
						coef = coefdata[ dp ];
						WR_C ( fp, coef, pix );
					}
					break;
					
				case 1:
					// simple band per band 'collection'
					if ( out_format == CIMG ) { // write pgm header
						fprintf( fp, "P5\n%i %i\n255\n", bx_c, cofsize / bx_c );
					}
					// actual conversion starts here
					for ( bp = 0; bp < ns_c; bp++ ) {
						for ( dp = 0; dp < cofsize; dp += ns_c ) {
							coef = coefdata[ dp + bp ];
							WR_C ( fp, coef, pix );
						}
					}
					break;
					
				case 2:
					// band per band schematic 'square collection'
					if ( out_format == CIMG ) { // write pgm header
						fprintf( fp, "P5\n%i %i\n255\n", cofwidth, cofheight );
					}
					// actual conversion starts here
					dp = 0;
					for ( bp = 0; bp < ns_c; ) {
						for ( hp = 0; hp < ( cofwidth * ny_c ); hp += ns_c ) {
							coef = coefdata[ dp + bp + hp ];
							WR_C ( fp, coef, pix );
						}
						if ( ( ++bp % nx_c ) == 0 ) {
							dp += ( cofwidth * ny_c );
							if ( dp >= cofsize )
								dp = 0;
							else
								bp -= nx_c;
						}
					}
					break;
					
				case 3:
					// block per block schematic 'uncollection'
					if ( out_format == CIMG ) { // write pgm header
						fprintf( fp, "P5\n%i %i\n255\n", cofwidth, cofheight );
					}
					// actual conversion starts here
					for ( vp = 0; vp < cofheight; vp++ ) {
						for ( hp = 0; hp < cofwidth; hp++ ) {
							bp = ( ( vp % ny_c ) * nx_c ) + ( hp % nx_c );
							dp = ( ( vp / ny_c ) * ( cofwidth * ny_c ) ) + ( ( hp / nx_c ) * ns_c );
							coef = coefdata[ dp + bp ];
							WR_C ( fp, coef, pix );
						}
					}
					break;
			}			
			break;
			
		case CSV:
			// this function shouldn't end up here - so it's an error
			sprintf( errormessage, "severe program logic error" );
			errorlevel = 2;
			return false;
			
	}
	
	// check for errors
	if ( ferror( fp ) ) {
		fclose( fp );
		sprintf( errormessage, "write error, possibly drive is full" );
		errorlevel = 2;		
		return false;
	}
	
	// close file
	fclose( fp );
	
	
	return true;
}


/* -----------------------------------------------
	read one column (tdas)
	----------------------------------------------- */
bool read_column( void )
{	
	signed CTYPE*  coldata;
	unsigned char* chklist;
	int numrows = 0;
	int nummiss = 0;
	
	int x_stp = curr_col->x_stp;
	int y_stp = curr_col->y_stp;
	int x_pos = curr_col->x_pos;
	int y_pos = curr_col->y_pos;
	int u_pos = curr_col->u_pos;
	int v_pos = curr_col->v_pos;
	int x_til = curr_col->x_til;
	int y_til = curr_col->y_til;
	int u_til = curr_col->u_til;
	int v_til = curr_col->v_til;
	
	int ix, iy;
	int iu, iv;
	int cf, ch;
	
	int ci = 0;
	int x, y;
	
	
	// calculate column size
	numrows = ( ( bx_c + ( x_stp - 1 ) ) / x_stp ) * ( ( by_c + ( y_stp - 1 ) ) / y_stp );
	// allocate memory
	coldata = (signed CTYPE*) calloc( numrows, sizeof( CTYPE ) );
	if ( coldata == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	chklist = (unsigned char*) calloc( numrows, sizeof( char ) );
	if ( chklist == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// calculate correction factor
	cf = ( x_til - x_pos + 1 ) * ( y_til - y_pos + 1 ) * ( u_til - u_pos + 1 ) * ( v_til - v_pos + 1 );
	ch = cf / 2;
	
	// go through coefficient data, fill column
	for ( y = 0; y < by_c; y += y_stp ) {
		for ( x = 0; x < bx_c; x += x_stp ) {
			coldata[ ci ] = 0;
			// treatment for image edges
			if ( ( ( x + x_pos ) < 0 ) || ( ( x + x_pos ) >= bx_c ) ||
				 ( ( y + y_pos ) < 0 ) || ( ( y + y_pos ) >= by_c ) ||
				 ( ( x + x_til ) < 0 ) || ( ( x + x_til ) >= bx_c ) ||
				 ( ( y + y_til ) < 0 ) || ( ( y + y_til ) >= by_c ) ) {				
				chklist[ ci ] = 0;
				nummiss++;
			}
			else {
				// accumulate coefficients
				for ( ix = x + x_pos; ix <= x + x_til; ix++ )
				for ( iy = y + y_pos; iy <= y + y_til; iy++ )
				for ( iu = u_pos; iu <= u_til; iu++ )
				for ( iv = v_pos; iv <= v_til; iv++ )
					coldata[ ci ] += coefdata[ CFPOS( ix, iy, iu, iv ) ];
				// calculate average
				coldata[ ci ] = ( coldata[ ci ] > 0 ) ? ( coldata[ ci ] + ch ) / cf : ( coldata[ ci ] - ch ) / cf;
				chklist[ ci ] = 1;
			}
			ci++;			
		}
	}
	
	// put coldata / chklist / numrows in place
	curr_col->coldata = coldata;
	curr_col->chklist = chklist;
	curr_col->numrows = numrows;
	curr_col->nummiss = nummiss;
	
	
	return true;
}


/* -----------------------------------------------
	write csv file
	----------------------------------------------- */
	
bool write_csv( void )
{
	FILE* fp;
	char tmp[256];
	
	signed CTYPE** coldata;
	unsigned char** chklist;
	char** names;
	int *t_nx, *t_ny;
	int *x_stp, *y_stp;
	int *x_pos, *y_pos;
	int *u_pos, *v_pos;
	int *x_til, *y_til;
	int *u_til, *v_til;
	int *numrows;
	int *nummiss;
	
	unsigned char* global_chklist;
	int nummiss_glb = 0;
	int numrows_max = 0;
	int numcols = 0;
	int col, row;
	
	// go back to first column	
	while ( curr_col->prev != NULL )
		curr_col = curr_col->prev;
	
	// count columns
	while ( true ) {		
		numcols++;
		numrows_max = ( curr_col->numrows > numrows_max ) ? curr_col->numrows : numrows_max;
		if ( curr_col->next != NULL )
			curr_col = curr_col->next;
		else break;
	}
	
	// allocate all the memory needed
	coldata = new signed CTYPE*[ numcols ];
	chklist = new unsigned char*[ numcols ];
	names = new char*[ numcols ];
	t_nx  = new int[ numcols ];
	t_ny  = new int[ numcols ];
	x_stp = new int[ numcols ];;
	y_stp = new int[ numcols ];
	x_pos = new int[ numcols ];
	y_pos = new int[ numcols ];
	u_pos = new int[ numcols ];
	v_pos = new int[ numcols ];
	x_til = new int[ numcols ];
	y_til = new int[ numcols ];
	u_til = new int[ numcols ];
	v_til = new int[ numcols ];
	numrows = new int[ numcols ];
	nummiss = new int[ numcols ];
	
	// allocate memory for global checklist (could get big)
	global_chklist = ( unsigned char* ) calloc ( numrows_max, sizeof( char ) );
	if ( global_chklist == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// go back to first column	
	while ( curr_col->prev != NULL )
		curr_col = curr_col->prev;
	
	// collect data
	for ( col = 0; col < numcols; col++ ) {
		coldata[ col ] = curr_col->coldata;
		chklist[ col ] = curr_col->chklist;
		names[ col ] = curr_col->name;
		t_nx[ col ]  = curr_col->t_nx;
		t_ny[ col ]  = curr_col->t_ny;
		x_stp[ col ] = curr_col->x_stp;
		y_stp[ col ] = curr_col->y_stp;
		x_pos[ col ] = curr_col->x_pos;
		y_pos[ col ] = curr_col->y_pos;
		u_pos[ col ] = curr_col->u_pos;
		v_pos[ col ] = curr_col->v_pos;
		x_til[ col ] = curr_col->x_til;
		y_til[ col ] = curr_col->y_til;
		u_til[ col ] = curr_col->u_til;
		v_til[ col ] = curr_col->v_til;
		numrows[ col ] = curr_col->numrows;
		nummiss[ col ] = curr_col->nummiss;
		if ( curr_col->next != NULL )
			curr_col = curr_col->next;
		else break;
	}
	
	// calculate global checklist
	for ( row = 0; row < numrows_max; row++ ) {
		global_chklist[ row ] = 1;
		for ( col = 0; col < numcols; col++ ) {
			if ( numrows[ col ] > row )
				global_chklist[ row ] &= chklist[ col ][ row ];
		}
		if ( global_chklist[ row ] == 0 )
			nummiss_glb++;
	}
	
	// check for same size if needed
	if ( chk_ssiz ) {
		for ( col = 0; col < numcols; col++ ) {
			if ( numrows_max != numrows[ col ] ) {
				sprintf( errormessage, "numbers of rows don't match" );
				errorlevel = 2;
				return false;
			}
		}
	}
	
	// open csv file for output
	fp = fopen( outfilename, "wb" );
	if ( fp == NULL ){
		sprintf( errormessage, FWR_ERRMSG, outfilename );
		errorlevel = 2;
		return false;
	}
	
	// write column numbers
	fprintf( fp, "# " );
	if ( inc_desc )	{
		fprintf( fp, "$column%c ", CSV_SEP );
		for ( col = 0; col < numcols; col++ )
			fprintf( fp, "$%i%c ", col + 2, CSV_SEP );
	}
	else {
		for ( col = 0; col < numcols; col++ )
			fprintf( fp, "$%i%c ", col + 1, CSV_SEP );
	}
	fprintf( fp, "\n" );
	
	// write names if needed
	if ( inc_name ) {
		fprintf( fp, "# " );
		if ( inc_desc )	fprintf( fp, "name%c ", CSV_SEP );
		for ( col = 0; col < numcols; col++ )
			fprintf( fp, "'%s'%c ", names[ col ], CSV_SEP );
		fprintf( fp, "\n" );
	}
	
	// write block sizes
	fprintf( fp, "# " );
	if ( inc_desc )	fprintf( fp, "block x/y%c ", CSV_SEP );
	for ( col = 0; col < numcols; col++ )
		fprintf( fp, "'%ix%i'%c ", t_nx[ col ], t_ny[ col ], CSV_SEP );
	fprintf( fp, "\n" );
	
	// write step sizes
	fprintf( fp, "# " );
	if ( inc_desc )	fprintf( fp, "step x/y%c ", CSV_SEP );
	for ( col = 0; col < numcols; col++ )
		fprintf( fp, "'(%i/%i)'%c ", x_stp[ col ], y_stp[ col ], CSV_SEP );
	fprintf( fp, "\n" );
	
	// write get position x/y
	fprintf( fp, "# " );
	if ( inc_desc )	fprintf( fp, "pos x/y (rel)%c ", CSV_SEP );
	for ( col = 0; col < numcols; col++ ) {
		if ( ( x_pos[ col ] == x_til[ col ] ) && ( y_pos[ col ] == y_til[ col] ) )
			fprintf( fp, "'(%i/%i)'%c ", x_pos[ col ], y_pos[ col ], CSV_SEP );
		else
			fprintf( fp, "'(%i/%i)...(%i/%i)'%c ", x_pos[ col ], y_pos[ col ], x_til[ col ], y_til[ col ], CSV_SEP );
	}
	fprintf( fp, "\n" );
	
	// write get position u/v
	fprintf( fp, "# " );
	if ( inc_desc )	fprintf( fp, "pos u/v (abs)%c ", CSV_SEP );
	for ( col = 0; col < numcols; col++ ) {
		if ( ( u_pos[ col ] == u_til[ col ] ) && ( v_pos[ col ] == v_til[ col] ) )
			fprintf( fp, "'(%i/%i)'%c ", u_pos[ col ], v_pos[ col ], CSV_SEP );
		else
			fprintf( fp, "'(%i/%i)...(%i/%i)'%c ", u_pos[ col ], v_pos[ col ], u_til[ col ], v_til[ col ], CSV_SEP );
	}
	fprintf( fp, "\n" );
	
	// write number of lines
	fprintf( fp, "# " );
	if ( inc_desc )	fprintf( fp, "num elements%c ", CSV_SEP );
	for ( col = 0; col < numcols; col++ )
		fprintf( fp, "'%i (-%i)'%c ", numrows[ col ], nummiss[ col ], CSV_SEP );
	if ( inc_chkl ) fprintf( fp, "'%i (-%i)'%c ", numrows_max, nummiss_glb, CSV_SEP );
	fprintf( fp, "\n" );
	
	// write empty line
	fprintf( fp, "\n" );
	
	// begin writing data
	for ( row = 0; row < numrows_max; row++ ) {
		// skip incomplete lines / comment out if needed
		if ( !global_chklist[ row ] ) {
			if ( !inc_incp )
				continue;
			else if ( cmt_incp )
				fprintf( fp, "# " );
		}
		if ( inc_desc )	fprintf( fp, "ln%i%c ", row, CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			if ( ( numrows[ col ] <= row ) || ( chklist[ col ][ row ] == 0 ) ) {
				fprintf( fp, "%c%c ", MISS_REP, CSV_SEP );
			}
			else fprintf( fp, "%i%c ", coldata[ col ][ row ], CSV_SEP );
		}
		if ( inc_chkl ) fprintf( fp, "%i%c ", global_chklist[ row ], CSV_SEP );
		fprintf( fp, "\n" );
	}
	
	// close file
	fclose( fp );
	
	if ( wrt_gplt ) // this code is for writing a gnuplot .plt file
	{
		// open pltfile for output
		fp = fopen( tmpfilename, "wb" );
		if ( fp == NULL ){
			sprintf( errormessage, FWR_ERRMSG, tmpfilename );
			errorlevel = 2;
			return false;
		}
		
		// write headline
		fprintf( fp, "# gnuplot file for analysis of file \"%s\" written by tranDCT v%i.%i\n", inpfilename, prgversion / 10, prgversion % 10  );
		fprintf( fp, "# this might need some editing\n" );
		fprintf( fp, "\n" );
		
		// write gnuplot commands
		fprintf( fp, "set grid\n" );
		fprintf( fp, "set nokey\n" );
		fprintf( fp, "\n" );
		fprintf( fp, "set title \"trandDCT coefficient dependency analysis of <%s>\"\n", inpfilename );
		if ( !inc_name ) {
			fprintf( fp, "set xlabel \"x = %i; y = %i; u = %i; v = %i \"\n", x_pos[ 0 ], y_pos[ 0 ], u_pos[ 0 ], v_pos[ 0 ] );
			fprintf( fp, "set ylabel \"x = %i; y = %i; u = %i; v = %i \"\n", x_pos[ 1 ], y_pos[ 1 ], u_pos[ 1 ], v_pos[ 1 ] );
		}
		else if ( ( names[ 0 ][ 0 ] == '"' ) && ( names[ 1 ][ 0 ] == '"' ) ) {
			fprintf( fp, "set xlabel %s\n", names[ 0 ] );
			fprintf( fp, "set ylabel %s\n", names[ 1 ] );
		}
		else {
			fprintf( fp, "set xlabel \"%s\"\n", names[ 0 ] );
			fprintf( fp, "set ylabel \"%s\"\n", names[ 1 ] );
		}
		fprintf( fp, "set size square\n" );
		fprintf( fp, "set zeroaxis -1\n" );
		fprintf( fp, "\n" );
		
		// calculate a best-fit line
		fprintf( fp, "g(x) = a*x + b\n" );
		fprintf( fp, "fit g(x) '%s' using %s via a, b\n", outfilename, ( inc_desc ) ? "($2):($3)" : "($1):($2)" );
		
		// this is the actual plot command
		fprintf( fp, "plot '%s' using %s with dots lw 3, g(x) lw 3\n", outfilename, ( inc_desc ) ? "($2):($3)" : "($1):($2)" );
		fprintf( fp, "\n" );
		
		// the pause command - needed to actually see the plot on screen
		fprintf( fp, "pause -1\n" );
		fprintf( fp, "\n" );
		
		// some info about columns in the dat file
		// write description
		fprintf( fp, "# the following data is found in \"%s\" and can be used for this plot:\n", outfilename );
		
		// write column numbers
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "$column", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "$%i", col + ( ( inc_desc ) ? 2 : 1 ) );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write names if needed
		if ( inc_name ) {
			fprintf( fp, "# " );
			fprintf( fp, "%15s%c ", "name", CSV_SEP );
			for ( col = 0; col < numcols; col++ )
				fprintf( fp, "%9s%c ", names[ col ], CSV_SEP );
			fprintf( fp, "\n" );
		}
		
		// write block sizes
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "block x/y", CSV_SEP );		
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "%ix%i", t_nx[ col ], t_ny[ col ] ); 
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write step sizes
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "step x/y", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "(%i/%i)", x_stp[ col ], y_stp[ col ] );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write get position x/y
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "pos x/y (rel)", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {			
			sprintf( tmp, "(%i/%i)", x_pos[ col ], y_pos[ col ] );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write get position u/v
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "pos u/v (abs)", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "(%i/%i)", u_pos[ col ], v_pos[ col ] );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write number of elements
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "num elements", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "%i", numrows[ col ] );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// write number of missing elements
		fprintf( fp, "# " );
		fprintf( fp, "%15s%c ", "num missing", CSV_SEP );
		for ( col = 0; col < numcols; col++ ) {
			sprintf( tmp, "%i", nummiss[ col ] );
			fprintf( fp, "%9s%c ", tmp, CSV_SEP );
		}
		fprintf( fp, "\n" );
		
		// close file
		fclose( fp );
	}
		
	// free memory
	free( global_chklist );
	
	return true;
}


/* -----------------------------------------------
	set each variable to its initial value
	----------------------------------------------- */

bool reset_buffers( void )
{
	// free memory from allocated data
	if ( imgdata  != NULL ) free ( imgdata  );
	if ( coefdata != NULL ) free ( coefdata );
	
	// free memory from precalculated tables if needed
	if ( icos_idct_nx != NULL ) free ( icos_idct_nx );
	if ( icos_idct_ny != NULL ) free ( icos_idct_ny );
	if ( icos_fdct_nx != NULL ) free ( icos_fdct_nx );
	if ( icos_fdct_ny != NULL ) free ( icos_fdct_ny );
	
	// free memory from fast tables if needed
	if ( icos_idct_fst != NULL ) free ( icos_idct_fst );
	if ( icos_fdct_fst != NULL ) free ( icos_fdct_fst );
	
	// set pointers NULL
	imgdata  = NULL;
	coefdata = NULL;
	
	// set pointers NULL
	icos_idct_nx = NULL;
	icos_idct_ny = NULL;
	icos_fdct_nx = NULL;
	icos_fdct_ny = NULL;
	
	// set pointer NULL
	icos_idct_fst = NULL;
	icos_fdct_fst = NULL;
	
	// reset variables
	imgsize   = 0;		// size of image
	imgwidth  = 0;		// width of image
	imgheight = 0;		// height of image

	cofsize   = 0;		// size of coefficient image;
	cofwidth  = 0;		// width of coefficient image
	cofheight = 0;		// heigt of coefficient image

	bs_c      = 0;		// block count;
	bx_c      = 0;		// block count x (input)
	by_c      = 0;		// block count y (input)

	ns_c      = 0;		// dct block size
	nx_c      = 0;		// dct block size x (input)
	ny_c      = 0;		// dct block size y (input)
	
	// partial clean up analysis data if needed
	if ( curr_col != NULL ) {
		// go back to first column
		while ( curr_col->prev != NULL )
			curr_col = curr_col->prev;
		// clean up
		while ( true ) {
			if ( curr_col->coldata != NULL ) free ( curr_col->coldata );
			if ( curr_col->chklist != NULL ) free ( curr_col->chklist );
			curr_col->coldata = NULL;
			curr_col->chklist = NULL;
			curr_col->numrows = 0;
			curr_col->nummiss = 0;
			// on to next column
			if ( curr_col->next != NULL )
				curr_col = curr_col->next;
			else break;
		}
	}
	
	return true;
}

/* ----------------------- End of main functions -------------------------- */

/* ----------------------- Begin ofDCT specific functions -------------------------- */


/* -----------------------------------------------
	precalculate some values for FDCT/IDCT
	----------------------------------------------- */
bool prepare_dct( int nx, int ny )
{
	int iu, iv;
	int iy, ix;
	int sx, sy;
	int su, sv;
	int ti;
	
	
	// set info variables
	nx_c      = nx;
	ny_c      = ny;
	ns_c      = nx_c * ny_c;
	bx_c      = DCT_BCH( nx_c );
	by_c      = DCT_BCV( ny_c );
	bs_c      = DCT_BCS( nx_c, ny_c );
	cofwidth  = DCT_IMGW( nx_c );
	cofheight = DCT_IMGH( ny_c );
	cofsize   = DCT_IMGS( nx_c, ny_c );
	
	// decide which idct/fdct function to use
	idct_2d = ( nx_c > ny_c ) ? idct_2d_bnx : idct_2d_bny;
	fdct_2d = ( nx_c > ny_c ) ? fdct_2d_bnx : fdct_2d_bny;
	
	// free memory from precalculated tables if needed
	if ( icos_idct_nx != NULL ) free ( icos_idct_nx );
	if ( icos_idct_ny != NULL ) free ( icos_idct_ny );
	if ( icos_fdct_nx != NULL ) free ( icos_fdct_nx );
	if ( icos_fdct_ny != NULL ) free ( icos_fdct_ny );
	
	// alloc memory for precalculated tables
	icos_idct_nx = (float*) calloc( nx_c * nx_c, sizeof( float ) );
	icos_idct_ny = (float*) calloc( ny_c * ny_c, sizeof( float ) );
	icos_fdct_nx = (float*) calloc( nx_c * nx_c, sizeof( float ) );
	icos_fdct_ny = (float*) calloc( ny_c * ny_c, sizeof( float ) );	
	
	// check for out of memory
	if ( ( icos_idct_nx == NULL ) || ( icos_idct_ny == NULL ) ||
		 ( icos_fdct_nx == NULL ) || ( icos_fdct_ny == NULL ) ) 
	{
		return false;
	}
	
	// precalculate tables
	// idct / nx table
	for ( ix = 0; ix < nx_c; ix++ ) {
		for ( iu = 0; iu < nx_c; iu++ ) {
			icos_idct_nx [ iu + ( ix * nx_c ) ] = ( C_DCT ( iu ) * COS_DCT( ix, iu, nx_c ) ) / DCT_SCALE;
		}
	}
	
	// idct / ny table
	for ( iy = 0; iy < ny_c; iy++ ) {
		for ( iv = 0; iv < ny_c; iv++ ) {
			icos_idct_ny [ iv + ( iy * ny_c ) ] = ( C_DCT ( iv ) * COS_DCT( iy, iv, ny_c ) ) / DCT_SCALE;
		}
	}
	
	// fdct / nx table
	for ( iu = 0; iu < nx_c; iu++ ) {
		for ( ix = 0; ix < nx_c; ix++ ) {
			icos_fdct_nx [ ix + ( iu * nx_c ) ] = ( C_DCT ( iu ) * COS_DCT( ix, iu, nx_c ) * DCT_SCALE ) / nx_c;
		}
	}
	
	// fdct / ny table
	for ( iv = 0; iv < ny_c; iv++ ) {
		for ( iy = 0; iy < ny_c; iy++ ) {
			icos_fdct_ny [ iy + ( iv * ny_c ) ] = ( C_DCT ( iv ) * COS_DCT( iy, iv, ny_c ) * DCT_SCALE ) / ny_c;
		}
	}
	
	// precalculation of fast DCT tables...
	if ( FAST_DCT )
	{
		// free memory if needed
		if ( icos_idct_fst != NULL ) free ( icos_idct_fst );
		if ( icos_fdct_fst != NULL ) free ( icos_fdct_fst );
		
		// alloc memory
		icos_idct_fst = (float*) calloc( nx_c * nx_c * ny_c * ny_c, sizeof( float ) );
		icos_fdct_fst = (float*) calloc( nx_c * nx_c * ny_c * ny_c, sizeof( float ) );
		
		// check for out of memory
		if ( ( icos_idct_fst == NULL ) || ( icos_fdct_fst == NULL ) )
			return false;		
		
		// idct / fast table
		ti = 0;
		for ( iy = 0; iy < ny_c; iy++ ) {
			for ( ix = 0; ix < nx_c; ix++ ) {
				sx = ( ix * nx_c );
				sy = ( iy * ny_c );
				for ( iv = 0; iv < ny_c; iv++ ) {
					for ( iu = 0; iu < nx_c; iu++ ) {
						icos_idct_fst[ ti++ ] = icos_idct_nx[ sx + iu ] * icos_idct_ny[ sy + iv ]; 
					}
				}
			}
		}
		
		// fdct / fast table
		ti = 0;
		for ( iv = 0; iv < ny_c; iv++ ) {
			for ( iu = 0; iu < nx_c; iu++ ) {
				su = ( iu * nx_c );
				sv = ( iv * ny_c );
				for ( iy = 0; iy < ny_c; iy++ ) {
					for ( ix = 0; ix < nx_c; ix++ ) {
						icos_fdct_fst[ ti++ ] = icos_fdct_nx[ su + ix ] * icos_fdct_ny[ sv + iy ]; 
					}
				}
			}
		}
		
		// use fast idct/fdct instead of the standard ones
		idct_2d = idct_2d_fst;
		fdct_2d = fdct_2d_fst;
	}
	
	
	return true;
}


/* -----------------------------------------------
	inverse DCT transform using precalc tables (bigger nx)
	----------------------------------------------- */
float idct_2d_bnx( signed CTYPE* F, int ix, int iy )
{
	float idctlin, idct;
	int sx, sy;
	int tx, ty;
	int iu, iv;
	
	
	// calculate table indexes
	sx = ( ix * nx_c );
	sy = ( iy * ny_c );
	
	// begin transform
	ty = sy;
	idct = 0;
	for ( iv = 0; iv < ny_c; iv++ ) {
		tx = sx;
		idctlin = 0;
		for ( iu = 0; iu < nx_c; iu++ ) {
			idctlin = idctlin + F[ iu + ( nx_c * iv ) ] * icos_idct_nx[ tx++ ];
		}
		idct = idct + idctlin * icos_idct_ny[ ty++ ];
	}

	
	return idct;
}


/* -----------------------------------------------
	inverse DCT transform using precalc tables (bigger ny)
	----------------------------------------------- */
float idct_2d_bny( signed CTYPE* F, int ix, int iy )
{
	float idctlin, idct;
	int sx, sy;
	int tx, ty;
	int iu, iv;
	
	// calculate table indexes
	sx = ( ix * nx_c );
	sy = ( iy * ny_c );
	
	// begin transform
	tx = sx;
	idct = 0;
	for ( iu = 0; iu < nx_c; iu++ ) {
		ty = sy;
		idctlin = 0;
		for ( iv = 0; iv < ny_c; iv++ ) {
			idctlin = idctlin + F[ iu + ( nx_c * iv ) ] * icos_idct_ny[ ty++ ];
		}
		idct = idct + idctlin * icos_idct_nx[ tx++ ];
	}

	
	return idct;
}


/* -----------------------------------------------
	forward DCT transform using precalc tables (bigger nx)
	----------------------------------------------- */
float fdct_2d_bnx( unsigned char* f, int iu, int iv )
{
	float fdctlin, fdct;
	int su, sv;
	int tu, tv;
	int ix, iy;
	
	
	// calculate table indexes
	su = ( iu * nx_c );
	sv = ( iv * ny_c );
	
	// begin transform
	tv = sv;
	fdct = 0;
	for ( iy = 0; iy < ny_c; iy++ ) {
		tu = su;
		fdctlin = 0;
		for ( ix = 0; ix < nx_c; ix++ ) {
			fdctlin = fdctlin + f[ ix + ( nx_c * iy ) ] * icos_fdct_nx[ tu++ ];
		}
		fdct = fdct + fdctlin * icos_fdct_ny[ tv++ ];
	}
	
	
	return fdct;
}


/* -----------------------------------------------
	forward DCT transform using precalc tables (bigger ny)
	----------------------------------------------- */
float fdct_2d_bny( unsigned char* f, int iu, int iv )
{
	float fdctlin, fdct;
	int su, sv;
	int tu, tv;
	int ix, iy;
	
	// calculate table indexes
	su = ( iu * nx_c );
	sv = ( iv * ny_c );
	
	// begin transform
	tu = su;
	fdct = 0;
	for ( ix = 0; ix < nx_c; ix++ ) {
		tv = sv;
		fdctlin = 0;
		for ( iy = 0; iy < ny_c; iy++ ) {
			fdctlin = fdctlin + f[ ix + ( nx_c * iy ) ] * icos_fdct_ny[ tv++ ];
		}
		fdct = fdct + fdctlin * icos_fdct_nx[ tu++ ];
	}
	
	
	return fdct;
}


/* -----------------------------------------------
	inverse DCT transform using precalc tables (fast)
	----------------------------------------------- */
float idct_2d_fst( signed CTYPE* F, int ix, int iy )
{
	float idct;
	int ixy;
	int i;
	
	
	// calculate start index
	ixy = ( ( iy * nx_c ) + ix ) * ns_c;
	
	// begin transform
	idct = 0;
	for ( i = 0; i < ns_c; i++ )
		// idct += F[ i ] * icos_idct_fst[ ixy + i ];
		idct += F[ i ] * icos_idct_fst[ ixy++ ];
	
	
	return idct;	
}


/* -----------------------------------------------
	forward DCT transform using precalc tables (fast)
	----------------------------------------------- */
float fdct_2d_fst( unsigned char* f, int iu, int iv )
{
	float fdct;
	int iuv;
	int i;
	
	
	// calculate start index
	iuv = ( ( iv * nx_c ) + iu ) * ns_c;
	
	// begin transform
	fdct = 0;
	for ( i = 0; i < ns_c; i++ )
		// fdct += f[ i ] * icos_fdct_fst[ iuv + i ];
		fdct += f[ i ] * icos_fdct_fst[ iuv++ ];
	
	
	return fdct;	
}


/* -----------------------------------------------
	inverse DCT block transform
	----------------------------------------------- */

void  block_idct_2d( signed CTYPE* F, unsigned char* f )
{
	float idct;
	int ix, iy;
	int tmp;
	int i = 0;
	
	
	for ( iy = 0; iy < ny_c; iy++ ) {
		for ( ix = 0; ix < nx_c; ix++ ) {
			idct = (*idct_2d)( F, ix, iy );
			tmp = ROUND_F( idct );
			f[ i++ ] = CLAMPED( 0, 255, tmp );
			// f[ ix + ( nx_c * iy ) ] = CLAMPED( 0, 255, tmp );
		}
	}
}


/* -----------------------------------------------
	forward DCT block transform
	----------------------------------------------- */
void  block_fdct_2d( unsigned char* f, signed CTYPE* F )
{
	float fdct;
	int iu, iv;
	int tmp;
	int i = 0;
	
	
	for ( iv = 0; iv < ny_c; iv++ ) {
		for ( iu = 0; iu < nx_c; iu++ ) {
			fdct = (*fdct_2d)( f, iu, iv );
			tmp = ROUND_F( fdct );
			F[ i++ ] = tmp;
			// F[ iu + ( nx_c * iv ) ] = tmp;
		}
	}
}

/* ----------------------- End ofDCT specific functions -------------------------- */

/* ----------------------- Begin of miscellaneous helper functions -------------------------- */

/* -----------------------------------------------
	fetches first word from text, returns pos after it
	----------------------------------------------- */
char* fetch_command( char* text, char* cmd )
{
	int begin = -1;
	int end   = -1;
	int i;
	
	
	// return immediately if string is empty
	if ( text == NULL )
		return NULL;
	
	// seek begin of command
	for ( i = 0; text[ i ] != '\0'; i++ ) {
		if ( ( text[ i ] != ' ' ) && ( text[ i ] != '\n' ) && ( text[ i ] != '\r' ) ) {
			begin = i;
			break;
		}
	}
	
	// check begin
	if ( begin < 0 )
		return NULL;
	
	// seek end of command
	for ( i = begin; text[ i ] != '\0'; i++ ) {
		if ( ( text[ i ] == ' ' ) || ( text[ i ] == '\n' ) || ( text[ i ] == '\r' ) ) {
			break;
		}
	}
	
	// set end
	end = i;
	
	// copy string to command
	strncpy( cmd, text + begin, end - begin );
	cmd[ end - begin ] = '\0';
	
	// return last position
	return text + end;
}

/* -----------------------------------------------
	creates filename, callocs memory for it
	----------------------------------------------- */	
char* create_filename( const char* base, const char* extension )
{
	int len = strlen(base);
	int tol = 8;
	char* filename = (char*) calloc( len + tol, sizeof( char ) );
	
	set_extension( filename, base, extension );
	
	return filename;
}

/* -----------------------------------------------
	changes extension of filename
	----------------------------------------------- */	
void set_extension( char* destination, const char* origin, const char* extension )
{
	int i;
	
	int dotpos = 0;
	int length = strlen( origin );

	// find position of dot in filename
	for ( i = 0; i < length; i++ ) {
		if ( origin[i] == '.' ) {
			dotpos = i;
		}
	}
	
	if ( !dotpos ){
		dotpos = length;
	}
	
	strncpy( destination, origin, dotpos );
	destination[ dotpos ] = '.';
	strcpy( destination + dotpos + 1, extension );
}

/* -----------------------------------------------
	adds underscore after filename
	----------------------------------------------- */	
void add_underscore( char* filename )
{
	int i;
	
	int dotpos = 0;
	int length = strlen( filename );
	char* tmpname;
	
	// copy filename to tmpname
	tmpname = new char[ length + 1 ];
	strcpy( tmpname, filename );

	// find position of dot in filename
	for ( i = 0; i < length; i++ ) {
		if ( filename[ i ] == '.' ) {
			dotpos = i;
		}
	}
	
	// add underscore before extension
	if ( !dotpos ) {
		sprintf( filename, "%s_", tmpname );
	}
	else {
		strncpy( filename, tmpname, dotpos );
		sprintf( filename + dotpos, "_%s", tmpname + dotpos );
	}
}

/* ----------------------- End of miscellaneous helper functions -------------------------- */
