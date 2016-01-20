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

#include "bitops.h"
#include "tables.h"

#define QUANT(cmp,bpos) ( cmpnfo[cmp].qtable[ bpos ] )
#define MAX_V(cmp,bpos) ( ( freqmax[bpos] + QUANT(cmp,bpos) - 1 ) /  QUANT(cmp,bpos) )

#define ENVLI(s,v)		( ( v > 0 ) ? v : ( v - 1 ) + ( 1 << s ) )
#define DEVLI(s,n)		( ( n >= ( 1 << (s - 1) ) ) ? n : n + 1 - ( 1 << s ) )
#define E_ENVLI(s,v)	( v - ( 1 << s ) )
#define E_DEVLI(s,n)	( n + ( 1 << s ) )

#define ABS(v1)			( (v1 < 0) ? -v1 : v1 )
#define ABSDIFF(v1,v2)	( (v1 > v2) ? (v1 - v2) : (v2 - v1) )
#define IPOS(w,v,h)		( ( v * w ) + h )
#define NPOS(n1,n2,p)	( ( ( p / n1 ) * n2 ) + ( p % n1 ) )
#define ROUND_F(v1)		( (v1 < 0) ? (int) (v1 - 0.5) : (int) (v1 + 0.5) )
#define B_SHORT(v1,v2)	( ( ((int) v1) << 8 ) + ((int) v2) )
#define CLAMPED(l,h,v)	( ( v < l ) ? l : ( v > h ) ? h : v )

#define CSV_NAME		"JPEGinfo.csv"
#define EPS				(1e-8)
#define BARLEN			36

#define MEM_ERRMSG	"out of memory error"
#define FRD_ERRMSG	"can not read file: %s"
#define FWR_ERRMSG	"can not write file: %s"


/* -----------------------------------------------
	struct & enum declarations
	----------------------------------------------- */

enum F_TYPE {   JPEG = 0, UNK = 1		};

struct componentInfo {
	unsigned short* qtable; // quantization table
	int qtblno; // no of quantization table
	int huffdc; // no of huffman table (DC)
	int huffac; // no of huffman table (AC)
	int sfv; // sample factor vertical
	int sfh; // sample factor horizontal	
	int mbs; // blocks in mcu		
	int bcv; // block count vertical (interleaved)
	int bch; // block count horizontal (interleaved)
	int bc;  // block count (all) (interleaved)
	int ncv; // block count vertical (non interleaved)
	int nch; // block count horizontal (non interleaved)
	int nc;  // block count (all) (non interleaved)
	int sid; // statistical identity
	int jid; // jpeg internal id
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
bool check_stdin( void );
bool read_jpeg( void );
bool write_info( void );
bool reset_buffers( void );


/* -----------------------------------------------
	function declarations: jpeg-specific
	----------------------------------------------- */

bool setup_imginfo_jpg( void );
bool parse_jfif_jpg( unsigned char type, unsigned int len, unsigned char* segment );

/* -----------------------------------------------
	function declarations: JPG info
	----------------------------------------------- */

bool check_std_hufftables( void );
void calc_quality_range( int tno, float* low, float* high );
/* 09-03-08  Se */
/* 09-03-08  St */
int calcijgq( float scf );
bool detquant( int colc, int ijgqua );


/* -----------------------------------------------
	function declarations: miscelaneous helpers
	----------------------------------------------- */

char* create_filename( char* base, char* extension );
void set_extension( char* destination, char* origin, char* extension );
void add_underscore( char* filename );
void progress_bar( int current, int last );


/* -----------------------------------------------
	global variables: data storage
	----------------------------------------------- */

unsigned short qtables[4][64];		// quantization tables
unsigned char* hdrdata  =   NULL;   // header data
int            hdrs     =    0  ;   // size of header


/* -----------------------------------------------
	global variables: info about image
	----------------------------------------------- */

// seperate info for each color component
componentInfo cmpnfo[ 4 ];

int cmpc        = 0; // component count
int imgwidth    = 0; // width of image
int imgheight   = 0; // height of image

int sfhm        = 0; // max horizontal sample factor
int sfvm        = 0; // max verical sample factor
int mcuv        = 0; // mcus per line
int mcuh        = 0; // mcus per collumn
int mcuc        = 0; // count of mcus


/* -----------------------------------------------
	global variables: info about current scan
	----------------------------------------------- */

int cs_cmpc      =   0  ; // component count in current scan
int cs_cmp[ 4 ]  = { 0 }; // component numbers  in current scan
int cs_from      =   0  ; // begin - band of current scan ( inclusive )
int cs_to        =   0  ; // end - band of current scan ( inclusive )
int cs_sah       =   0  ; // successive approximation bit pos high
int cs_sal       =   0  ; // successive approximation bit pos low
int rsti         =   0  ;   // restart interval
	

/* -----------------------------------------------
	global variables: info about files
	----------------------------------------------- */
	
char*  jpgfilename;			// name of JPEG file
F_TYPE filetype;			// type of input file

int    jpgfilesize;			// size of JPEG file
int    jpegtype = 0;		// type of JPEG coding: 0->unknown, 1->sequential, 2->progressive

char** filelist = NULL;		// list of files to process 
int    file_cnt = 0;		// count of files in list
int    file_no  = 0;		// number of current file


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

int  verbosity  = -1;		// level of verbosity
bool overwrite  = false;	// overwrite files yes / no
int  err_tresh  = 2;		// error threshold ( proceed on warnings yes (2) / no (1) )

FILE*  msgout   = stdout;	// stream for output of messages
bool   pipe_on  = false;	// use stdin/stdout instead of filelist


/* -----------------------------------------------
	global variables: info about program
	----------------------------------------------- */

const unsigned char prgversion   = 10;
static const char*  subversion   = "b";
static const char*  versiondate  = "2009-04-24";
static const char*  author       = "Matthias Stirner/Gerhard Seelmann";


/* -----------------------------------------------
	main-function
	----------------------------------------------- */

int main( int argc, char** argv )
{	
	sprintf( statusmessage, "no statusmessage specified" );
	sprintf( errormessage, "no errormessage specified" );
	
	clock_t begin, end;
	
	int error_cnt   = 0;
	int warn_cnt    = 0;	
	int acc_jpgsize = 0;
	
	int speed, bpms;
	
	errorlevel = 0;
	
	
	// read options from command line
	initialize_options( argc, argv );
	
	// write program info to screen
	fprintf( msgout,  "\n--> jpeginfo v%i.%i%s (%s) by %s <--\n\n",
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
			if ( errorlevel < err_tresh ) {
				acc_jpgsize += jpgfilesize;
			}
		}
	}
	end = clock();
		
	// show statistics
	fprintf( msgout,  "\n\n-> %i file(s) processed, %i error(s), %i warning(s)\n",
		file_cnt, error_cnt, warn_cnt );
	if ( ( file_cnt > error_cnt ) && ( verbosity > 0 ) ) {
		speed = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC ); 
		bpms  = ( speed > 0 ) ? ( acc_jpgsize / speed ) : acc_jpgsize;
		
		fprintf( msgout,  " --------------------------------- \n" );
		fprintf( msgout,  " time taken        : %8i msec\n", speed );
		fprintf( msgout,  " avrg. byte per ms : %8i byte\n", bpms );
		fprintf( msgout,  " --------------------------------- \n" );
	}
	
	
	return file_cnt;
}

/* ----------------------- Begin of main interface functions -------------------------- */


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
		if ( sscanf( (*argv), "-v%i", &tmp_val ) == 1 ) {
			verbosity = tmp_val;
			verbosity = ( verbosity < -1 ) ? -1 : verbosity;
			verbosity = ( verbosity > 2 ) ? 2 : verbosity;			
		}
		else if ( strcmp((*argv), "-vp") == 0 ) {
			verbosity = -1;
		}
		else if ( strcmp((*argv), "-") == 0 ) {			
			msgout = stderr;
			// set binary mode for stdin & stdout
			#if !defined( unix )				
				setmode( fileno( stdin ), O_BINARY );
				setmode( fileno( stdout ), O_BINARY );
			#endif
			// use "-" as placeholder for stdin
			*(tmp_flp++) = "-";
		}
		else {
			// if argument is not switch, it's a filename
			*(tmp_flp++) = *argv;
		}		
	}
	
	// count number of files (or filenames) in filelist
	for ( file_cnt = 0; filelist[ file_cnt ] != NULL; file_cnt++ );	
}


/* -----------------------------------------------
	processes one file
	----------------------------------------------- */

void process_file( void )
{
	clock_t begin, end;
	
	
	errorfunction = NULL;
	errorlevel = 0;
	jpgfilesize = 0;
	
	
	// compare file name, set pipe if needed
	if ( strcmp( filelist[ file_no ], "-" ) == 0 ) {
		pipe_on = true;
		filelist[ file_no ] = "STDIN";
	}
	
	if ( verbosity == -1 ) { // progress bar UI
		fprintf( msgout, "Processing file %2i of %2i ", file_no + 1, file_cnt );
		progress_bar( file_no, file_cnt );
		fprintf( msgout, "\r" );
		
		// check input file and determine filetype
		execute( ( pipe_on ) ? check_stdin : check_file );	
		
		// main function routine
		begin = clock();
		
		if ( filetype == JPEG )	{
			execute( read_jpeg );
			execute( write_info );
		}
		// reset buffers
		reset_buffers();
		
		end = clock();	
		
		// if this is the last file, update progress bar one last time
		if ( file_no + 1 == file_cnt ) {
			// update progress message
			fprintf( msgout, "Processed %2i of %2i files ", file_no + 1, file_cnt );
			progress_bar( 1, 1 );
			fprintf( msgout, "\r" );
		}
	}
	else { // 'old' UI
		char* actionmsg  = NULL;
		char* errtypemsg = NULL;
		int speed, bpms;
	
		fprintf( msgout,  "\nProcessing file %i of %i \"%s\" -> ",
			file_no + 1, file_cnt, filelist[ file_no ] );
		if ( verbosity > 1 )
			fprintf( msgout,  "\n----------------------------------------" );
		
		// check input file and determine filetype
		execute( ( pipe_on ) ? check_stdin : check_file );
		
		// get specific action message
		if ( filetype == UNK ) actionmsg = "unknown filetype";
		else actionmsg = "Writing info";
		
		if ( verbosity < 2 ) fprintf( msgout, "%s -> ", actionmsg );
		
		
		// main function routine
		begin = clock();
		
		if ( filetype == JPEG )	{
			execute( read_jpeg );
			execute( write_info );
		}
		// reset buffers
		reset_buffers();
		
		end = clock();
		
		
		// speed and compression ratio calculation
		speed = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );
		bpms  = ( speed > 0 ) ? ( jpgfilesize / speed ) : jpgfilesize;

		
		switch ( verbosity )
		{
			case 0:
				if ( errorlevel < err_tresh )
					fprintf( msgout,  "OK" );
				else fprintf( msgout,  "ERROR" );
				if ( errorlevel > 0 )
					fprintf( msgout,  "\n" );
				break;
			
			case 1:
				if ( errorlevel < err_tresh ) fprintf( msgout,  "DONE\n" );
				else fprintf( msgout,  "ERROR\n" );
				break;
			
			case 2:
				fprintf( msgout,  "\n----------------------------------------\n" );
				if ( errorlevel < err_tresh ) fprintf( msgout,  "-> %s OK\n", actionmsg );
				break;
		}
		
		switch ( errorlevel )
		{
			case 0:
				errtypemsg = "none";
				break;
				
			case 1:
				if ( errorlevel < err_tresh )
					errtypemsg = "warning (ignored)";
				else
					errtypemsg = "warning (skipped file)";
				break;
			
			case 2:
				errtypemsg = "fatal error";
				break;
		}
		
		if ( errorlevel > 0 )
		{
			get_status( errorfunction );
			fprintf( stderr, " %s -> %s:\n", statusmessage, errtypemsg  );
			fprintf( stderr, " %s\n", errormessage );
			if ( verbosity > 1 )
				fprintf( stderr, " (in file \"%s\")\n", filelist[ file_no ] );
		}
		if ( (verbosity > 0) && (errorlevel < err_tresh) )
		{		
			fprintf( msgout,  " time taken  : %7i msec\n", speed );
			fprintf( msgout,  " byte per ms : %7i byte\n", bpms );
		}
		
		if (verbosity > 1)
			fprintf( msgout,  "\n" );
	}
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
		if ( verbosity == 2 ) {
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
			if ( verbosity == 2 ) fprintf( msgout,  "%6ims",
				(int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC ) );
		}
		else {
			errorfunction = function;
			if ( verbosity == 2 ) fprintf( msgout,  "%8s", "ERROR" );
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
		sprintf( statusmessage, "Determining filetype" );
	}
	else if ( function == *check_stdin ) {
		sprintf( statusmessage, "Determining filetype" );
	}
	else if ( function == *read_jpeg ) {
		sprintf( statusmessage, "Reading header & image data" );
	}
	else if ( function == *write_info ) {
		sprintf( statusmessage, "Writing info to csv" );
	}
	else if ( function == *reset_buffers ) {
		sprintf( statusmessage, "Resetting program" );
	}
}


/* -----------------------------------------------
	shows help in case of wrong input
	----------------------------------------------- */
	
void show_help( void )
{
	fprintf( msgout, "\n" );
	fprintf( msgout, "Usage: jpeginfo [switches] [filename(s)]");
	fprintf( msgout, "\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-v?]    set level of verbosity (max: 2) (def: 0)\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, "Example: \"jpeginfo -v1 baboon.jpg\"\n" );
}

/* ----------------------- End of main interface functions -------------------------- */

/* ----------------------- Begin of main functions -------------------------- */


/* -----------------------------------------------
	check file and determine filetype
	----------------------------------------------- */
	
bool check_file( void )
{
	FILE* fp;
	unsigned char fileid[ 2 ] = { 0, 0 };
	int namelen = strlen( filelist[ file_no ] ) + 5;
	
	
	// free memory from filenames if needed
	if ( jpgfilename != NULL ) free( jpgfilename );
	// alloc memory for filenames
	jpgfilename = (char*) calloc( namelen, sizeof( char ) );
	
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
	
	// immediately return error if 2 bytes can't be read
	if ( fread( fileid, 1, 2, fp ) != 2 ) {
		filetype = UNK;
		fclose( fp );
		sprintf( errormessage, "file doesn't contain enough data" );
		errorlevel = 2;
		return false;
	}
	
	// check file id, determine filetype
	if ( (fileid[0] == 0xFF) && (fileid[1] == 0xD8) ) {
		// file is JPEG
		filetype = JPEG;
		// copy filename
		strcpy( jpgfilename, filelist[ file_no ] );
	}
	else {
		// file is neither
		filetype = UNK;
		fclose( fp );
		sprintf( errormessage, "filetype of file \"%s\" is unknown", filelist[ file_no ] );
		errorlevel = 2;
		return false;		
	}
	
	// close file
	fclose ( fp );
	
	
	return true;
}


/* -----------------------------------------------
	check stdin and determine filetype
	----------------------------------------------- */
	
bool check_stdin( void )
{	
	unsigned char fileid[ 2 ] = { 0, 0 };
	
	
	// checken ob file an stdin ( ob puffer gefüllt? ) -> TODO!
	if ( fread( fileid, 1, 2, stdin ) != 2 ) {
		filetype = UNK;
		sprintf( errormessage, "file doesn't contain enough data" );
		errorlevel = 2;
		return false;
	}
	
	// check file id, determine filetype
	if ( ( fileid[0] == 0xFF) && (fileid[1] == 0xD8) ) {
		// file is JPEG
		filetype = JPEG;
		// create filenames
		jpgfilename = "JPG file from stdin";
	}
	else {
		// file is neither
		filetype = UNK;
		sprintf( errormessage, "filetype is unknown" );
		errorlevel = 2;
		return false;		
	}
	
	
	return true;
}


/* -----------------------------------------------
	Read in header & image data
	----------------------------------------------- */
	
bool read_jpeg( void )
{
	unsigned char* segment = NULL; // storage for current segment
	unsigned int   ssize = 1024; // current size of segment array
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned char  tmp;	
	
	abytewriter* hdrw;
	FILE *fp = NULL;
	
	
	if ( !pipe_on ) {
		// open file for input	
		fp = fopen( jpgfilename, "rb" );
		if ( fp == NULL ) {
			sprintf( errormessage, FRD_ERRMSG, jpgfilename );
			errorlevel = 2;
			return false;
		}
		// skip SOI
		fseek( fp, 2, SEEK_SET );
	}
	else {
		// use stdin
		fp = stdin;
	}
	
	// start headerwriter
	hdrw = new abytewriter( 4096 );
	hdrs = 0; // size of header data, start with 0
	
	// alloc memory for segment data
	segment = ( unsigned char* ) calloc( ssize, sizeof( char ) );
	if ( segment == NULL ) {
		if ( !pipe_on ) fclose( fp );
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// JPEG reader loop
	while ( true ) {		
		if ( type == 0xDA ) { // if last marker was sos
			// switch to huffman data reading mode
			while ( true ) {
				// read byte from imagedata
				if ( fread( &tmp, 1, 1, fp ) == 0 )
					break; // break of EOF encountered
				// check if 0xFF encountered
				if ( tmp == 0xFF ) {					
					fread( &tmp, 1, 1, fp ); // read next byte & check
					if ( ( tmp != 0x00 ) && ( ( tmp < 0xD0 ) || ( tmp > 0xD7 ) ) ) {
						// switch back to header parsing for all markers other than FILL (0x00) or RST (0xD0...0xD7) 
						segment[ 0 ] = 0xFF;
						segment[ 1 ] = tmp;
						break;
					}
				}
			}
		}
		else {
			// read in next marker
			if ( fread( segment, 1, 2, fp ) != 2 ) break;
			if ( segment[ 0 ] != 0xFF ) {
				// ugly fix for incorrect marker segment sizes
				sprintf( errormessage, "size mismatch in marker segment FF %2X", type );
				errorlevel = 2;
				if ( type == 0xFE ) { //  if last marker was COM try again
					if ( fread( segment, 1, 2, fp ) != 2 ) break;
					if ( segment[ 0 ] == 0xFF ) errorlevel = 1;
				}
				if ( errorlevel == 2 ) {
					if ( !pipe_on ) fclose( fp );
					delete ( hdrw );
					free ( segment );
					return false;
				}
			}
		}
		
		// read segment type
		type = segment[ 1 ];
		
		// if EOI is encountered make a quick exit
		if ( type == 0xD9 ) {
			// get pointer for header data & size
			hdrdata  = hdrw->getptr();
			hdrs     = hdrw->getpos();
			// everything is done here now
			break;			
		}
		
		// read in next segments' length and check it
		if ( fread( segment + 2, 1, 2, fp ) != 2 ) break;
		len = 2 + B_SHORT( segment[ 2 ], segment[ 3 ] );
		
		// realloc segment data if needed
		if ( ssize < len ) {
			segment = ( unsigned char* ) realloc( segment, len );
			if ( segment == NULL ) {
				if ( !pipe_on ) fclose( fp );
				sprintf( errormessage, MEM_ERRMSG );
				errorlevel = 2;
				delete ( hdrw );
				return false;
			}
			ssize = len;
		}
		
		// read rest of segment, store back in header writer
		if ( fread( ( segment + 4 ), 1, ( len - 4 ), fp ) !=
			( unsigned short ) ( len - 4 ) ) break;
		hdrw->write_n( segment, len );
	}
	// JPEG reader loop end
	
	// free writers
	delete ( hdrw );
	
	// check if everything went OK
	if ( hdrs == 0 ) {
		if ( !pipe_on ) fclose( fp );
		sprintf( errormessage, "unexpected end of data encountered" );
		errorlevel = 2;
		return false;
	}
	
	// free segment
	free( segment );
	
	// get filesize
	fseek( fp, 0, SEEK_END );
	jpgfilesize = ftell( fp );
	
	// close file
	if ( !pipe_on )
		fclose( fp );
	
	// parse header for image info
	if ( !setup_imginfo_jpg() ) {
		return false;
	}
	
	
	return true;
}


/* -----------------------------------------------
	write info to csv file
	----------------------------------------------- */

bool write_info( void )
{
	float ijg_quality_l[ 4 ] = { 0.0 };
	float ijg_quality_h[ 4 ] = { 0.0 };
	float comp_ratio = 0;
	bool new_file = true;
	bool std_htables;
	
	FILE* fp;	
	int cmp;
	int i;
	int ijgqu0, ijgqu1;
	char cssstr[ 8 ]; // string for color subsampling
	
	
	fp = fopen( CSV_NAME, "r" );
	if ( fp != NULL ) {
		new_file = false;
		fclose( fp );
	}
	fp = fopen( CSV_NAME, "a" );
	if ( fp == NULL ) {
		sprintf( errormessage, FWR_ERRMSG, CSV_NAME );
		errorlevel = 2;
		return false;
	}
	
	// gather additional info
	std_htables = check_std_hufftables();
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		for ( i = 0; i < 4; i++ )
			if ( cmpnfo[ cmp ].qtable == qtables[ i ] ) break;
		if ( i < 4 ) calc_quality_range( i, &(ijg_quality_l[cmp]), &(ijg_quality_h[cmp]) );
	}
	
	// Total/Lum Pixel ratio
	cmpnfo[1].sfh = (cmpc > 1) ? (cmpnfo[1].sfh) : (0); 
	cmpnfo[1].sfv = (cmpc > 1) ? (cmpnfo[1].sfv) : (0); 
	comp_ratio = 1.0*(cmpnfo[0].sfh*cmpnfo[0].sfv + cmpnfo[1].sfh*cmpnfo[1].sfv) / (cmpnfo[0].sfh*cmpnfo[0].sfv); 
	comp_ratio = comp_ratio * imgwidth * imgheight /(jpgfilesize - hdrs);
	
	// subsampling string
	strcpy( cssstr, "unk" );
	if ( ( cmpnfo[1].sfh == 1 ) && ( cmpnfo[1].sfv == 1 ) ) {
		if ( ( cmpnfo[0].sfh == 1 ) && ( cmpnfo[0].sfv == 1 ) ) strcpy( cssstr, "4:4:4" );
		if ( ( cmpnfo[0].sfh == 2 ) && ( cmpnfo[0].sfv == 1 ) ) strcpy( cssstr, "4:4:0" );
		if ( ( cmpnfo[0].sfh == 1 ) && ( cmpnfo[0].sfv == 2 ) ) strcpy( cssstr, "4:2:2" );
		if ( ( cmpnfo[0].sfh == 2 ) && ( cmpnfo[0].sfv == 2 ) ) strcpy( cssstr, "4:2:0" );
	}
	if ( cmpc == 1 ) strcpy( cssstr, "4:0:0" );	
	
	// heuristc determination ijg quality via quatization scaling factor	
	ijgqu0 = ( cmpc > 0 ) ? ( ( cmpnfo[0].qtable[58] == 1 ) ? 100 : calcijgq( cmpnfo[0].qtable[58] / 100.0 ) ) : 0;
	ijgqu1 = ( cmpc > 1 ) ? ( ( cmpnfo[1].qtable[58] == 1 ) ? 100 : calcijgq( cmpnfo[1].qtable[58] / 99.0 ) ) : 0;
	// further treatment if quantizer step size is at end of range
	if ( cmpnfo[0].qtable[58] >= 255 ) ijgqu0 = calcijgq( cmpnfo[0].qtable[5] / 10.0 );
	if ( cmpnfo[0].qtable[ 5] >= 255 ) ijgqu0 = 1;
	if ( cmpc > 1 ) {
		if ( cmpnfo[1].qtable[58] >= 255 ) ijgqu1 = calcijgq( cmpnfo[1].qtable[0] / 17.0 );
		if ( ( cmpnfo[1].qtable[ 0] >= 255 ) && ( cmpnfo[0].qtable[58] >= 255 ) ) ijgqu1 = ijgqu0;
	}
	
	if ( new_file ) { // write collumn titles only if file is empty
		fprintf( fp, "File name; " );
		fprintf( fp, "Width; " );
		fprintf( fp, "Height; " );
		fprintf( fp, "File size; " );
		fprintf( fp, "Header size; " );
		fprintf( fp, "Col. Subs.; " );
		fprintf( fp, "IJG Q.(0); " );
		fprintf( fp, "IJG Q.(1); " );
		fprintf( fp, "Quant(0); " );
		fprintf( fp, "Quant(1); " );
		fprintf( fp, "Img. compr.; " );
		fprintf( fp, "Huff tables; " );
		fprintf( fp, "Jpeg mode; " );
		
		if (verbosity > 0) {
			fprintf( fp, "No. of comp.; " );
			fprintf( fp, "samples (0); " );
			fprintf( fp, "samples (1); " );
			fprintf( fp, "samples (2); " );
			fprintf( fp, "samples (3); " );
			fprintf( fp, "quality (0); " );
			fprintf( fp, "quality (1); " );
			fprintf( fp, "quality (2); " );
			fprintf( fp, "quality (3); " );
			fprintf( fp, "quantDC (0); " );
			fprintf( fp, "quantDC (1); " );
			fprintf( fp, "quantDC (2); " );
			fprintf( fp, "quantDC (3); " );
		}
		fprintf( fp, "\n" );
	}
	
	// write data to file
	fprintf( fp, "%s; ", jpgfilename );
	fprintf( fp, "%i; ", imgwidth );
	fprintf( fp, "%i; ", imgheight );
	fprintf( fp, "%i; ", jpgfilesize );
	fprintf( fp, "%i; ", hdrs );
	fprintf( fp, "%s; ", cssstr );		
	fprintf( fp, "%i; ", ijgqu0 );
	fprintf( fp, "%i; ", ijgqu1 );	
	fprintf( fp, "%s; ", ( cmpc > 0 ) ? ( ( detquant(0, ijgqu0) ) ? "def" : "cst" ) : ("---") );
	fprintf( fp, "%s; ", ( cmpc > 1 ) ? ( ( detquant(1, ijgqu1) ) ? "def" : "cst" ) : ("---") );
	fprintf( fp, "%5.2f; ", comp_ratio );
	fprintf( fp, "%s; ", ( std_htables ) ? "std" : "opt" );
	fprintf( fp, "%s; ", ( jpegtype == 0 ) ? "unk" : ( jpegtype == 1 ) ? "seq" : "prg");

	if ( verbosity > 0 ) {
		fprintf( fp, "%i; ", cmpc );
		fprintf( fp, "%ix%i; ", ( cmpc > 0 ) ? cmpnfo[0].sfh : 0, ( cmpc > 0 ) ? cmpnfo[0].sfv : 0 );
		fprintf( fp, "%ix%i; ", ( cmpc > 1 ) ? cmpnfo[1].sfh : 0, ( cmpc > 1 ) ? cmpnfo[1].sfv : 0 );
		fprintf( fp, "%ix%i; ", ( cmpc > 2 ) ? cmpnfo[2].sfh : 0, ( cmpc > 2 ) ? cmpnfo[2].sfv : 0 );
		fprintf( fp, "%ix%i; ", ( cmpc > 3 ) ? cmpnfo[3].sfh : 0, ( cmpc > 3 ) ? cmpnfo[3].sfv : 0 );
		fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 0 ) ? ijg_quality_l[0] * 100.0 : 0, ( cmpc > 0 ) ? ijg_quality_h[0] * 100.0 : 0 );
		fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 1 ) ? ijg_quality_l[1] * 100.0 : 0, ( cmpc > 1 ) ? ijg_quality_h[1] * 100.0 : 0 );
		fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 2 ) ? ijg_quality_l[2] * 100.0 : 0, ( cmpc > 2 ) ? ijg_quality_h[2] * 100.0 : 0 );
		fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 3 ) ? ijg_quality_l[3] * 100.0 : 0, ( cmpc > 3 ) ? ijg_quality_h[3] * 100.0 : 0 );
		fprintf( fp, "%i; ", ( cmpc > 0 ) ? cmpnfo[0].qtable[0] : 0 );
		fprintf( fp, "%i; ", ( cmpc > 1 ) ? cmpnfo[1].qtable[0] : 0 );	
		fprintf( fp, "%i; ", ( cmpc > 2 ) ? cmpnfo[2].qtable[0] : 0 );
		fprintf( fp, "%i; ", ( cmpc > 3 ) ? cmpnfo[3].qtable[0] : 0 );
	}
	fprintf( fp, "\n" );
	
	fclose( fp );
	
	
	return true;
}


/* -----------------------------------------------
	set each variable to it initial value
	----------------------------------------------- */

bool reset_buffers( void )
{
	int cmp, bpos;
	int i;
	
	
	// -- free buffers --
	
	// free buffers & set pointers NULL
	if ( hdrdata  != NULL ) free ( hdrdata );
	hdrdata   = NULL;
	
	
	// -- set variables --
	
	// preset componentinfo
	for ( cmp = 0; cmp < 4; cmp++ ) {
		cmpnfo[ cmp ].sfv = -1;
		cmpnfo[ cmp ].sfh = -1;
		cmpnfo[ cmp ].mbs = -1;
		cmpnfo[ cmp ].bcv = -1;
		cmpnfo[ cmp ].bch = -1;
		cmpnfo[ cmp ].bc  = -1;
		cmpnfo[ cmp ].ncv = -1;
		cmpnfo[ cmp ].nch = -1;
		cmpnfo[ cmp ].nc  = -1;
		cmpnfo[ cmp ].sid = -1;
		cmpnfo[ cmp ].jid = -1;
		cmpnfo[ cmp ].qtable = NULL;
		cmpnfo[ cmp ].qtblno = -1;
		cmpnfo[ cmp ].huffdc = -1;
		cmpnfo[ cmp ].huffac = -1;
	}
	
	// preset imgwidth / imgheight / component count 
	imgwidth  = 0;
	imgheight = 0;
	cmpc      = 0;
	
	// preset mcu info variables / restart interval
	sfhm      = 0;
	sfvm      = 0;
	mcuc      = 0;
	mcuh      = 0;
	mcuv      = 0;
	rsti      = 0;
	
	// reset quantization  tables
	for ( i = 0; i < 4; i++ ) {
		for ( bpos = 0; bpos < 64; bpos++ )
			qtables[ i ][ bpos ] = 0;
	}
	
	// preset jpegtype
	jpegtype  = 0;
	
	
	return true;
}

/* ----------------------- End of main functions -------------------------- */

/* ----------------------- Begin of JPEG specific functions -------------------------- */


/* -----------------------------------------------
	Parses header for imageinfo
	----------------------------------------------- */
bool setup_imginfo_jpg( void )
{
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header
	
	int cmp;
	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		// do not parse DHT & DRI
		if ( ( type != 0xDA ) && ( type != 0xC4 ) && ( type != 0xDD ) ) {
			if ( !parse_jfif_jpg( type, len, &( hdrdata[ hpos ] ) ) )
				return false;
		}
		hpos += len;
	}
	
	// check if information is complete
	if ( cmpc == 0 ) {
		sprintf( errormessage, "header contains incomplete information" );
		errorlevel = 2;
		return false;
	}
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		if ( ( cmpnfo[cmp].sfv == 0 ) ||
			 ( cmpnfo[cmp].sfh == 0 ) ||
			 ( cmpnfo[cmp].qtable == NULL ) ||
			 ( cmpnfo[cmp].qtable[0] == 0 ) ||
			 ( jpegtype == 0 ) ) {
			sprintf( errormessage, "header information is incomplete" );
			errorlevel = 2;
			return false;
		}
	}
	
	// do all remaining component info calculations
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		if ( cmpnfo[ cmp ].sfh > sfhm ) sfhm = cmpnfo[ cmp ].sfh;
		if ( cmpnfo[ cmp ].sfv > sfvm ) sfvm = cmpnfo[ cmp ].sfv;
	}
	mcuv = ( int ) ceil( (float) imgheight / (float) ( 8 * sfhm ) );
	mcuh = ( int ) ceil( (float) imgwidth  / (float) ( 8 * sfvm ) );
	mcuc  = mcuv * mcuh;
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		cmpnfo[ cmp ].mbs = cmpnfo[ cmp ].sfv * cmpnfo[ cmp ].sfh;		
		cmpnfo[ cmp ].bcv = mcuv * cmpnfo[ cmp ].sfh;
		cmpnfo[ cmp ].bch = mcuh * cmpnfo[ cmp ].sfv;
		cmpnfo[ cmp ].bc  = cmpnfo[ cmp ].bcv * cmpnfo[ cmp ].bch;
		cmpnfo[ cmp ].ncv = ( int ) ceil( (float) imgheight * 
							( (float) cmpnfo[ cmp ].sfh / ( 8.0 * sfhm ) ) );
		cmpnfo[ cmp ].nch = ( int ) ceil( (float) imgwidth * 
							( (float) cmpnfo[ cmp ].sfv / ( 8.0 * sfvm ) ) );
		cmpnfo[ cmp ].nc  = cmpnfo[ cmp ].ncv * cmpnfo[ cmp ].nch;
	}
	
	// decide components' statistical ids
	if ( cmpc <= 3 ) {
		for ( cmp = 0; cmp < cmpc; cmp++ ) cmpnfo[ cmp ].sid = cmp;
	}
	else {
		for ( cmp = 0; cmp < cmpc; cmp++ ) cmpnfo[ cmp ].sid = 0;
	}
	
	
	return true;
}


/* -----------------------------------------------
	Parse routines for JFIF segments
	----------------------------------------------- */
bool parse_jfif_jpg( unsigned char type, unsigned int len, unsigned char* segment )
{
	unsigned int hpos = 4; // current position in segment, start after segment header
	int lval, rval; // temporary variables
	int cmp;
	int i;
	
	
	switch ( type )
	{
		case 0xC4: // DHT segment
			// build huffman trees & codes
			// nothing to do here!
			return true;
		
		case 0xDB: // DQT segment
			// copy quantization tables to internal memory
			while ( hpos < len ) {
				lval = LBITS( segment[ hpos ], 4 );
				rval = RBITS( segment[ hpos ], 4 );
				if ( (lval < 0) || (lval >= 2) ) break;
				if ( (rval < 0) || (rval >= 4) ) break;
				hpos++;				
				if ( lval == 0 ) { // 8 bit precision
					for ( i = 0; i < 64; i++ ) {
						qtables[ rval ][ i ] = ( unsigned short ) segment[ hpos + i ];
						if ( qtables[ rval ][ i ] == 0 ) break;
					}
					hpos += 64;
				}
				else { // 16 bit precision
					for ( i = 0; i < 64; i++ ) {
						qtables[ rval ][ i ] =
							B_SHORT( segment[ hpos + (2*i) ], segment[ hpos + (2*i) + 1 ] );
						if ( qtables[ rval ][ i ] == 0 ) break;
					}
					hpos += 128;
				}
			}
			
			if ( hpos != len ) {
				// if we get here, something went wrong
				sprintf( errormessage, "size mismatch in dqt marker" );
				errorlevel = 2;
				return false;
			}
			return true;
			
		case 0xDD: // DRI segment
			// define restart interval
			rsti = B_SHORT( segment[ hpos ], segment[ hpos + 1 ] );			
			return true;
			
		case 0xDA: // SOS segment
			// prepare next scan
			cs_cmpc = segment[ hpos ];
			if ( cs_cmpc > cmpc ) {
				sprintf( errormessage, "%i components in scan, only %i are allowed",
							cs_cmpc, cmpc );
				errorlevel = 2;
				return false;
			}
			hpos++;
			for ( i = 0; i < cs_cmpc; i++ ) {
				for ( cmp = 0; ( segment[ hpos ] != cmpnfo[ cmp ].jid ) && ( cmp < cmpc ); cmp++ );
				if ( cmp == cmpc ) {
					sprintf( errormessage, "component id mismatch in start-of-scan" );
					errorlevel = 2;
					return false;
				}
				cs_cmp[ i ] = cmp;
				cmpnfo[ cmp ].huffdc = LBITS( segment[ hpos + 1 ], 4 );
				cmpnfo[ cmp ].huffac = RBITS( segment[ hpos + 1 ], 4 );
				if ( ( cmpnfo[ cmp ].huffdc < 0 ) || ( cmpnfo[ cmp ].huffdc >= 4 ) ||
					 ( cmpnfo[ cmp ].huffac < 0 ) || ( cmpnfo[ cmp ].huffac >= 4 ) ) {
					sprintf( errormessage, "huffman table number mismatch" );
					errorlevel = 2;
					return false;
				}
				hpos += 2;
			}
			cs_from = segment[ hpos + 0 ];
			cs_to   = segment[ hpos + 1 ];
			cs_sah  = LBITS( segment[ hpos + 2 ], 4 );
			cs_sal  = RBITS( segment[ hpos + 2 ], 4 );
			// check for errors
			if ( ( cs_from > cs_to ) || ( cs_from > 63 ) || ( cs_to > 63 ) ) {
				sprintf( errormessage, "spectral selection parameter out of range" );
				errorlevel = 2;
				return false;
			}
			if ( ( cs_sah >= 12 ) || ( cs_sal >= 12 ) ) {
				sprintf( errormessage, "successive approximation parameter out of range" );
				errorlevel = 2;
				return false;
			}
			return true;
		
		case 0xC0: // SOF0 segment
			// coding process: baseline DCT
			
		case 0xC1: // SOF1 segment
			// coding process: extended sequential DCT
		
		case 0xC2: // SOF2 segment
			// coding process: progressive DCT
			
			// set JPEG coding type
			if ( type == 0xC2 )
				jpegtype = 2;
			else
				jpegtype = 1;
				
			// check data precision, only 8 bit is allowed
			lval = segment[ hpos ];
			if ( lval != 8 ) {
				sprintf( errormessage, "%i bit data precision is not supported", lval );
				errorlevel = 2;
				return false;
			}
			
			// image size, height & component count
			imgheight = B_SHORT( segment[ hpos + 1 ], segment[ hpos + 2 ] );
			imgwidth  = B_SHORT( segment[ hpos + 3 ], segment[ hpos + 4 ] );
			cmpc      = segment[ hpos + 5 ];
			if ( cmpc > 4 ) {
				sprintf( errormessage, "image has %i components, max 4 are supported", cmpc );
				errorlevel = 2;
				return false;
			}
			
			hpos += 6;
			// components contained in image
			for ( cmp = 0; cmp < cmpc; cmp++ ) {
				cmpnfo[ cmp ].jid = segment[ hpos ];
				cmpnfo[ cmp ].sfv = LBITS( segment[ hpos + 1 ], 4 );
				cmpnfo[ cmp ].sfh = RBITS( segment[ hpos + 1 ], 4 );				
				cmpnfo[ cmp ].qtable = qtables[ segment[ hpos + 2 ] ];
				cmpnfo[ cmp ].qtblno = segment[ hpos + 2 ];
				hpos += 3;
			}
			
			return true;
		
		case 0xC3: // SOF3 segment
			// coding process: lossless sequential
			sprintf( errormessage, "sof3 marker found, image is coded lossless" );
			errorlevel = 2;
			return false;
		
		case 0xC5: // SOF5 segment
			// coding process: differential sequential DCT
			sprintf( errormessage, "sof5 marker found, image is coded diff. sequential" );
			errorlevel = 2;
			return false;
		
		case 0xC6: // SOF6 segment
			// coding process: differential progressive DCT
			sprintf( errormessage, "sof6 marker found, image is coded diff. progressive" );
			errorlevel = 2;
			return false;
		
		case 0xC7: // SOF7 segment
			// coding process: differential lossless
			sprintf( errormessage, "sof7 marker found, image is coded diff. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xC9: // SOF9 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof9 marker found, image is coded arithm. sequential" );
			errorlevel = 2;
			return false;
			
		case 0xCA: // SOF10 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof10 marker found, image is coded arithm. progressive" );
			errorlevel = 2;
			return false;
			
		case 0xCB: // SOF11 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof11 marker found, image is coded arithm. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xCD: // SOF13 segment
			// coding process: arithmetic differntial sequential DCT
			sprintf( errormessage, "sof13 marker found, image is coded arithm. diff. sequential" );
			errorlevel = 2;
			return false;
			
		case 0xCE: // SOF14 segment
			// coding process: arithmetic differential progressive DCT
			sprintf( errormessage, "sof14 marker found, image is coded arithm. diff. progressive" );
			errorlevel = 2;
			return false;
		
		case 0xCF: // SOF15 segment
			// coding process: arithmetic differntial lossless
			sprintf( errormessage, "sof15 marker found, image is coded arithm. diff. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xE0: // APP0 segment	
		case 0xE1: // APP1 segment
		case 0xE2: // APP2 segment
		case 0xE3: // APP3 segment
		case 0xE4: // APP4 segment
		case 0xE5: // APP5 segment
		case 0xE6: // APP6 segment
		case 0xE7: // APP7 segment
		case 0xE8: // APP8 segment
		case 0xE9: // APP9 segment
		case 0xEA: // APP10 segment
		case 0xEB: // APP11 segment
		case 0xEC: // APP12segment
		case 0xED: // APP13 segment
		case 0xEE: // APP14 segment
		case 0xEF: // APP15 segment
		case 0xFE: // COM segment
			// do nothing - return true
			return true;
			
		case 0xD0: // RST0 segment
		case 0xD1: // RST1segment
		case 0xD2: // RST2 segment
		case 0xD3: // RST3 segment
		case 0xD4: // RST4 segment
		case 0xD5: // RST5 segment
		case 0xD6: // RST6 segment
		case 0xD7: // RST7 segment
			// return errormessage - RST is out of place here
			sprintf( errormessage, "rst marker found out of place" );
			errorlevel = 2;
			return false;
		
		case 0xD8: // SOI segment
			// return errormessage - start-of-image is out of place here
			sprintf( errormessage, "soi marker found out of place" );
			errorlevel = 2;
			return false;
		
		case 0xD9: // EOI segment
			// return errormessage - end-of-image is out of place here
			sprintf( errormessage, "eoi marker found out of place" );
			errorlevel = 2;
			return false;
			
		default: // unknown marker segment
			// return warning
			sprintf( errormessage, "unknown marker found: FF %2X", type );
			errorlevel = 1;
			return true;
	}
}

/* ----------------------- End of JPEG specific functions -------------------------- */

/* ----------------------- Begin of JPEG information functions -------------------------- */

/* -----------------------------------------------
	check if standard huffman tables are used
	----------------------------------------------- */
bool check_std_hufftables( void )
{
	unsigned char type = 0x00; // type of current marker segment
	int len  = 0; // length of current marker segment
	int hpos = 0; // position in header
	int fpos; // end of marker position
	int skip; // bytes to skip
	int spos; // sub position
	int i;
	
	
	// doesn't apply to coding modes other than sequential
	if ( jpegtype != 1 ) return false;
	
	// search for DHT (0xFFC4) marker segments	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		if ( type == 0xC4 ) { // for DHT
			fpos = hpos + len; // reassign length to end position
			hpos += 4; // skip marker & length
			while ( hpos < fpos ) {			
				hpos++;				
				// table found - compare with each of the four standard tables		
				for ( i = 0; i < 4; i++ ) {
					for ( spos = 0; spos < std_huff_lengths[ i ]; spos++ ) {
						if ( hdrdata[ hpos + spos ] != std_huff_tables[ i ][ spos ] )
							break;
					}
					// check if comparison ok
					if ( spos == std_huff_lengths[ i ] )
						break;
				}
				if ( i == 4 ) return false;
				skip = 16;
				for ( i = 0; i < 16; i++ )		
					skip += ( int ) hdrdata[ hpos + i ];				
				hpos += skip;
			}
		}
		else { // skip segment
			hpos += len;
		}		
	}	
	
	
	return true;
}


/* -----------------------------------------------
	calculate current quality range for table #
	----------------------------------------------- */
void calc_quality_range( int tno, float* low, float* high )
{
	float s = 0;
	float sl = 0;
	float sh = 65536;
	int ref = -1; // # of reference table
	int cmp;
	int i;
	
	
	// decide reference table
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		if ( cmpnfo[ cmp ].qtable == qtables[ tno ] ) {
			ref = cmpnfo[ cmp ].sid;
			ref = CLAMPED( 0, 1, ref );
			break;
		}
	}
	
	// check if no reference was found
	if ( ref == -1 ) {
		(*low) = 0;
		(*high) = 0;
		return;
	}
	
	// compare tables, calculate average scaling factor s
	for ( i = 0; i < 64; i++ ) {
		s = ( ( float ) qtables[ tno ][ i ] ) / ( ( float ) std_qtables[ ref ][ unzigzag[i] ] );
		sl = ( s > sl ) ? s : sl;
		sh = ( s < sh ) ? s : sh;
	}
	
	// calculate qualities based on s
	(*low)  = ( sl < 1.0 ) ? 1.0 - (sl/2) : 1.0 / (2*sl);
	(*high) = ( sh < 1.0 ) ? 1.0 - (sh/2) : 1.0 / (2*sh);
}


/* -----------------------------------------------
	calculate IJG Quality from scalefactor
	----------------------------------------------- */
int calcijgq( float scf )
{
	if ( ( scf >= 1.0 - EPS ) && ( scf <= 50.0 + EPS ) )
		return (int) ( 50.0 / scf + 0.5 );
	else if ( ( scf >= 0.01 - EPS ) && ( scf <= 1.0 + EPS ) )
		return (int) ( 100.0 * ( 1 - scf / 2.0 ) + 0.5 );
	else return 0;
}


/* -----------------------------------------------
	determine whether qtables are standard or custom
	----------------------------------------------- */
bool detquant( int colc, int ijgqua )
{
	int iscf;
	int qs[ 64 ];
	int i;

	if ( ijgqua <= 0 ) ijgqua = 1;
	else if ( ijgqua > 100 ) ijgqua = 100;

	if (ijgqua < 50) iscf = 5000 / ijgqua;
	else iscf = 200 - ijgqua*2;

	for ( i = 0; i < 64; i++ ) {
		qs[ i ] = ( std_qtables[ colc ][ unzigzag[i] ] * iscf + 50) / 100;
		if ( qs[ i ] <   1 ) qs[ i ] =   1;
		if ( qs[ i ] > 255 ) qs[ i ] = 255;
	}

	/*****
	printf("\n");
	printf("qsij = %2d %2d %2d %2d %2d\n", qs00, qs10, qs01, qs11, qs77);
	printf("cqij = %2d %2d %2d %2d %2d\n", colcnfo[colc].qtable[0], \
	 colcnfo[colc].qtable[1], \
	 colcnfo[colc].qtable[2], \
	 colcnfo[colc].qtable[4], \
	 colcnfo[colc].qtable[63] );
	printf("bool = %2d %2d %2d %2d %2d\n", colcnfo[colc].qtable[ 0] == qs00, \
	 colcnfo[colc].qtable[ 1] == qs10, \
	 colcnfo[colc].qtable[ 2] == qs01, \
	 colcnfo[colc].qtable[ 4] == qs11, \
	 colcnfo[colc].qtable[63] == qs77);
	****/

	if ( ijgqua == 100 ) {
		for ( i = 0; i < 64; i++ )
			if ( cmpnfo[colc].qtable[ i ] != 1 ) break;
		if ( i != 64 ) return false;
	}
	else {
		for ( i = 0; i < 64; i++ )
			if ( cmpnfo[colc].qtable[ i ] != qs[ i ] ) break;
		if ( i != 64 ) return false;
	}
	
	return true;
}


/* ----------------------- End of JPEG information functions -------------------------- */

/* ----------------------- Begin of miscellaneous helper functions -------------------------- */


/* -----------------------------------------------
	creates filename, callocs memory for it
	----------------------------------------------- */	
char* create_filename( char* base, char* extension )
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
void set_extension( char* destination, char* origin, char* extension )
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

/* -----------------------------------------------
	displays progress bar on screen
	----------------------------------------------- */
void progress_bar( int current, int last )
{
	int barpos = ( ( current * BARLEN ) + ( last / 2 ) ) / last;
	int i;
	
	
	// generate progress bar
	fprintf( msgout, "[" );
	#if !defined( UNIX )
	for ( i = 0; i < barpos; i++ )
		fprintf( msgout, "\xFE" );
	#else
	for ( i = 0; i < barpos; i++ )
		fprintf( msgout, "X" );
	#endif
	for (  ; i < BARLEN; i++ )
		fprintf( msgout, " " );
	fprintf( msgout, "]" );
}

/* ----------------------- End of miscellaneous helper functions -------------------------- */

/* ----------------------- End of file -------------------------- */
