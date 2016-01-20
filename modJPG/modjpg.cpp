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
#include "htables.h"


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

#define CSV_NAME	"JPEGinfo.csv"
#define MEM_ERRMSG	"out of memory error"
#define FRD_ERRMSG	"could not write file / file writeprotected: %s"
#define FWR_ERRMSG	"could not read file / file not found: %s"


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

struct huffCodes {
	unsigned short cval[ 256 ];
	unsigned short clen[ 256 ];
	unsigned short max_eobrun;
};

struct huffTree {
	unsigned short l[ 256 ];
	unsigned short r[ 256 ];
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
bool merge_jpeg( void );
bool decode_jpeg( void );
bool recode_jpeg( void );
bool mod_jpeg( void );
bool write_info( void );
bool reset_buffers( void );


/* -----------------------------------------------
	function declarations: jpeg-specific
	----------------------------------------------- */

bool setup_imginfo_jpg( void );
bool parse_jfif_jpg( unsigned char type, unsigned int len, unsigned char* segment );

int decode_block_seq( abitreader* huffr, huffTree* dctree, huffTree* actree, short* block );
int encode_block_seq( abitwriter* huffw, huffCodes* dctbl, huffCodes* actbl, short* block );

int decode_dc_prg_fs( abitreader* huffr, huffTree* dctree, short* block );
int encode_dc_prg_fs( abitwriter* huffw, huffCodes* dctbl, short* block );
int decode_ac_prg_fs( abitreader* huffr, huffTree* actree, short* block,
						int* eobrun, int from, int to );
int encode_ac_prg_fs( abitwriter* huffw, huffCodes* actbl, short* block,
						int* eobrun, int from, int to );

int decode_dc_prg_sa( abitreader* huffr, short* block );
int encode_dc_prg_sa( abitwriter* huffw, short* block );
int decode_ac_prg_sa( abitreader* huffr, huffTree* actree, short* block,
						int* eobrun, int from, int to );
int encode_ac_prg_sa( abitwriter* huffw, abytewriter* storw, huffCodes* actbl,
						short* block, int* eobrun, int from, int to );

int decode_eobrun_sa( abitreader* huffr, short* block, int* eobrun, int from, int to );
int encode_eobrun( abitwriter* huffw, huffCodes* actbl, int* eobrun );
int encode_crbits( abitwriter* huffw, abytewriter* storw );

int next_huffcode( abitreader *huffw, huffTree *ctree );
int next_mcupos( int* mcu, int* cmp, int* csc, int* sub, int* dpos, int* rstw );
int next_mcuposn( int* cmp, int* dpos, int* rstw );
int skip_eobrun( int* cmp, int* dpos, int* rstw, int* eobrun );

void build_huffcodes( unsigned char *clen, unsigned char *cval,
				huffCodes *hc, huffTree *ht );

/* -----------------------------------------------
	function declarations: JPG modifying
	----------------------------------------------- */

bool rebuild_header( void );
bool remove_rst( void );
bool recalc_qtable_div( float df, int tno );
bool reinsert_qtables( void );
bool check_std_hufftables( void );
float calc_quality( int tno );
void calc_quality_range( int tno, float* low, float* high );


/* -----------------------------------------------
	function declarations: miscelaneous helpers
	----------------------------------------------- */

char* create_filename( char* base, char* extension );
void set_extension( char* destination, char* origin, char* extension );
void add_underscore( char* filename );
int median_int( int* values, int size );
float median_float( float* values, int size );


/* -----------------------------------------------
	global variables: data storage
	----------------------------------------------- */

unsigned short qtables[4][64];				// quantization tables
huffCodes      hcodes[2][4];				// huffman codes
huffTree       htrees[2][4];				// huffman decoding trees
unsigned char  htset[2][4];					// 1 if huffman table is set

unsigned char* grbgdata			= 	NULL;	// garbage data
unsigned char* hdrdata          =   NULL;   // header data
unsigned char* huffdata         =   NULL;   // huffman coded data
int            hufs             =    0  ;   // size of huffman data
int            hdrs             =    0  ;   // size of header
int            grbs             =    0  ;   // size of garbage

unsigned int*  rstp             =   NULL;   // restart markers positions in huffdata
unsigned int*  scnp             =   NULL;   // scan start positions in huffdata
int            rstc             =    0  ;   // count of restart markers
int            scnc             =    0  ;   // count of scans
int            rsti             =    0  ;   // restart interval
char           padbit           =    -1 ;   // padbit (for huffman coding)

signed short*  colldata[4][64]  = { NULL }; // collection sorted DCT coefficients


/* -----------------------------------------------
	global variables: info about image
	----------------------------------------------- */

// seperate info for each color component
componentInfo cmpnfo[ 4 ];
float ijg_q0[ 4 ] = { -1 };
float ijg_q1[ 4 ] = { -1 };

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
	

/* -----------------------------------------------
	global variables: info about files
	----------------------------------------------- */
	
char*  jpgfilename;			// name of JPEG file
char*  modfilename;			// name of modified file
F_TYPE filetype;			// type of input file

int    jpgfilesize;			// size of JPEG file
int    modfilesize;			// size of modified file
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

int  verbosity  = 0;		// level of verbosity
bool overwrite  = false;	// overwrite files yes / no
int  err_tresh  = 1;		// error threshold ( proceed on warnings yes (2) / no (1) )

FILE*  msgout   = stdout;	// stream for output of messages
bool   pipe_on  = false;	// use stdin/stdout instead of filelist


/* -----------------------------------------------
	global variables: modification steps
	----------------------------------------------- */

bool  write_csv  = false;	// write info csv yes / no
bool  disc_meta  = false;	// discard meta-info yes / no
bool  disc_grbg  = false;	// discard garbage yes/no
bool  rem_rst    = false;	// remove restart markers yes/no
float div_qu[4]  = 			// quantization redividing steps
			{ 1, 1, 1, 1 };


/* -----------------------------------------------
	global variables: info about program
	----------------------------------------------- */

const unsigned char prgversion   = 9;
static const char*  subversion   = "";
static const char*  versiondate  = "12/05/2011";
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
	
	int acc_jpgsize = 0;
	int acc_modsize = 0;
	
	int speed, bpms;
	float cr;
	
	errorlevel = 0;
	
	
	// read options from command line
	initialize_options( argc, argv );
	
	// write program info to screen
	fprintf( msgout,  "\n--> modJPG v%i.%i%s (%s) by %s <--\n\n",
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
				acc_modsize += modfilesize;
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
		cr    = ( acc_jpgsize > 0 ) ? ( 100.0 * acc_modsize / acc_jpgsize ) : 0;
		
		fprintf( msgout,  " --------------------------------- \n" );
		fprintf( msgout,  " time taken        : %8i msec\n", speed );
		fprintf( msgout,  " avrg. byte per ms : %8i byte\n", bpms );
		fprintf( msgout,  " avrg. comp. ratio : %8.2f %%\n", cr );
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
	float tmp_flt;
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
			verbosity = ( verbosity < 0 ) ? 0 : verbosity;
			verbosity = ( verbosity > 2 ) ? 2 : verbosity;			
		}
		else if ( strcmp((*argv), "-o" ) == 0 ) {
			overwrite = true;
		}
		else if ( strcmp((*argv), "-p" ) == 0 ) {
			err_tresh = 2;
		}
		else if ( strcmp((*argv), "-csv" ) == 0 ) {
			write_csv = true;
		}
		else if ( sscanf( (*argv), "-dq%i,%f", &tmp_val, & tmp_flt ) == 2 ) {
			if ( ( tmp_val >= 0 ) && ( tmp_val < 4 ) )
				div_qu[ tmp_val ] = ( tmp_flt < 0 ) ? 0 : tmp_flt;
		}
		else if ( sscanf( (*argv), "-dq%f", &tmp_flt ) == 1 ) {
			for ( i = 0; i < 4; i++ )
				div_qu[ i ] = ( tmp_flt < 0 ) ? 0 : tmp_flt;
		}
		else if ( strcmp((*argv), "-dm" ) == 0 ) {
			disc_meta = true;
		}
		else if ( strcmp((*argv), "-dg" ) == 0 ) {
			disc_grbg = true;
		}
		else if ( strcmp((*argv), "-rr" ) == 0 ) {
			rem_rst = true;
		}
		else if ( strcmp((*argv), "-") == 0 ) {			
			msgout = stderr;
			// set binary mode for stdin & stdout
			#if !defined( unix )				
				setmode( fileno( stdin ), O_BINARY );
				setmode( fileno( stdout ), O_BINARY );
			#endif
			// use "-" as placeholder for stdin
			*(tmp_flp++) = (char*) "-";
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
	char* actionmsg  = NULL;
	char* errtypemsg = NULL;
	int speed, bpms;
	float cr;
	int i;
	
	
	errorfunction = NULL;
	errorlevel = 0;
	jpgfilesize = 0;
	modfilesize = 0;	
	
	
	// compare file name, set pipe if needed
	if ( strcmp( filelist[ file_no ], "-" ) == 0 ) {
		pipe_on = true;
		filelist[ file_no ] = (char*) "STDIN";
	}
		
	fprintf( msgout,  "\nProcessing file %i of %i \"%s\" -> ",
				file_no + 1, file_cnt, filelist[ file_no ] );
	if ( verbosity > 1 )
		fprintf( msgout,  "\n----------------------------------------" );
	
	// check input file and determine filetype
	execute( ( pipe_on ) ? check_stdin : check_file );
	
	// get specific action message
	if ( filetype == UNK ) actionmsg = (char*) "unknown filetype";
	else if ( !write_csv ) actionmsg = (char*) "Modifying";
	else actionmsg = (char*) "Writing info";
	
	if ( verbosity < 2 ) fprintf( msgout, "%s -> ", actionmsg );
	
	
	// main function routine
	begin = clock();
	
	if ( filetype == JPEG )
	{
		execute( read_jpeg );
		execute( decode_jpeg );
		if ( !write_csv ) {
			execute( mod_jpeg );
			execute( recode_jpeg );
			execute( merge_jpeg );
		}
		else {
			execute( write_info );
		}
	}
	// reset buffers
	reset_buffers();
	
	end = clock();
	
	
	// speed and compression ratio calculation
	speed = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );
	bpms  = ( speed > 0 ) ? ( jpgfilesize / speed ) : jpgfilesize;
	cr    = ( jpgfilesize > 0 ) ? ( 100.0 * modfilesize / jpgfilesize ) : 0;

	
	switch ( verbosity )
	{
		case 0:
			if ( errorlevel < err_tresh )
				if ( !write_csv ) fprintf( msgout,  "%.2f%%", cr );
				else fprintf( msgout,  "OK" );
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
		if ( verbosity > 1 )
			fprintf( stderr, " (in file \"%s\")\n", filelist[ file_no ] );
	}
	if ( (verbosity > 0) && (errorlevel < err_tresh) && (!write_csv) )
	{
		for ( i = 0; i < 4; i++ ) {
			if ( ijg_q0[ i ] > 0 ) {
				fprintf( msgout,  " qtable %i Q%% : %4.1f%% -> %4.1f%%\n", i, ijg_q0[ i ] * 100, ijg_q1[ i ] * 100 );
			}
		}
		fprintf( msgout,  " time taken  : %7i msec\n", speed );
		fprintf( msgout,  " byte per ms : %7i byte\n", bpms );
		fprintf( msgout,  " comp. ratio : %7.2f %%\n", cr );
	}
	
	if (verbosity > 1)
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
	else if ( function == *merge_jpeg ) {
		sprintf( statusmessage, "Merging header & image data" );
	}
	else if ( function == *decode_jpeg ) {
		sprintf( statusmessage, "Decompressing JPEG image data" );
	}
	else if ( function == *recode_jpeg ) {
		sprintf( statusmessage, "Recompressing JPEG image data" );
	}
	else if ( function == *mod_jpeg ) {
		sprintf( statusmessage, "Modifying JPEG data" );
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
	fprintf( msgout, "Usage: modjpg [switches] [filename(s)]");
	fprintf( msgout, "\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-v?]    set level of verbosity (max: 2) (def: 0)\n" );
	fprintf( msgout, " [-o]     overwrite existing files\n" );
	fprintf( msgout, " [-p]     proceed on warnings\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-dm]    discard meta-info\n" );
	fprintf( msgout, " [-dg]    discard garbage data\n" );
	fprintf( msgout, " [-rr]    remove restart markers\n" );
	fprintf( msgout, " [-dq?]   requantize coefficients with factor ?\n" );
	fprintf( msgout, " [-dq?,?] requantize table ? with factor ?\n" );
	fprintf( msgout, " [-csv]   write info to '%s'\n", CSV_NAME );
	fprintf( msgout, "\n" );
	fprintf( msgout, "Examples: \"modjpg -v1 -o -dm baboon.jpg\"\n" );
	fprintf( msgout, "          \"modjpg -dm -dg -dq3 *.jpg\"\n" );	
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
	if ( modfilename != NULL ) free( modfilename );
	// alloc memory for filenames
	jpgfilename = (char*) calloc( namelen, sizeof( char ) );
	modfilename = (char*) calloc( namelen, sizeof( char ) );
	
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
		// create filenames
		strcpy( jpgfilename, filelist[ file_no ] );
		strcpy( modfilename, filelist[ file_no ] );
		add_underscore( modfilename );
		if ( !overwrite ) {
			while ( access( modfilename, 0 ) == 0 ) {
				namelen += sizeof( char );
				modfilename = (char*) realloc( modfilename, namelen );
				add_underscore( modfilename );
			}
		}
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
		jpgfilename = (char*) "JPG file from stdin";
		modfilename = (char*) "JPG file to stdout";
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
	
	abytewriter* huffw;	
	abytewriter* hdrw;
	abytewriter* grbgw;	
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
	
	// start huffman writer
	huffw = new abytewriter( 0 );
	hufs  = 0; // size of image data, start with 0
	
	// alloc memory for segment data first
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
					break;				
				// check if 0xFF encountered
				if ( tmp == 0xFF ) {					
					fread( &tmp, 1, 1, fp ); // read next byte & check
					if ( tmp == 0x00 )
						// no zeroes needed -> ignore 0x00. write 0xFF
						huffw->write( 0xFF );
					else if ( ( tmp >= 0xD0 ) && ( tmp <= 0xD7 ) ) { // restart marker
						// leave marker out of data, but check next bytes
						fread( &tmp, 1, 1, fp ); // read next byte & check
						if ( tmp == 0xFF ) {							
							fread( &tmp, 1, 1, fp ); // check again
							if ( tmp == 0x00 )
								huffw->write( 0xFF );
							else { // put out warning and continue
								segment[ 0 ] = 0xFF;
								segment[ 1 ] = tmp;
								sprintf( errormessage, "RST marker found at end of scan" );
								errorlevel = 1;
								break;
							}
						}
						else // write current byte
							huffw->write( tmp );
					}
					else { // in all other cases leave it to the header parser routines
						segment[ 0 ] = 0xFF;
						segment[ 1 ] = tmp;
						break;
					}
				}
				else { // write current byte
					huffw->write( tmp );
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
					delete ( huffw );
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
			// get pointer for huffman data & size
			huffdata = huffw->getptr();
			hufs     = huffw->getpos();
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
				delete ( huffw );
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
	delete ( huffw );
	
	// check if everything went OK
	if ( ( hdrs == 0 ) || ( hufs == 0 ) ) {
		if ( !pipe_on ) fclose( fp );
		sprintf( errormessage, "unexpected end of data encountered" );
		errorlevel = 2;
		return false;
	}
	
	// store garbage after EOI if needed
	grbs = fread( &tmp, 1, 1, fp );	
	if ( grbs > 0 ) {
		// sprintf( errormessage, "data after EOI - last bytes: FF D9 %2X", tmp );
		// errorlevel = 1;		
		grbgw = new abytewriter( 1024 );
		grbgw->write( tmp );
		while( true ) {
			len = fread( segment, 1, ssize, fp );
			if ( len == 0 ) break;
			grbgw->write_n( segment, len );
		}
		grbgdata = grbgw->getptr();
		grbs     = grbgw->getpos();
		delete ( grbgw );
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
	Merges header & image data to jpeg
	----------------------------------------------- */
	
bool merge_jpeg( void )
{
	FILE* fp;
	
	unsigned char SOI[ 2 ] = { 0xFF, 0xD8 }; // SOI segment
	unsigned char EOI[ 2 ] = { 0xFF, 0xD9 }; // EOI segment
	unsigned char mrk = 0xFF; // marker start
	unsigned char stv = 0x00; // 0xFF stuff value
	unsigned char rst = 0xD0; // restart marker
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // current position in header
	unsigned int   ipos = 0; // current position in imagedata
	unsigned int   rpos = 0; // current restart marker position
	unsigned int   cpos = 0; // in scan corrected rst marker position
	unsigned int   scan = 1; // number of current scan
	unsigned int   tmp; // temporary storage variable
	
	
	if ( !pipe_on ) {
		// open jpeg for output
		fp = fopen( modfilename, "wb" );
		if ( fp == NULL ){
			sprintf( errormessage, FWR_ERRMSG, jpgfilename );
			errorlevel = 2;
			return false;
		}
	}
	else {
		// use stdout
		fp = stdout;
	}
	
	// write SOI
	fwrite( SOI, 1, 2, fp );
	
	// JPEG writing loop
	while ( true )
	{
		// store current header position
		tmp = hpos;
		
		// seek till start-of-scan
		for ( type = 0x00; type != 0xDA; ) {
			if ( ( int ) hpos >= hdrs ) break;
			type = hdrdata[ hpos + 1 ];
			len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
			hpos += len;
		}
		
		// write header data to file
		fwrite( hdrdata + tmp, 1, ( hpos - tmp ), fp );
		
		// get out if last marker segment type was not SOS
		if ( type != 0xDA ) break;
		
		
		// (re)set corrected rst pos
		cpos = 0;
		
		// write & expand huffman coded image data
		for ( ipos = scnp[ scan - 1 ]; ipos < scnp[ scan ]; ipos++ ) {
			// write current byte
			fwrite( huffdata + ipos, 1, 1, fp );
			// check current byte, stuff if needed
			if ( huffdata[ ipos ] == 0xFF )
				fwrite( &stv, 1, 1, fp );
			// insert restart markers if needed
			if ( rstp != NULL ) {
				if ( ipos == rstp[ rpos ] ) {
					rst = 0xD0 + ( cpos % 8 );
					fwrite( &mrk, 1, 1, fp );
					fwrite( &rst, 1, 1, fp );
					rpos++; cpos++;
				}
			}
		}

		// proceed with next scan
		scan++;
	}
	
	// write EOI
	fwrite( EOI, 1, 2, fp );
	
	// write garbage if needed
	if ( grbs > 0 )
		fwrite( grbgdata, 1, grbs, fp );
	
	// errormessage if write error
	if ( ferror( fp ) ) {
		if ( !pipe_on ) fclose( fp );
		sprintf( errormessage, "write error, possibly drive is full" );
		errorlevel = 2;		
		return false;
	}
	
	// get filesize
	fflush( fp );
	modfilesize = ftell( fp );
	
	// close file
	if ( !pipe_on )
		fclose( fp );	
	
	
	return true;
}


/* -----------------------------------------------
	JPEG decoding routine
	----------------------------------------------- */

bool decode_jpeg( void )
{
	abitreader* huffr; // bitwise reader for image data
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // current position in header
	
	int lastdc[ 4 ]; // last dc for each component
	short block[ 64 ]; // store block for coeffs
	int peobrun; // previous eobrun
	int eobrun; // run of eobs
	int rstw; // restart wait counter
	
	int cmp, bpos, dpos;
	int mcu, sub, csc;
	int eob, sta;
	
	
	// open huffman coded image data for input in abitreader
	huffr = new abitreader( huffdata, hufs );
	
	// preset count of scans
	scnc = 0;
	
	// JPEG decompression loop
	while ( true )
	{
		// seek till start-of-scan, parse only DHT, DRI and SOS
		for ( type = 0x00; type != 0xDA; ) {
			if ( ( int ) hpos >= hdrs ) break;
			type = hdrdata[ hpos + 1 ];
			len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
			if ( ( type == 0xC4 ) || ( type == 0xDA ) || ( type == 0xDD ) ) {
				if ( !parse_jfif_jpg( type, len, &( hdrdata[ hpos ] ) ) ) {
					return false;
				}
			}
			hpos += len;
		}
		
		// get out if last marker segment type was not SOS
		if ( type != 0xDA ) break;
		
		// check if huffman tables are available
		for ( csc = 0; csc < cs_cmpc; csc++ ) {
			cmp = cs_cmp[ csc ];
			if ( ( ( cs_sal == 0 ) && ( htset[ 0 ][ cmpnfo[cmp].huffdc ] == 0 ) ) ||
				 ( ( cs_sah >  0 ) && ( htset[ 1 ][ cmpnfo[cmp].huffac ] == 0 ) ) ) {
				sprintf( errormessage, "huffman table missing in scan%i", scnc );
				delete huffr;
				errorlevel = 2;
				return false;
			}
		}
		
		
		// intial variables set for decoding
		cmp  = cs_cmp[ 0 ];
		csc  = 0;
		mcu  = 0;
		sub  = 0;
		dpos = 0;
		
		// JPEG imagedata decoding routines
		while ( true )
		{			
			// (re)set last DCs for diff coding
			lastdc[ 0 ] = 0;
			lastdc[ 1 ] = 0;
			lastdc[ 2 ] = 0;
			lastdc[ 3 ] = 0;
			
			// (re)set status
			sta = 0;
			
			// (re)set eobrun
			eobrun  = 0;
			peobrun = 0;
			
			// (re)set rst wait counter
			rstw = rsti;
			
			// decoding for interleaved data
			if ( cs_cmpc > 1 )
			{				
				if ( jpegtype == 1 ) {
					// ---> sequential interleaved decoding <---
					while ( sta == 0 ) {
						// decode block
						eob = decode_block_seq( huffr,
							&(htrees[ 0 ][ cmpnfo[cmp].huffdc ]),
							&(htrees[ 1 ][ cmpnfo[cmp].huffdc ]),
							block );
						
						// fix dc
						block[ 0 ] += lastdc[ cmp ];
						lastdc[ cmp ] = block[ 0 ];
						
						// copy to colldata
						for ( bpos = 0; bpos < eob; bpos++ )
							colldata[ cmp ][ bpos ][ dpos ] = block[ bpos ];
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else if ( cs_sah == 0 ) {
					// ---> progressive interleaved DC decoding <---
					// ---> succesive approximation first stage <---
					while ( sta == 0 ) {
						sta = decode_dc_prg_fs( huffr,
							&(htrees[ 0 ][ cmpnfo[cmp].huffdc ]),
							block );
						
						// fix dc for diff coding
						colldata[cmp][0][dpos] = block[0] + lastdc[ cmp ];
						lastdc[ cmp ] = colldata[cmp][0][dpos];
						
						// bitshift for succesive approximation
						colldata[cmp][0][dpos] <<= cs_sal;
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else {
					// ---> progressive interleaved DC decoding <---
					// ---> succesive approximation later stage <---					
					while ( sta == 0 ) {
						// decode next bit
						sta = decode_dc_prg_sa( huffr,
							block );
						
						// shift in next bit
						colldata[cmp][0][dpos] += block[0] << cs_sal;
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
			}
			else // decoding for non interleaved data
			{
				if ( jpegtype == 1 ) {
					// ---> sequential non interleaved decoding <---
					while ( sta == 0 ) {
						// decode block
						eob = decode_block_seq( huffr,
							&(htrees[ 0 ][ cmpnfo[cmp].huffdc ]),
							&(htrees[ 1 ][ cmpnfo[cmp].huffdc ]),
							block );
						
						// fix dc
						block[ 0 ] += lastdc[ cmp ];
						lastdc[ cmp ] = block[ 0 ];
						
						// copy to colldata
						for ( bpos = 0; bpos < eob; bpos++ )
							colldata[ cmp ][ bpos ][ dpos ] = block[ bpos ];
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = next_mcuposn( &cmp, &dpos, &rstw );
					}
				}
				else if ( cs_to == 0 ) {					
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved DC decoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							sta = decode_dc_prg_fs( huffr,
								&(htrees[ 0 ][ cmpnfo[cmp].huffdc ]),
								block );
								
							// fix dc for diff coding
							colldata[cmp][0][dpos] = block[0] + lastdc[ cmp ];
							lastdc[ cmp ] = colldata[cmp][0][dpos];
							
							// bitshift for succesive approximation
							colldata[cmp][0][dpos] <<= cs_sal;
							
							// check for errors, increment dpos otherwise
							if ( sta != -1 )
								sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
					else {
						// ---> progressive non interleaved DC decoding <---
						// ---> succesive approximation later stage <---
						while( sta == 0 ) {
							// decode next bit
							sta = decode_dc_prg_sa( huffr,
								block );
							
							// shift in next bit
							colldata[cmp][0][dpos] += block[0] << cs_sal;
							
							// check for errors, increment dpos otherwise
							if ( sta != -1 )
								sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
				}
				else {
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved AC decoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							// decode block
							eob = decode_ac_prg_fs( huffr,
								&(htrees[ 1 ][ cmpnfo[cmp].huffac ]),
								block, &eobrun, cs_from, cs_to );
							
							// check for non optimal coding
							if ( ( eob == cs_from ) && ( eobrun > 0 ) &&
								( peobrun > 0 ) && ( peobrun <
								hcodes[ 1 ][ cmpnfo[cmp].huffac ].max_eobrun - 1 ) ) {
								sprintf( errormessage,
									"reconstruction of non optimal coding not supported" );
								errorlevel = 1;
							}
							
							// copy to colldata
							for ( bpos = cs_from; bpos < eob; bpos++ )
								colldata[ cmp ][ bpos ][ dpos ] = block[ bpos ] << cs_sal;
							
							// check for errors
							if ( eob < 0 ) sta = -1;
							else sta = skip_eobrun( &cmp, &dpos, &rstw, &eobrun );
							
							// proceed only if no error encountered
							if ( sta == 0 )
								sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
					else {
						// ---> progressive non interleaved AC decoding <---
						// ---> succesive approximation later stage <---
						while ( sta == 0 ) {
							// copy from colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								block[ bpos ] = colldata[ cmp ][ bpos ][ dpos ];
							
							if ( eobrun == 0 ) {
								// decode block (long routine)
								eob = decode_ac_prg_sa( huffr,
									&(htrees[ 1 ][ cmpnfo[cmp].huffac ]),
									block, &eobrun, cs_from, cs_to );
								
								// check for non optimal coding
								if ( ( eob == cs_from ) && ( eobrun > 0 ) &&
									( peobrun > 0 ) && ( peobrun <
									hcodes[ 1 ][ cmpnfo[cmp].huffac ].max_eobrun - 1 ) ) {
									sprintf( errormessage,
										"reconstruction of non optimal coding not supported" );
									errorlevel = 1;
								}
								
								// store eobrun
								peobrun = eobrun;
							}
							else {
								// decode block (short routine)
								eob = decode_eobrun_sa( huffr,
									block, &eobrun, cs_from, cs_to );
							}
								
							// copy back to colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								colldata[ cmp ][ bpos ][ dpos ] += block[ bpos ] << cs_sal;
							
							// proceed only if no error encountered
							if ( eob < 0 ) sta = -1;
							else sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
				}
			}			
			
			// unpad huffman reader / check padbit
			if ( padbit != -1 ) {
				if ( padbit != huffr->unpad( padbit ) ) {
					sprintf( errormessage, "inconsistent use of padbits" );
					padbit = 1;
					errorlevel = 1;
				}
			}
			else {
				padbit = huffr->unpad( padbit );
			}
			
			// evaluate status
			if ( sta == -1 ) { // status -1 means error
				sprintf( errormessage, "decode error in scan%i / mcu%i",
					scnc, ( cs_cmpc > 1 ) ? mcu : dpos );
				delete huffr;
				errorlevel = 2;
				return false;
			}
			else if ( sta == 2 ) { // status 2/3 means done
				scnc++; // increment scan counter
				break; // leave decoding loop, everything is done here
			}
			// else if ( sta == 1 ); // status 1 means restart - so stay in the loop
		}
	}
	
	// check for unneeded data
	if ( !huffr->eof ) {
		sprintf( errormessage, "unneeded data found after coded image data" );
		errorlevel = 1;
	}
			
	// clean up
	delete( huffr );
	
	
	return true;
}


/* -----------------------------------------------
	JPEG encoding routine
	----------------------------------------------- */

bool recode_jpeg( void )
{
	abitwriter*  huffw; // bitwise writer for image data
	abytewriter* storw; // bytewise writer for storage of correction bits
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // current position in header
		
	int lastdc[ 4 ]; // last dc for each component
	short block[ 64 ]; // store block for coeffs
	int eobrun; // run of eobs
	int rstw; // restart wait counter
	
	int cmp, bpos, dpos;
	int mcu, sub, csc;
	int eob, sta;
	int tmp;
	
	
	// open huffman coded image data in abitwriter
	huffw = new abitwriter( 0 );
	huffw->fillbit = padbit;
	
	// init storage writer
	storw = new abytewriter( 0 );
	
	// preset count of scans and restarts
	scnc = 0;
	rstc = 0;
	
	// JPEG decompression loop
	while ( true )
	{
		// seek till start-of-scan, parse only DHT, DRI and SOS
		for ( type = 0x00; type != 0xDA; ) {
			if ( ( int ) hpos >= hdrs ) break;
			type = hdrdata[ hpos + 1 ];
			len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
			if ( ( type == 0xC4 ) || ( type == 0xDA ) || ( type == 0xDD ) ) {
				if ( !parse_jfif_jpg( type, len, &( hdrdata[ hpos ] ) ) ) {
					return false;
				}
				hpos += len;
			}
			else {
				hpos += len;
				continue;
			}			
		}
		
		// get out if last marker segment type was not SOS
		if ( type != 0xDA ) break;
		
		
		// (re)alloc scan positons array
		if ( scnp == NULL ) scnp = ( unsigned int* ) calloc( scnc + 2, sizeof( int ) );
		else scnp = ( unsigned int* ) realloc( scnp, ( scnc + 2 ) * sizeof( int ) );
		if ( scnp == NULL ) {
			sprintf( errormessage, MEM_ERRMSG );
			errorlevel = 2;
			return false;
		}
		
		// (re)alloc restart marker positons array if needed
		if ( rsti > 0 ) {
			tmp = rstc + ( ( cs_cmpc > 1 ) ?
				( mcuc / rsti ) : ( cmpnfo[ cs_cmp[ 0 ] ].bc / rsti ) );
			if ( rstp == NULL ) rstp = ( unsigned int* ) calloc( tmp + 1, sizeof( int ) );
			else rstp = ( unsigned int* ) realloc( rstp, ( tmp + 1 ) * sizeof( int ) );
			if ( rstp == NULL ) {
				sprintf( errormessage, MEM_ERRMSG );
				errorlevel = 2;
				return false;
			}
		}		
		
		// intial variables set for encoding
		cmp  = cs_cmp[ 0 ];
		csc  = 0;
		mcu  = 0;
		sub  = 0;
		dpos = 0;
		
		// store scan position
		scnp[ scnc ] = huffw->getpos();
		
		// JPEG imagedata encoding routines
		while ( true )
		{
			// (re)set last DCs for diff coding
			lastdc[ 0 ] = 0;
			lastdc[ 1 ] = 0;
			lastdc[ 2 ] = 0;
			lastdc[ 3 ] = 0;
			
			// (re)set status
			sta = 0;
			
			// (re)set eobrun
			eobrun = 0;
			
			// (re)set rst wait counter
			rstw = rsti;
			
			// encoding for interleaved data
			if ( cs_cmpc > 1 )
			{				
				if ( jpegtype == 1 ) {
					// ---> sequential interleaved encoding <---
					while ( sta == 0 ) {
						// copy from colldata
						for ( bpos = 0; bpos < 64; bpos++ )
							block[ bpos ] = colldata[ cmp ][ bpos ][ dpos ];
						
						// diff coding for dc
						block[ 0 ] -= lastdc[ cmp ];
						lastdc[ cmp ] = colldata[ cmp ][ 0 ][ dpos ];
						
						// encode block
						eob = encode_block_seq( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffac ]),
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							block );
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else if ( cs_sah == 0 ) {
					// ---> progressive interleaved DC encoding <---
					// ---> succesive approximation first stage <---
					while ( sta == 0 ) {
						// diff coding & bitshifting for dc 
						tmp = colldata[ cmp ][ 0 ][ dpos ] >> cs_sal;
						block[ 0 ] = tmp - lastdc[ cmp ];
						lastdc[ cmp ] = tmp;
						
						// encode dc
						sta = encode_dc_prg_fs( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffdc ]),
							block );
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else {
					// ---> progressive interleaved DC encoding <---
					// ---> succesive approximation later stage <---
					while ( sta == 0 ) {
						// fetch bit from current bitplane
						block[ 0 ] = BITN( colldata[ cmp ][ 0 ][ dpos ], cs_sal );
						
						// encode dc correction bit
						sta = encode_dc_prg_sa( huffw, block );
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
			}
			else // encoding for non interleaved data
			{
				if ( jpegtype == 1 ) {
					// ---> sequential non interleaved encoding <---
					while ( sta == 0 ) {
						// copy from colldata
						for ( bpos = 0; bpos < 64; bpos++ )
							block[ bpos ] = colldata[ cmp ][ bpos ][ dpos ];
						
						// diff coding for dc
						block[ 0 ] -= lastdc[ cmp ];
						lastdc[ cmp ] = colldata[ cmp ][ 0 ][ dpos ];
						
						// encode block
						eob = encode_block_seq( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffac ]),
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							block );
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = next_mcuposn( &cmp, &dpos, &rstw );	
					}
				}
				else if ( cs_to == 0 ) {
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved DC encoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							// diff coding & bitshifting for dc 
							tmp = colldata[ cmp ][ 0 ][ dpos ] >> cs_sal;
							block[ 0 ] = tmp - lastdc[ cmp ];
							lastdc[ cmp ] = tmp;
							
							// encode dc
							sta = encode_dc_prg_fs( huffw,
								&(hcodes[ 0 ][ cmpnfo[cmp].huffdc ]),
								block );							
							
							// check for errors, increment dpos otherwise
							if ( sta != -1 )
								sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
					else {
						// ---> progressive non interleaved DC encoding <---
						// ---> succesive approximation later stage <---
						while ( sta == 0 ) {
							// fetch bit from current bitplane
							block[ 0 ] = BITN( colldata[ cmp ][ 0 ][ dpos ], cs_sal );
							
							// encode dc correction bit
							sta = encode_dc_prg_sa( huffw, block );
							
							// next mcupos if no error happened
							if ( sta != -1 )
								sta = next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
				}
				else {
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved AC encoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							// copy from colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								block[ bpos ] =
									FDIV2( colldata[ cmp ][ bpos ][ dpos ], cs_sal );
							
							// encode block
							eob = encode_ac_prg_fs( huffw,
								&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
								block, &eobrun, cs_from, cs_to );
							
							// check for errors, proceed if no error encountered
							if ( eob < 0 ) sta = -1;
							else sta = next_mcuposn( &cmp, &dpos, &rstw );
						}						
						
						// encode remaining eobrun
						encode_eobrun( huffw,
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							&eobrun );
					}
					else {
						// ---> progressive non interleaved AC encoding <---
						// ---> succesive approximation later stage <---
						while ( sta == 0 ) {
							// copy from colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								block[ bpos ] =
									FDIV2( colldata[ cmp ][ bpos ][ dpos ], cs_sal );
							
							// encode block
							eob = encode_ac_prg_sa( huffw, storw,
								&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
								block, &eobrun, cs_from, cs_to );
							
							// check for errors, proceed if no error encountered
							if ( eob < 0 ) sta = -1;
							else sta = next_mcuposn( &cmp, &dpos, &rstw );
						}						
						
						// encode remaining eobrun
						encode_eobrun( huffw,
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							&eobrun );
							
						// encode remaining correction bits
						encode_crbits( huffw, storw );
					}
				}
			}
			
			// pad huffman writer
			huffw->pad( padbit );
			
			// evaluate status
			if ( sta == -1 ) { // status -1 means error
				sprintf( errormessage, "encode error in scan%i / mcu%i",
					scnc, ( cs_cmpc > 1 ) ? mcu : dpos );
				delete huffw;
				errorlevel = 2;
				return false;
			}
			else if ( sta == 2 ) { // status 2 means done
				scnc++; // increment scan counter
				break; // leave decoding loop, everything is done here
			}
			else if ( sta == 1 ) { // status 1 means restart
				if ( rsti > 0 ) // store rstp & stay in the loop
					rstp[ rstc++ ] = huffw->getpos() - 1;
			}
		}
	}
	
	// safety check for error in huffwriter
	if ( huffw->error ) {
		delete huffw;
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// get data into huffdata
	huffdata = huffw->getptr();
	hufs = huffw->getpos();	
	delete huffw;
	
	// remove storage writer
	delete storw;
	
	// store last scan & restart positions
	scnp[ scnc ] = hufs;
	if ( rstp != NULL )
		rstp[ rstc ] = hufs;
	
	
	return true;
}


/* -----------------------------------------------
	modify JPEGs
	----------------------------------------------- */
	
bool mod_jpeg( void )
{
	bool std_htables;
	int i;
	
	// calculate quality (before)
	for ( i = 0; i < 4; i++ )
		ijg_q0[ i ] = calc_quality( i );
	
	// check if standard tables are used
	std_htables = check_std_hufftables();
	
	// discard meta-info
	if ( disc_meta ) {
		rebuild_header();
	}
	// discard garbage
	if ( disc_grbg ) {
		if ( grbgdata != NULL ) free ( grbgdata );
		grbgdata = NULL;
		grbs = 0;
	}
	// recalculate quantization
	for ( i = 0; i < 4; i++ ) {
		// if ( div_qu[ i ] <= 1 ) continue;
		recalc_qtable_div( div_qu[ i ], i );
		if ( !std_htables ) {
			sprintf( errormessage, "requantizing is only safe if standard huffman tables are used" );
			errorlevel = 1;
		}
	}
	reinsert_qtables();
	
	// calculate quality (after)
	for ( i = 0; i < 4; i++ )
		ijg_q1[ i ] = calc_quality( i );
	
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
		comp_ratio += cmpnfo[ cmp ].bc;
		for ( i = 0; i < 4; i++ )
			if ( cmpnfo[ cmp ].qtable == qtables[ i ] ) break;
		if ( i < 4 ) calc_quality_range( i, &(ijg_quality_l[cmp]), &(ijg_quality_h[cmp]) );
	}
	comp_ratio = ( comp_ratio * 64 ) / jpgfilesize;
	
	if ( new_file ) { // write collumn titles only if file is empty
		fprintf( fp, "filename; " );
		fprintf( fp, "width; " );
		fprintf( fp, "height; " );
		fprintf( fp, "file size; " );
		fprintf( fp, "header size; " );
		fprintf( fp, "comp. ratio; " );
		fprintf( fp, "jpeg mode; " );
		fprintf( fp, "htables; " );
		fprintf( fp, "cmp#; " );
		fprintf( fp, "samples (0); " ); ;
		fprintf( fp, "samples (1); " ); ;
		#if defined( VERBOSE_CSV )
		fprintf( fp, "samples (2); " ); ;
		fprintf( fp, "samples (3); " ); ;
		#endif
		fprintf( fp, "quality (0); " ); ;
		fprintf( fp, "quality (1); " ); ;
		#if defined( VERBOSE_CSV )
		fprintf( fp, "quality (2); " ); ;
		fprintf( fp, "quality (3); " ); ;
		#endif		
		#if defined( VERBOSE_CSV )
		fprintf( fp, "quantDC (0); " ); ;
		fprintf( fp, "quantDC (1); " ); ;
		fprintf( fp, "quantDC (2); " ); ;
		fprintf( fp, "quantDC (3); " ); ;
		#endif
		fprintf( fp, "\n" );
	}
	
	// write data to file
	fprintf( fp, "%s; ", jpgfilename );
	fprintf( fp, "%i; ", imgwidth );
	fprintf( fp, "%i; ", imgheight );
	fprintf( fp, "%i; ", jpgfilesize );
	fprintf( fp, "%i; ", hdrs );
	fprintf( fp, "%.0f to 1; ", comp_ratio );
	fprintf( fp, "%s; ", ( jpegtype == 0 ) ? "unk" : ( jpegtype == 1 ) ? "seq" : "prg"  );
	fprintf( fp, "%s; ", ( std_htables ) ? "std" : "opt" );
	fprintf( fp, "%i; ", cmpc );
	fprintf( fp, "%ix%i; ", ( cmpc > 0 ) ? cmpnfo[0].sfh : 0, ( cmpc > 0 ) ? cmpnfo[0].sfv : 0 );
	fprintf( fp, "%ix%i; ", ( cmpc > 1 ) ? cmpnfo[1].sfh : 0, ( cmpc > 1 ) ? cmpnfo[1].sfv : 0 );
	#if defined( VERBOSE_CSV )
	fprintf( fp, "%ix%i; ", ( cmpc > 2 ) ? cmpnfo[2].sfh : 0, ( cmpc > 2 ) ? cmpnfo[2].sfv : 0 );
	fprintf( fp, "%ix%i; ", ( cmpc > 3 ) ? cmpnfo[3].sfh : 0, ( cmpc > 3 ) ? cmpnfo[3].sfv : 0 );
	#endif
	fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 0 ) ? ijg_quality_l[0] * 100.0 : 0, ( cmpc > 0 ) ? ijg_quality_h[0] * 100.0 : 0 );
	fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 1 ) ? ijg_quality_l[1] * 100.0 : 0, ( cmpc > 1 ) ? ijg_quality_h[1] * 100.0 : 0 );
	#if defined( VERBOSE_CSV )
	fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 2 ) ? ijg_quality_l[2] * 100.0 : 0, ( cmpc > 2 ) ? ijg_quality_h[2] * 100.0 : 0 );
	fprintf( fp, "%.1f%%...%.1f%%; ", ( cmpc > 3 ) ? ijg_quality_l[3] * 100.0 : 0, ( cmpc > 3 ) ? ijg_quality_h[3] * 100.0 : 0 );
	#endif
	#if defined( VERBOSE_CSV )
	fprintf( fp, "%i; ", ( cmpc > 0 ) ? cmpnfo[0].qtable[0] : 0 );
	fprintf( fp, "%i; ", ( cmpc > 1 ) ? cmpnfo[1].qtable[0] : 0 );	
	fprintf( fp, "%i; ", ( cmpc > 2 ) ? cmpnfo[2].qtable[0] : 0 );
	fprintf( fp, "%i; ", ( cmpc > 3 ) ? cmpnfo[3].qtable[0] : 0 );
	#endif
	fprintf( fp, "\n" );
	
	
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
	if ( huffdata != NULL ) free ( huffdata );
	if ( grbgdata != NULL ) free ( grbgdata );
	if ( rstp     != NULL ) free ( rstp );
	if ( scnp     != NULL ) free ( scnp );
	hdrdata   = NULL;
	huffdata  = NULL;
	grbgdata  = NULL;
	rstp      = NULL;
	scnp      = NULL;
	
	// free image arrays
	for ( cmp = 0; cmp < 4; cmp++ )	{
		for ( bpos = 0; bpos < 64; bpos++ ) {
			if (colldata[ cmp ][ bpos ] != NULL) free( colldata[cmp][bpos] );
			colldata[ cmp ][ bpos ] = NULL;
		}		
	}
	
	
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
	
	// reset quantization / huffman tables
	for ( i = 0; i < 4; i++ ) {
		htset[ 0 ][ i ] = 0;
		htset[ 1 ][ i ] = 0;
		for ( bpos = 0; bpos < 64; bpos++ )
			qtables[ i ][ bpos ] = 0;
	}
	
	// preset jpegtype
	jpegtype  = 0;
	
	// reset padbit
	padbit = -1;
	
	
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
	
	int cmp, bpos;
	int mem_usage = hdrs + hufs + grbs; // hack
	
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
	
	// alloc memory for further operations
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		// alloc memory for colls
		for ( bpos = 0; bpos < 64; bpos++ ) {
			colldata[cmp][bpos] = (short int*) calloc ( cmpnfo[cmp].bc, sizeof( short ) );
			if (colldata[cmp][bpos] == NULL) {
				// sprintf( errormessage, MEM_ERRMSG );
				sprintf( errormessage, "out of mem (%i -> %i)", mem_usage, mem_usage + ( cmpnfo[cmp].bc * sizeof( short ) ) ); // hack
				errorlevel = 2;
				return false;
			}
			mem_usage += ( cmpnfo[cmp].bc * sizeof( short ) ); // hack
		}
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
	int skip;
	int cmp;
	int i;
	
	
	switch ( type )
	{
		case 0xC4: // DHT segment
			// build huffman trees & codes
			while ( hpos < len ) {
				lval = LBITS( segment[ hpos ], 4 );
				rval = RBITS( segment[ hpos ], 4 );
				if ( ((lval < 0) || (lval >= 2)) || ((rval < 0) || (rval >= 4)) )
					break;
					
				hpos++;
				// build huffman codes & trees
				build_huffcodes( &(segment[ hpos + 0 ]), &(segment[ hpos + 16 ]),
					&(hcodes[ lval ][ rval ]), &(htrees[ lval ][ rval ]) );
				htset[ lval ][ rval ] = 1;
				
				skip = 16;
				for ( i = 0; i < 16; i++ )		
					skip += ( int ) segment[ hpos + i ];				
				hpos += skip;
			}
			
			if ( hpos != len ) {
				// if we get here, something went wrong
				sprintf( errormessage, "size mismatch in dht marker" );
				errorlevel = 2;
				return false;
			}
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


/* -----------------------------------------------
	sequential block decoding routine
	----------------------------------------------- */
int decode_block_seq( abitreader* huffr, huffTree* dctree, huffTree* actree, short* block )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int eob = 64;
	int bpos;
	int hc;
	
	
	// decode dc
	hc = next_huffcode( huffr, dctree );
	if ( hc < 0 ) return -1; // return error
	else s = ( unsigned char ) hc;
	n = huffr->read( s );	
	block[ 0 ] = DEVLI( s, n );
	
	// decode ac
	for ( bpos = 1; bpos < 64; )
	{
		// decode next
		hc = next_huffcode( huffr, actree );
		// analyse code
		if ( hc > 0 ) {
			z = LBITS( hc, 4 );
			s = RBITS( hc, 4 );
			n = huffr->read( s );
			if ( ( z + bpos ) >= 64 )
				return -1; // run is to long
			while ( z > 0 ) { // write zeroes
				block[ bpos++ ] = 0;
				z--;
			}
			block[ bpos++ ] = ( short ) DEVLI( s, n ); // decode cvli
		}
		else if ( hc == 0 ) { // EOB
			eob = bpos;
			// while( bpos < 64 ) // fill remaining block with zeroes
			//	block[ bpos++ ] = 0;
			break;
		}
		else {
			return -1; // return error
		}
	}
	
	
	// return position of eob
	return eob;
}


/* -----------------------------------------------
	sequential block encoding routine
	----------------------------------------------- */
int encode_block_seq( abitwriter* huffw, huffCodes* dctbl, huffCodes* actbl, short* block )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int bpos;
	int hc;
	int tmp;
	
	
	// encode DC	
	tmp = ABS( block[ 0 ] );
	BITLEN( s, tmp );
	n = ENVLI( s, block[ 0 ] );
	huffw->write( dctbl->cval[ s ], dctbl->clen[ s ] );
	huffw->write( n, s );
	
	// encode AC
	z = 0;
	for ( bpos = 1; bpos < 64; bpos++ )
	{
		// if nonzero is encountered
		if ( block[ bpos ] != 0 ) {
			// write remaining zeroes
			while ( z >= 16 ) {
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				z -= 16;
			}			
			// vli encode
			tmp = ABS( block[ bpos ] );
			BITLEN( s, tmp );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// reset zeroes
			z = 0;
		}
		else { // increment zero counter
			z++;
		}
	}
	// write eob if needed
	if ( z > 0 )
		huffw->write( actbl->cval[ 0x00 ], actbl->clen[ 0x00 ] );
		
	
	return 64 - z;
}


/* -----------------------------------------------
	progressive DC decoding routine
	----------------------------------------------- */
int decode_dc_prg_fs( abitreader* huffr, huffTree* dctree, short* block )
{
	unsigned short n;
	unsigned char  s;
	int hc;
	
	
	// decode dc
	hc = next_huffcode( huffr, dctree );
	if ( hc < 0 ) return -1; // return error
	else s = ( unsigned char ) hc;
	n = huffr->read( s );	
	block[ 0 ] = DEVLI( s, n );
	
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive DC encoding routine
	----------------------------------------------- */
int encode_dc_prg_fs( abitwriter* huffw, huffCodes* dctbl, short* block )
{
	unsigned short n;
	unsigned char  s;
	int tmp;
	
	
	// encode DC	
	tmp = ABS( block[ 0 ] );
	BITLEN( s, tmp );
	n = ENVLI( s, block[ 0 ] );
	huffw->write( dctbl->cval[ s ], dctbl->clen[ s ] );
	huffw->write( n, s );
	
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive AC decoding routine
	----------------------------------------------- */
int decode_ac_prg_fs( abitreader* huffr, huffTree* actree, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int eob = to + 1;
	int bpos;
	int hc;
	int l;
	int r;
	
	
	// check eobrun
	if ( (*eobrun) > 0 ) {
		for ( bpos = from; bpos <= to; )
			block[ bpos ] = 0;
		(*eobrun)--;
		return from;
	}
	
	// decode ac
	for ( bpos = from; bpos <= to; )
	{
		// decode next
		hc = next_huffcode( huffr, actree );
		if ( hc < 0 ) return -1;
		l = LBITS( hc, 4 );
		r = RBITS( hc, 4 );
		// analyse code
		if ( ( l == 15 ) || ( r > 0 ) ) { // decode run/level combination
			z = l;
			s = r;
			n = huffr->read( s );
			if ( ( z + bpos ) > to )
				return -1; // run is to long			
			while ( z > 0 ) { // write zeroes
				block[ bpos++ ] = 0;
				z--;
			}			
			block[ bpos++ ] = ( short ) DEVLI( s, n ); // decode cvli
		}
		else { // decode eobrun
			eob = bpos;
			s = l;
			n = huffr->read( s );
			(*eobrun) = E_DEVLI( s, n );			
			// while( bpos <= to ) // fill remaining block with zeroes
			//	block[ bpos++ ] = 0;
			(*eobrun)--; // decrement eobrun ( for this one )
			break;
		}
	}
	
	
	// return position of eob
	return eob;
}


/* -----------------------------------------------
	progressive AC encoding routine
	----------------------------------------------- */
int encode_ac_prg_fs( abitwriter* huffw, huffCodes* actbl, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int bpos;
	int hc;
	int tmp;
	
	// encode AC
	z = 0;
	for ( bpos = from; bpos <= to; bpos++ )
	{
		// if nonzero is encountered
		if ( block[ bpos ] != 0 ) {
			// encode eobrun
			encode_eobrun( huffw, actbl, eobrun );
			// write remaining zeroes
			while ( z >= 16 ) {
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				z -= 16;
			}			
			// vli encode
			tmp = ABS( block[ bpos ] );
			BITLEN( s, tmp );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// reset zeroes
			z = 0;
		}
		else { // increment zero counter
			z++;
		}
	}
	
	// check eob, increment eobrun if needed
	if ( z > 0 ) {
		(*eobrun)++;
		// check eobrun, encode if needed
		if ( (*eobrun) == actbl->max_eobrun )
			encode_eobrun( huffw, actbl, eobrun );
		return 1 + to - z;		
	}
	else {
		return 1 + to;
	}
}


/* -----------------------------------------------
	progressive DC SA decoding routine
	----------------------------------------------- */
int decode_dc_prg_sa( abitreader* huffr, short* block )
{
	// decode next bit of dc coefficient
	block[ 0 ] = huffr->read( 1 );
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive DC SA encoding routine
	----------------------------------------------- */
int encode_dc_prg_sa( abitwriter* huffw, short* block )
{
	// enocode next bit of dc coefficient
	huffw->write( block[ 0 ], 1 );
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive AC SA decoding routine
	----------------------------------------------- */
int decode_ac_prg_sa( abitreader* huffr, huffTree* actree, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	signed char    z;
	signed char    v;
	int bpos = from;
	int eob = to;
	int hc;
	int l;
	int r;
	
	
	// decode AC succesive approximation bits
	if ( (*eobrun) == 0 )
	while ( bpos <= to )
	{
		// decode next
		hc = next_huffcode( huffr, actree );
		if ( hc < 0 ) return -1;
		l = LBITS( hc, 4 );
		r = RBITS( hc, 4 );
		// analyse code
		if ( ( l == 15 ) || ( r > 0 ) ) { // decode run/level combination
			z = l;
			s = r;
			if ( s == 0 ) v = 0;
			else if ( s == 1 ) {
				n = huffr->read( 1 );
				v = ( n == 0 ) ? -1 : 1; // fast decode vli
			}
			else return -1; // decoding error
			// write zeroes / write correction bits
			while ( true ) {
				if ( block[ bpos ] == 0 ) { // skip zeroes / write value
					if ( z > 0 ) z--;
					else {
						block[ bpos++ ] = v;
						break;
					}
				}
				else { // read correction bit
					n = huffr->read( 1 );
					block[ bpos ] = ( block[ bpos ] > 0 ) ? n : -n;
				}
				if ( bpos++ >= to ) return -1; // error check					
			}
		}
		else { // decode eobrun
			eob = bpos;
			s = l;
			n = huffr->read( s );
			(*eobrun) = E_DEVLI( s, n );
			break;
		}
	}
	
	// read after eob correction bits
	if ( (*eobrun) > 0 ) {
		for ( ; bpos <= to; bpos++ ) {
			if ( block[ bpos ] != 0 ) {
				n = huffr->read( 1 );
				block[ bpos ] = ( block[ bpos ] > 0 ) ? n : -n;
			}
		}
		// decrement eobrun
		(*eobrun)--;
	}
	
	// return eob
	return eob;
}


/* -----------------------------------------------
	progressive AC SA encoding routine
	----------------------------------------------- */
int encode_ac_prg_sa( abitwriter* huffw, abytewriter* storw, huffCodes* actbl, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int eob = from;
	int bpos;
	int hc;
	int tmp;
	
	// check if block contains any newly nonzero coefficients and find out position of eob
	for ( bpos = to; bpos >= from; bpos-- )	{
		if ( ( block[ bpos ] == 1 ) || ( block[ bpos ] == -1 ) ) {
			eob = bpos + 1;
			break;
		}
	}
	
	// encode eobrun if needed
	if ( ( eob > from ) && ( (*eobrun) > 0 ) ) {
		encode_eobrun( huffw, actbl, eobrun );
		encode_crbits( huffw, storw );
	}
	
	// encode AC
	z = 0;
	for ( bpos = from; bpos < eob; bpos++ )
	{
		// if zero is encountered
		if ( block[ bpos ] == 0 ) {
			z++; // increment zero counter
			if ( z == 16 ) { // write zeroes if needed
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				encode_crbits( huffw, storw );
				z = 0;
			}
		}
		// if nonzero is encountered
		else if ( ( block[ bpos ] == 1 ) || ( block[ bpos ] == -1 ) ) {
			// vli encode
			tmp = ABS( block[ bpos ] );
			BITLEN( s, tmp );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// write correction bits
			encode_crbits( huffw, storw );
			// reset zeroes
			z = 0;
		}
		else { // store correction bits
			n = block[ bpos ] & 0x1;
			storw->write( n );
		}
	}
	
	// fast processing after eob
	for ( ;bpos <= to; bpos++ )
	{
		if ( block[ bpos ] != 0 ) { // store correction bits
			n = block[ bpos ] & 0x1;
			storw->write( n );
		}
	}
	
	// check eob, increment eobrun if needed
	if ( eob <= to ) {
		(*eobrun)++;	
		// check eobrun, encode if needed
		if ( (*eobrun) == actbl->max_eobrun ) {
			encode_eobrun( huffw, actbl, eobrun );
			encode_crbits( huffw, storw );		
		}
	}	
	
	// return eob
	return eob;
}


/* -----------------------------------------------
	run of EOB SA decoding routine
	----------------------------------------------- */
int decode_eobrun_sa( abitreader* huffr, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	int bpos;
	
	
	// fast eobrun decoding routine for succesive approximation
	for ( bpos = from; bpos <= to; bpos++ ) {
		if ( block[ bpos ] != 0 ) {
			n = huffr->read( 1 );
			block[ bpos ] = ( block[ bpos ] > 0 ) ? n : -n;
		}
	}
	
	// decrement eobrun
	(*eobrun)--;
	
	
	return 0;
}


/* -----------------------------------------------
	run of EOB encoding routine
	----------------------------------------------- */
int encode_eobrun( abitwriter* huffw, huffCodes* actbl, int* eobrun )
{
	unsigned short n;
	unsigned char  s;
	int hc;
	
	
	if ( (*eobrun) > 0 ) {
		while ( (*eobrun) > actbl->max_eobrun ) {
			huffw->write( actbl->cval[ 0xE0 ], actbl->clen[ 0xE0 ] );
			huffw->write( E_ENVLI( 14, 32767 ), 14 );
			(*eobrun) -= actbl->max_eobrun;
		}
		BITLEN( s, (*eobrun) );
		s--;
		n = E_ENVLI( s, (*eobrun) );
		hc = ( s << 4 );
		huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
		huffw->write( n, s );
		(*eobrun) = 0;
	}

	
	return 0;
}


/* -----------------------------------------------
	correction bits encoding routine
	----------------------------------------------- */
int encode_crbits( abitwriter* huffw, abytewriter* storw )
{	
	unsigned char* data;
	int len;
	int i;
	
	
	// peek into data from abytewriter	
	len = storw->getpos();
	if ( len == 0 ) return 0;
	data = storw->peekptr();
	
	// write bits to huffwriter
	for ( i = 0; i < len; i++ )
		huffw->write( data[ i ], 1 );
	
	// reset abytewriter, discard data
	storw->reset();
	
	
	return 0;
}


/* -----------------------------------------------
	returns next code (from huffman-tree & -data)
	----------------------------------------------- */
int next_huffcode( abitreader *huffw, huffTree *ctree )
{	
	int node = 0;
	
	
	while ( node < 256 ) {
		node = ( huffw->read( 1 ) == 1 ) ?
				ctree->r[ node ] : ctree->l[ node ];
		if ( node == 0 ) break;
	}
	
	return ( node - 256 );
}


/* -----------------------------------------------
	calculates next position for MCU
	----------------------------------------------- */
int next_mcupos( int* mcu, int* cmp, int* csc, int* sub, int* dpos, int* rstw )
{
	int sta = 0; // status
	
	
	// increment all counts where needed
	if ( ( ++(*sub) ) >= cmpnfo[(*cmp)].mbs ) {
		(*sub) = 0;
		
		if ( ( ++(*csc) ) >= cs_cmpc ) {
			(*csc) = 0;
			(*cmp) = cs_cmp[ 0 ];
			(*mcu)++;
			if ( (*mcu) >= mcuc ) sta = 2;
			else if ( rsti > 0 )
				if ( --(*rstw) == 0 ) sta = 1;
		}
		else {
			(*cmp) = cs_cmp[(*csc)];
		}
	}
	
	// get correct position in image ( x & y )
	if ( cmpnfo[(*cmp)].sfh > 1 ) { // to fix mcu order
		(*dpos)  = ( (*mcu) / mcuh ) * cmpnfo[(*cmp)].sfh + ( (*sub) / cmpnfo[(*cmp)].sfv );
		(*dpos) *= cmpnfo[(*cmp)].bch;
		(*dpos) += ( (*mcu) % mcuh ) * cmpnfo[(*cmp)].sfv + ( (*sub) % cmpnfo[(*cmp)].sfv );
	}
	else if ( cmpnfo[(*cmp)].sfv > 1 ) {
		// simple calculation to speed up things if simple fixing is enough
		(*dpos) = ( (*mcu) * cmpnfo[(*cmp)].mbs ) + (*sub);
	}
	else {
		// no calculations needed without subsampling
		(*dpos) = (*mcu);
	}
	
	
	return sta;
}


/* -----------------------------------------------
	calculates next position (non interleaved)
	----------------------------------------------- */
int next_mcuposn( int* cmp, int* dpos, int* rstw )
{
	// increment position
	(*dpos)++;
	
	// fix for non interleaved mcu - horizontal
	if ( cmpnfo[(*cmp)].bch != cmpnfo[(*cmp)].nch ) {
		if ( (*dpos) % cmpnfo[(*cmp)].bch == cmpnfo[(*cmp)].nch )
			(*dpos) += ( cmpnfo[(*cmp)].bch - cmpnfo[(*cmp)].nch );
	}
	
	// fix for non interleaved mcu - vertical
	if ( cmpnfo[(*cmp)].bcv != cmpnfo[(*cmp)].ncv ) {
		if ( (*dpos) / cmpnfo[(*cmp)].bch == cmpnfo[(*cmp)].ncv )
			(*dpos) = cmpnfo[(*cmp)].bc;
	}
	
	// check position
	if ( (*dpos) >= cmpnfo[(*cmp)].bc ) return 2;
	else if ( rsti > 0 )
		if ( --(*rstw) == 0 ) return 1;
	

	return 0;
}


/* -----------------------------------------------
	skips the eobrun, calculates next position
	----------------------------------------------- */
int skip_eobrun( int* cmp, int* dpos, int* rstw, int* eobrun )
{
	if ( (*eobrun) > 0 ) // error check for eobrun
	{		
		// compare rst wait counter if needed
		if ( rsti > 0 ) {
			if ( (*eobrun) > (*rstw) )
				return -1;
			else
				(*rstw) -= (*eobrun);
		}
		
		// fix for non interleaved mcu - horizontal
		if ( cmpnfo[(*cmp)].bch != cmpnfo[(*cmp)].nch ) {
			(*dpos) += ( ( ( (*dpos) % cmpnfo[(*cmp)].bch ) + (*eobrun) ) /
						cmpnfo[(*cmp)].nch ) * ( cmpnfo[(*cmp)].bch - cmpnfo[(*cmp)].nch );
		}
		
		// fix for non interleaved mcu - vertical
		if ( cmpnfo[(*cmp)].bcv != cmpnfo[(*cmp)].ncv ) {
			if ( (*dpos) / cmpnfo[(*cmp)].bch >= cmpnfo[(*cmp)].ncv )
				(*dpos) += ( cmpnfo[(*cmp)].bcv - cmpnfo[(*cmp)].ncv ) *
						cmpnfo[(*cmp)].bch;
		}		
		
		// skip blocks 
		(*dpos) += (*eobrun);
		
		// reset eobrun
		(*eobrun) = 0;
		
		// check position
		if ( (*dpos) == cmpnfo[(*cmp)].bc ) return 2;
		else if ( (*dpos) > cmpnfo[(*cmp)].bc ) return -1;
		else if ( rsti > 0 ) 
			if ( (*rstw) == 0 ) return 1;
	}
	
	return 0;
}


/* -----------------------------------------------
	creates huffman-codes & -trees from dht-data
	----------------------------------------------- */
void build_huffcodes( unsigned char *clen, unsigned char *cval,	huffCodes *hc, huffTree *ht )
{
	int nextfree;	
	int code;
	int node;
	int i, j, k;
	
	
	// fill with zeroes
	memset( hc->clen, 0, 256 * sizeof( short ) );
	memset( hc->cval, 0, 256 * sizeof( short ) );
	memset( ht->l, 0, 256 * sizeof( short ) );
	memset( ht->r, 0, 256 * sizeof( short ) );
	
	// 1st part -> build huffman codes
	
	// creating huffman-codes	
	k = 0;
	code = 0;	
	
	// symbol-value of code is its position in the table
	for( i = 0; i < 16; i++ ) {
		for( j = 0; j < (int) clen[ i ]; j++ ) {
			hc->clen[ (int) cval[k] ] = 1 + i;
			hc->cval[ (int) cval[k] ] = code;
			
			k++;			
			code++;
		}		
		code = code << 1;
	}
	
	// find out eobrun max value
	hc->max_eobrun = 0;
	for ( i = 14; i >= 0; i-- ) {
		if ( hc->clen[ i << 4 ] > 0 ) {
			hc->max_eobrun = ( 2 << i ) - 1;
			break;
		}
	}
	
	// 2nd -> part use codes to build the coding tree
	
	// initial value for next free place
	nextfree = 1;

	// work through every code creating links between the nodes (represented through ints)
	for ( i = 0; i < 256; i++ )	{
		// (re)set current node
		node = 0;   		   		
		// go through each code & store path
		for ( j = hc->clen[ i ] - 1; j > 0; j-- ) {
			if ( BITN( hc->cval[ i ], j ) == 1 ) {
				if ( ht->r[ node ] == 0 )
					 ht->r[ node ] = nextfree++;
				node = ht->r[ node ];
			}
			else{
				if ( ht->l[ node ] == 0 )
					ht->l[ node ] = nextfree++;
				node = ht->l[ node ];
			}   					
		}
		// last link is number of targetvalue + 256
		if ( hc->clen[ i ] > 0 ) {
			if ( BITN( hc->cval[ i ], 0 ) == 1 )
				ht->r[ node ] = i + 256;
			else
				ht->l[ node ] = i + 256;
		}	   	
	}
}

/* ----------------------- End of JPEG specific functions -------------------------- */

/* ----------------------- Begin of JPEG modifying functions -------------------------- */


/* -----------------------------------------------
	JFIF header rebuilding routine
	----------------------------------------------- */
bool rebuild_header( void )
{	
	abytewriter* hdrw; // new header writer
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header	
	
	
	// start headerwriter
	hdrw = new abytewriter( 4096 );
	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		// discard any unneeded meta info
		if ( ( type == 0xDA ) || ( type == 0xC4 ) || ( type == 0xDB ) ||
			 ( type == 0xC0 ) || ( type == 0xC1 ) || ( type == 0xC2 ) ||
			 ( type == 0xDD ) ) {
			hdrw->write_n( &(hdrdata[ hpos ]), len );
		}
		hpos += len;
	}
	
	// replace current header with the new one
	free( hdrdata );
	hdrdata = hdrw->getptr();
	hdrs    = hdrw->getpos();
	delete( hdrw );
	
	
	return true;
}


/* -----------------------------------------------
	RST marker removing routine
	----------------------------------------------- */
bool remove_rst( void )
{
	abytewriter* hdrw; // new header writer
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header	
	
	
	// start headerwriter
	hdrw = new abytewriter( 4096 );
	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		// discard any unneeded meta info
		if ( type != 0xDD ) {
			hdrw->write_n( &(hdrdata[ hpos ]), len );
		}
		hpos += len;
	}
	
	// replace current header with the new one
	free( hdrdata );
	hdrdata = hdrw->getptr();
	hdrs    = hdrw->getpos();
	delete( hdrw );
	
	// reset restart interval data
	if ( rstp != NULL ) free ( rstp );
	rstp = NULL;
	rsti = 0;
	
	
	return true;
}


/* -----------------------------------------------
	recalculate table by division factor
	----------------------------------------------- */
bool recalc_qtable_div( float df, int tno )
{
	float qst_f;
	int qst_i;
	int i;
	
	
	// recalculate table - multiply by division factor
	for ( i = 0; i < 64; i++ ) {
		qst_f = qtables[ tno ][ i ];
		qst_f = qst_f * df;
		qst_i = ROUND_F( qst_f );
		qst_i = CLAMPED( 0, 65535, qst_i );
		qtables[ tno ][ i ] = qst_i;
	}
	
	
	return true;
}


/* -----------------------------------------------
	put qtables back in place
	----------------------------------------------- */
bool reinsert_qtables( void )
{
	unsigned char type = 0x00; // type of current marker segment
	unsigned int  len  = 0; // length of current marker segment
	unsigned int  hpos = 0; // position in header	
	unsigned int  fpos = 0; // end of marker position
	
	int cmp, dpos;
	float coef;
	float df;
	int qst0;
	int qst1;
	int p, n;
	int i;
	
	
	// seek for quantization tables in header
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		if ( type == 0xDB ) { // DQT segment / fix up quantization steps
			fpos = hpos + len; // reassign length to end position
			hpos += 4; // skip marker & length
			while ( hpos < fpos ) {
				p = LBITS( hdrdata[ hpos ], 4 ); // precision (0->8bit / 1->16bit)
				n = RBITS( hdrdata[ hpos ], 4 );
				hpos++;
				// table found
				if ( p == 1 ) { // 16 bit precision
					for ( i = 0; i < 64; i++ ) {
						// fetch steps old & new
						qst0 = B_SHORT( hdrdata[ hpos + (2*i) ], hdrdata[ hpos + (2*i) + 1 ] );
						qst1 = qtables[ n ][ i ];
						// calc division factor
						df  = ( qst1 != 0 ) ? ( float ) qst0 / ( float ) qst1 : 0;
						if ( df != 1 ) {
							// fix coefficients
							for ( cmp = 0; cmp < cmpc; cmp++ ) {
								if ( cmpnfo[ cmp ].qtblno == n ) {
									for ( dpos = 0; dpos < cmpnfo[ cmp ].bc; dpos++ ) {
										coef = colldata[ cmp ][ i ][ dpos ];
										coef = coef * df;
										colldata[ cmp ][ i ][ dpos ] = ROUND_F( coef );
									}
								}
							}
							// reinsert into header
							hdrdata[ hpos + (2*i) ]     = LBITS16( qst1, 8 );
							hdrdata[ hpos + (2*i) + 1 ] = RBITS16( qst1, 8 );
						}
						// fix up internal table
						// qtables[ n ][ i ] = qst1;
					}
					hpos += 128;
				}
				else { // for 8 bit precision
					for ( i = 0; i < 64; i++ ) {
						// fetch steps old & new
						qst0 = hdrdata[ hpos + i ];
						qst1 = ( qtables[ n ][ i ] < 256 ) ? qtables[ n ][ i ] : 255;
						// calc division factor
						df  = ( qst1 != 0 ) ? ( float ) qst0 / ( float ) qst1 : 0;
						if ( df != 1 ) {
							// fix coefficients
							for ( cmp = 0; cmp < cmpc; cmp++ ) {
								if ( cmpnfo[ cmp ].qtblno == n ) {
									for ( dpos = 0; dpos < cmpnfo[ cmp ].bc; dpos++ ) {
										coef = colldata[ cmp ][ i ][ dpos ];
										coef = coef * df;
										colldata[ cmp ][ i ][ dpos ] = ROUND_F( coef );
									}
								}
							}
							// reinsert into header
							hdrdata[ hpos + i ] = qst1;
						}
						// fix up internal table
						qtables[ n ][ i ] = qst1;
					}
					hpos += 64;
				}
			}			
		}
		else { // skip segment
			hpos += len;
		}
	}
	
	
	return true;
}


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
	calculate current quality for table #
	----------------------------------------------- */
float calc_quality( int tno )
{
	float s = 0;
	float q = 0;
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
	if ( ref == -1 )
		return -1;
	
	// compare tables, calculate average scaling factor s
	for ( i = 0; i < 64; i++ ) {
		s += ( ( float ) qtables[ tno ][ i ] ) / ( ( float ) std_qtables[ ref ][ jpeg_natural_order[i] ] );
	}
	s /= 64;
	
	// calculate quality based on s
	q = ( s < 1.0 ) ? 1.0 - (s/2) : 1.0 / (2*s);
	
	
	// return q
	return q;
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
		s = ( ( float ) qtables[ tno ][ i ] ) / ( ( float ) std_qtables[ ref ][ jpeg_natural_order[i] ] );
		sl = ( s > sl ) ? s : sl;
		sh = ( s < sh ) ? s : sh;
	}
	
	// calculate qualities based on s
	(*low)  = ( sl < 1.0 ) ? 1.0 - (sl/2) : 1.0 / (2*sl);
	(*high) = ( sh < 1.0 ) ? 1.0 - (sh/2) : 1.0 / (2*sh);
}


/* ----------------------- End of JPEG modifying functions -------------------------- */

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
	calculates median out of an integer array
	----------------------------------------------- */
int median_int( int* values, int size )
{
	int middle = ( size >> 1 );
	bool done;
	int swap;
	int i;
	
	
	// sort data first
	done = false;
	while ( !done ) {
		done = true;
		for ( i = 1; i < size; i++ )
		if ( values[ i ] < values[ i - 1 ] ) {
			swap = values[ i ];
			values[ i ] = values[ i - 1 ];
			values[ i - 1 ] = swap;
			done = false;
		}
	}
	
	// return median
	return ( ( size % 2 ) == 0 ) ?
		( values[ middle ] + values[ middle - 1 ] ) / 2 : values[ middle ];
}

/* -----------------------------------------------
	calculates median out of an float array
	----------------------------------------------- */
float median_float( float* values, int size )
{
	int middle = ( size >> 1 );
	bool done;
	float swap;
	int i;
	
	
	// sort data first
	done = false;
	while ( !done ) {
		done = true;
		for ( i = 1; i < size; i++ )
		if ( values[ i ] < values[ i - 1 ] ) {
			swap = values[ i ];
			values[ i ] = values[ i - 1 ];
			values[ i - 1 ] = swap;
			done = false;
		}
	}
	
	// return median	
	if ( ( size % 2 ) == 0 ) {
		return ( values[ middle ] + values[ middle - 1 ] ) / 2.0;
	}
	else
		return ( values[ middle ] );
}

/* ----------------------- End of miscellaneous helper functions -------------------------- */

/* ----------------------- End of file -------------------------- */
