#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <math.h>
#include <time.h>
#include <sys\stat.h>
#include <string>
#include <map>

#define MAX( x,y ) ( x < y ) ? y : x
#define MIN( x,y ) ( x > y ) ? y : x
#define ABSDIFF( x,y ) ( x > y ) ? ( x - y ) : ( y - x )
#define R( x ) h[0] - h[x]
#define FSIZTEST ( totalbytes % bpv ) ? " ( ! )" : ""
#define DSYMTEST ( kpv < 0.5 ) ? "" : " ( ! )"

using namespace std;

typedef map<unsigned int, unsigned int> freqtable;
typedef map<unsigned int, unsigned int>::iterator freqindex;


/* -----------------------------------------------
	declarations - functions ( main )
	----------------------------------------------- */

bool treatfile( void );
bool fillfreq( char* fbuff );
void calcmodel( void );
void calcsimples( void );
void calcfinals( void );
void calcmaxdiff( char* fbuff );

signed int convert_uval( unsigned int uval );
unsigned int get_uval( char* cval, int bsize );
unsigned int get_bytes( unsigned int uval, bool sgn, int pos, int bts );
char* extension( char* filename );

bool initcsv( void );
bool add2csv( void );
bool finalizecsv( void );

bool checkfile( char* filename );
int getfilesize( char* filename );

void showhelp( void );
void showresults( void );
void showfinalresults( void );
void showerrormessage( void );
void printcmb( FILE* outf, unsigned int uval, bool sgn );

	
/* -----------------------------------------------
	declarations - global variables
	----------------------------------------------- */

static char* headline = "--- hcheck v2.0 by Matthias Stirner ---";

static bool wait = false;
char errormessage[256];

char* filename;
double range; // highest possible value ( depends on bpv )
unsigned int bpv = 1; // bytes per value ( 0 < bpv <= 4 )
unsigned int ordmax = 1; // maximum order of model ( 0 < ordmax <= 4 )
unsigned int order; // current order of model
unsigned int filecount; // counted files

freqtable freqtbl; // frequency of values
unsigned int totalvalues; // total read values
unsigned int totalbytes; // total read bytes
double h[5];  // conditional entropy -> h[0] -> decision content; h[1] -> entropy 1st order; ...
double hc[5]; // combined entropy ( needed for calculations )
double mu; // average value ( unsigned )
double ms; // average value ( signed )
double pcomp; // possible compression %
unsigned int pcomb; // possible compression bytes
unsigned int k; // number of different symbols
unsigned int umin; // minimum value found ( unsigned )
unsigned int umax; // maximum value found ( unsigned )
unsigned int smin; // minimum value found ( signed )
unsigned int smax; // maximum value found ( signed )
unsigned int urng; // range ( unsigned )
unsigned int srng; // range ( signed ) 
unsigned int pmax; // value with maximum probability
unsigned int umdiff; // maximum difference ( unsigned )
unsigned int smdiff; // maximum difference ( signed )
unsigned int usdiff; // average difference ( unsigned )
unsigned int ssdiff; // average difference ( signed )
double kpv; // different symbols per symbol

unsigned int gumin = ( unsigned int ) -1; // global umin
unsigned int gumax = 0; // global umax
signed int gsmin = ( 1 << 31 ) - 1 ; // global smin
signed int gsmax = ( 1 << 31 ) * ( -1 ); // global smax
unsigned int urngmax = 0; // maximum range ( unsigned )
unsigned int srngmax = 0; // maximum range ( signed )
unsigned int gumdiff = 0; // global maximum difference ( unsigned )
unsigned int gsmdiff = 0; // global maximum difference ( signed )
unsigned int kmax = 0; // maximum different symbols

double hmin = 32;   char* hminfn; // minimum entropy
double hmax = 0;    char* hmaxfn; // maximum entropy
double h0min = 32;  char* h0minfn; // minimum content
double h0max = 0;   char* h0maxfn; // maximum content

double pcmin = 0;   char* pcminfn; // minimum possible compression
double pcmax = 100; char* pcmaxfn; // maximum possible compression
double pcsum = 0;   // sum of possible compression
unsigned int tbsum = 0; // sum of total bytes
unsigned int tvsum = 0; // sum of total values
unsigned int pbsum[5] = { 0,0,0,0,0 }; // sum of possible compressed size ( bytes )

char* csvfilename = "hcheck.csv";
char* csvsep = ";";
bool wcsv = false;
bool quiet = false;


/* -----------------------------------------------
	Reads in arguments, calls treatfile()
	----------------------------------------------- */
		
int main( int argc, char** argv )
{
	int argpos = 0;
	filecount = 0;
	char trash;
	
	printf( "\n%s\n", headline );
	
	
	if ( argc < 2 ) {
		showhelp();
		exit( 0 );
	}
	
	for ( argpos = 1; argpos < argc; argpos++ ) {
		if ( argv[argpos][0] == '-' ) {
   		if ( ( argv[argpos][1] == 'b' ) ||
   		    ( argv[argpos][1] == 'u' ) ||
   		    ( argv[argpos][1] == 's' ) ) {
	   		bpv = strtoul( argv[argpos] + 2, NULL, 10 );
   		}
   		else if ( argv[argpos][1] == 'o' ) {
	   		ordmax = strtoul( argv[argpos] + 2, NULL, 10 );
	   		ordmax++;
   		}
   		else if ( strcmp( argv[argpos] + 1, "csv" ) == 0 ) {
   			wcsv = true;
   		}
   		else if ( strcmp( argv[argpos] + 1, "quiet" ) == 0 ) {
   			quiet = true;
   		}
   		else{
	   		showhelp();
	   		exit( 0 );
	   	}
		}
	}
	
	if ( ( ordmax < 1 ) || ( ordmax > 4 ) ||
		( bpv < 1 ) || ( bpv > 4 ) ||
		( argv[argc - 1][0] == '-' ) ) {
		showhelp();
		return 0;
	}
	
	while ( ( ordmax * bpv ) > 4 ) ordmax--;
	range = pow( 2, bpv * 8 );
	
	
	if ( wcsv ) initcsv();
	
	for ( argpos = 1; argpos < argc; argpos++ ) {
		if ( argv[argpos][0] != '-' ) {
			filename = argv[argpos];
			if ( treatfile() ) filecount++;
		}
	}	
	
	if ( filecount > 0 ) showfinalresults();
	if ( wcsv ) ( filecount > 0 ) ? finalizecsv() : remove( csvfilename );
	
	if ( ( wait ) && ( !quiet ) ) {
		printf( "\n<please press any key>" );
		scanf( "%c", &trash );
	}
	
	
	return 1;
}


/* -----------------------------------------------
	Treats one file - puts out result, if ok
	----------------------------------------------- */
	
bool treatfile( void )
{
	FILE* file;
	char* fbuff;
	
	
	if ( checkfile( filename ) ) { // get file into buffer
		fbuff = new char[totalbytes + 3]; // 3 extra for combinations
		file = fopen( filename, "rb" );
		fread( fbuff, 1, totalbytes, file );
		fclose( file );
		fbuff[totalbytes + 0] = fbuff[0];
		fbuff[totalbytes + 1] = fbuff[1];
		fbuff[totalbytes + 2] = fbuff[2];
	}
	else{
		fbuff = new char[1];
		if ( !quiet ) showerrormessage();
		return false;		
	}
	
	order = 1;	
	calcmaxdiff( fbuff );
	while ( order <= ordmax ) {
		if ( !fillfreq( fbuff ) ) {
			if ( !quiet ) showerrormessage();
			freqtbl.clear();
			return false;
		}
		calcmodel();
		// simples can only be calculated for 1st order model
		if ( order == 1 ) calcsimples();
		
		if ( wcsv ) add2csv();
		if ( !quiet ) showresults();				
		order++;
	}
	
	calcfinals();
	freqtbl.clear();
	
	return true;
}


/* -----------------------------------------------
	Fills freqtable from filebuffer
	----------------------------------------------- */
	
bool fillfreq( char* fbuff )
{	
	unsigned int val;
	unsigned int i;
	
	if ( ( order * bpv ) > 4 ) {
		sprintf( errormessage, "cannot calculate %u. order model for %ubit values", order - 1, bpv*8 );
		return false;
	}	
	
	totalvalues = 0;
	k = 0;
	
	// clear freqtable
	freqtbl.clear();
	
	// fill freqtable
	for ( i = 0; i < totalbytes; i += bpv ) {
		val = get_uval( fbuff + i, order * bpv );
		totalvalues++;
		if ( freqtbl.find( val ) == freqtbl.end() ) {
			freqtbl[val] = 1;
			k++;
		}
		else freqtbl[val]++;
	}	
	
	// calculate different symbols per symbol
	kpv = ( ( double ) k ) / ( ( double ) totalvalues );
	
	if ( k > freqtbl.max_size() ) {
		sprintf( errormessage, "too many values from file \"%s\"", filename );
		return false;
	}
	if ( freqtbl.empty() ) {
		sprintf( errormessage, "couldn't read any values from file \"%s\"", filename );
		return false;
	}
	
	return true;
}


/* -----------------------------------------------
	Calculates everything of interest ( h,h0,r,... )
	----------------------------------------------- */
	
void calcmodel( void )
{
	freqindex fi = freqtbl.begin();	
	double prob;
	
	hc[order] = 0;
	pmax = ( *fi ).first;
		
	// maximum probability value/combination & combined entropy
	while ( fi != freqtbl.end() ) {
		prob = ( ( double ) ( *fi ).second ) / ( ( double ) totalvalues ); // probability ( of one value )
		pmax = ( freqtbl[pmax] > ( *fi ).second ) ? pmax : ( *fi ).first; // max prob value
		hc[order] -= ( prob * log2( prob ) ); // entropy
		fi++;		
	}
	
	// conditional entropy
	if ( order > 1 ) h[order] = hc[order] - hc[order - 1];
	else h[order] = hc[order];
	
	// possible compression ( % & bytes )
	pcomp = ( 12.5 * ( h[order] / ( ( double ) bpv ) ) );
	pcomb = ( unsigned long int ) ( pcomp * totalbytes / 100 );
	pbsum[order] += pcomb;
}


/* -----------------------------------------------
	Calculates simple things ( min/max/m )
	----------------------------------------------- */
void calcsimples( void )
{
	freqindex fi;
	
	mu = 0;
	ms = 0;
	umin = 0;
	umax = ( unsigned int ) ( range - 1 );
	smin = ( unsigned int ) ( range / 2 );
	smax = ( unsigned int ) ( range / 2 ) - 1;
	usdiff = 0;
	ssdiff = 0;
	
	// minimum value - unsigned
	fi = freqtbl.begin();
	umin = ( *fi ).first;
	
	// maximum value - unsigned
	fi = freqtbl.end(); fi--;
	umax = ( *fi ).first;
	
	// minimum value - signed
	if ( ( smin <= umin ) || ( smin > umax ) ) smin = umin; // signed smin < 0
	else{
		fi = freqtbl.end(); fi--;
		while ( ( *fi ).first >= smin ) fi--;	
		smin = ( *++fi ).first;				
	}
	
	// maximum value - signed
	if ( ( smax >= umax ) || ( smax < umin ) ) smax = umax; // signed smax > 0
	else{
		fi = freqtbl.begin();
		while ( ( *fi ).first <= smax ) fi++;
		smax = ( *--fi ).first;
	}
	
	// ranges of values ( signed/unsigned )
	urng = umax - umin;
	srng = convert_uval( smax ) - convert_uval( smin );
	
	// average values ( uval & sval )
	fi = freqtbl.begin();
	while ( fi != freqtbl.end() ) {
		mu += ( ( double ) ( *fi ).first ) * ( ( double ) ( *fi ).second );
		ms += ( ( double ) convert_uval( ( *fi ).first ) ) * ( ( double ) ( *fi ).second );
		fi++;
	}
	mu /= totalvalues;
	ms /= totalvalues;
	
	// standard difference
	/* IMPLEMENT
	fi = freqtbl.begin();
	while ( fi != freqtbl.end() ) {
		usdiff += ( ( double ) ( *fi ).first ) * ( double ) ( ( *fi ).second - mu );
		mu += ( ( double ) ( *fi ).first ) * ( ( double ) ( *fi ).second );
		ms += ( ( double ) convert_uval( ( *fi ).first ) ) * ( ( double ) ( *fi ).second );
		fi++;
	}
	*/
	
	// decision content
	h[0] = log2( k );
	
	gumin = MIN( umin, gumin ); // global umin
	gumax = MAX( umax, gumax ); // global umax
	gsmin = MIN( convert_uval( smin ), gsmin ); // global smin
	gsmax = MAX( convert_uval( smax ), gsmax ); // global smax
	urngmax = MAX( urngmax, urng ); // maximum range ( unsigned )
	srngmax = MAX( srngmax, srng ); // maximum range ( signed )
}


/* -----------------------------------------------
	Calculates things for final results
	----------------------------------------------- */
	
void calcfinals( void )
{
	// calculate kmax
	kmax = MAX( k, kmax );
	
	// increment tvsum
	tvsum += totalvalues;
		
	if ( h[ordmax] < hmin ) {
		hmin = h[ordmax];
		hminfn = new char[strlen( filename ) + 1];
		strcpy( hminfn, filename );
	}
	if ( h[ordmax] > hmax ) {
		hmax = h[ordmax];
		hmaxfn = new char[strlen( filename ) + 1];
		strcpy( hmaxfn, filename );
	}
	
	if ( h[0] < h0min ) {
		h0min = h[0];
		h0minfn = new char[strlen( filename ) + 1];
		strcpy( h0minfn, filename );
	}
	if ( h[0] > h0max ) {
		h0max = h[0];
		h0maxfn = new char[strlen( filename ) + 1];
		strcpy( h0maxfn, filename );
	}
	
	if ( pcomp > pcmin ) {
		pcmin = pcomp;
		pcminfn = new char[strlen( filename ) + 1];
		strcpy( pcminfn, filename );
	}
	if ( pcomp < pcmax ) {
		pcmax = pcomp;
		pcmaxfn = new char[strlen( filename ) + 1];
		strcpy( pcmaxfn, filename );
	}
	pcsum += pcomp;
}


/* -----------------------------------------------
	Calculates maximum difference between values
	----------------------------------------------- */
	
void calcmaxdiff( char* fbuff )
{	
	unsigned int middle = 1 << ( ( 8 * bpv ) - 1 );
	unsigned int uval[3];
	unsigned int sval[3];
	unsigned int i;
	
	uval[0] = get_uval( fbuff, bpv );
	sval[0] = uval[0] + ( ( uval[0] > middle ) ? -middle : middle );
	umdiff = 0;
	smdiff = 0;
	
	for ( i = bpv; i < totalbytes; i += bpv ) {
		uval[1] = get_uval( fbuff + i, bpv );
		sval[1] = uval[1] + ( ( uval[1] > middle ) ? -middle : middle );
		uval[2] = ABSDIFF( uval[0], uval[1] );
		sval[2] = ABSDIFF( sval[0], sval[1] );
		umdiff = MAX( umdiff, uval[2] );
		smdiff = MAX( smdiff, sval[2] );
		uval[0] = uval[1];
		sval[0] = sval[1];
	}
	
	gumdiff = MAX( gumdiff, umdiff ); // global maximum difference ( unsigned )
	gsmdiff = MAX( gsmdiff, smdiff ); // global maximum difference ( signed )
}


/* -----------------------------------------------
	Converts unsigned to signed
	----------------------------------------------- */
	
signed int convert_uval( unsigned int uval )
{
	signed int sval;
	
	if ( uval >= ( unsigned int ) ( range / 2 ) ) sval = ( signed int ) ( uval - range );
	else sval = uval;
	
	return sval;
}


/* -----------------------------------------------
	Gets uint from array of chars
	----------------------------------------------- */
	
unsigned int get_uval( char* cval, int bsize )
{
	unsigned int uval = 0;
	int i;
	
	for ( i = 0; i < (int) bsize; i++ ) {
		uval += ( ( ( unsigned char ) cval[i] ) << ( 8*i ) );
	}
	
	return uval;
}


/* -----------------------------------------------
	Returns bytes at position "pos"
	----------------------------------------------- */
	
unsigned int get_bytes( unsigned int uval, bool sgn, int pos, int bts )
{
	uval = uval >> ( 8 * pos * bts );
	uval = uval << ( 32 - ( 8 * bts ) );
	uval = uval >> ( 32 - ( 8 * bts ) );
	return uval;
}


/* -----------------------------------------------
	Returns pointer to extension
	----------------------------------------------- */
	
char* extension( char* filename )
{	
	int length = strlen( filename );
	int dotpos = length;
	int i;
	
	for ( i = length - 1; i > 0; i-- ) {
		if ( filename[i] == '.' ) {
			dotpos = i;
			break;
		}
	}
	
	return filename + dotpos + 1;
}


/* -----------------------------------------------
	Initilialize csv
	----------------------------------------------- */
	
bool initcsv( void )
{
	FILE *csvfile;
	time_t tsec;
	tm *tnow;
	char date[16];
	int i;
	
	time( &tsec );
	tnow = localtime( &tsec );
	sprintf( date, "%i-%i-%i", tnow->tm_mday, 1 + tnow->tm_mon, 1900 + tnow->tm_year );
	
	i = 0;
	csvfilename = new char[256];
	sprintf( csvfilename, "hcheck-%s-%i.csv", date, i );
	while ( access( csvfilename, F_OK ) == 0 ) {
		i++;
		sprintf( csvfilename, "hcheck-%s-%i.csv", date, i );
	}
		
	if ( !( csvfile = fopen( csvfilename, "w" ) ) ) return false;
	fprintf( csvfile, "filename%s", csvsep );
	fprintf( csvfile, "filetype%s", csvsep );
	fprintf( csvfile, "order%s", csvsep );
	fprintf( csvfile, "bytes%s", csvsep );
	fprintf( csvfile, "values%s", csvsep );
	fprintf( csvfile, "diff. values%s", csvsep );
	fprintf( csvfile, "bpv%s", csvsep );
	fprintf( csvfile, "umin%s", csvsep );
	fprintf( csvfile, "umax%s", csvsep );
	fprintf( csvfile, "smin%s", csvsep );
	fprintf( csvfile, "smax%s", csvsep );
	fprintf( csvfile, "range ( uval )%s", csvsep );
	fprintf( csvfile, "range ( sval )%s", csvsep );
	fprintf( csvfile, "max diff ( uval )%s", csvsep );
	fprintf( csvfile, "max diff ( sval )%s", csvsep );
	fprintf( csvfile, "max freq ( uval )%s", csvsep );
	fprintf( csvfile, "max freq ( sval )%s", csvsep );
	fprintf( csvfile, "max freq ( freq )%s", csvsep );
	fprintf( csvfile, "average uval%s", csvsep );
	fprintf( csvfile, "average sval%s", csvsep );	
	fprintf( csvfile, "decision content%s", csvsep );
	fprintf( csvfile, "cond. entropy%s", csvsep );
	fprintf( csvfile, "comb. entropy%s", csvsep );
	fprintf( csvfile, "redundancy%s", csvsep );
	fprintf( csvfile, "comp. size ( percent )%s", csvsep );
	fprintf( csvfile, "comp. size ( bytes )%s", csvsep );
	fprintf( csvfile, "\n" );
	fclose( csvfile );
	
	return true;
}


/* -----------------------------------------------
	Adds data to csvfile
	----------------------------------------------- */
	
bool add2csv( void )
{
	FILE* csvfile;
	
	if ( !( csvfile = fopen( csvfilename, "a" ) ) ) return false;
	fseek( csvfile, 0, SEEK_END );
	fprintf( csvfile, "%s%s", filename, csvsep );
	fprintf( csvfile, "%s%s", extension( filename ), csvsep );
	fprintf( csvfile, "%u%s", order - 1, csvsep );
	fprintf( csvfile, "%u%s", totalbytes, csvsep );
	fprintf( csvfile, "%u%s", totalvalues, csvsep );
	fprintf( csvfile, "%u%s", k, csvsep );
	fprintf( csvfile, "%u%s", bpv, csvsep );
	if ( order == 1 ) {
	fprintf( csvfile, "%u%s", umin, csvsep );
	fprintf( csvfile, "%u%s", umax, csvsep );
	fprintf( csvfile, "%i%s", convert_uval( smin ), csvsep );
	fprintf( csvfile, "%i%s", convert_uval( smax ), csvsep );
	fprintf( csvfile, "%u%s", urng, csvsep );
	fprintf( csvfile, "%u%s", srng, csvsep );
	fprintf( csvfile, "%u%s", umdiff, csvsep );
	fprintf( csvfile, "%u%s", smdiff, csvsep );
	fprintf( csvfile, "%u%s", pmax, csvsep );
	fprintf( csvfile, "%i%s", convert_uval( pmax ), csvsep );
	fprintf( csvfile, "%u%s", freqtbl[pmax], csvsep );
	fprintf( csvfile, "%f%s", mu, csvsep );
	fprintf( csvfile, "%f%s", ms, csvsep );
	}
	else{
	fprintf( csvfile, "%s%s%s%s%s%s%s%s", csvsep, csvsep, csvsep, csvsep,
	                                     csvsep, csvsep, csvsep, csvsep );
	printcmb( csvfile, pmax, false );
	fprintf( csvfile, "%s", csvsep );
	printcmb( csvfile, pmax, true );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%u%s", freqtbl[pmax], csvsep );
	fprintf( csvfile, "%s%s", csvsep, csvsep );
	}
	fprintf( csvfile, "%f%s", h[0], csvsep );
	fprintf( csvfile, "%f%s", h[order], csvsep );
	fprintf( csvfile, "%f%s", hc[order], csvsep );
	fprintf( csvfile, "%f%s", R( order ), csvsep );
	fprintf( csvfile, "%f%s", pcomp, csvsep );
	fprintf( csvfile, "%u%s", pcomb, csvsep );
	fprintf( csvfile, "\n" );
	fclose( csvfile );
	
	return true;
}


/* -----------------------------------------------
	Adds final results to csv
	----------------------------------------------- */
	
bool finalizecsv( void )
{
	FILE* csvfile;
	int i;
	
	if ( !( csvfile = fopen( csvfilename, "a" ) ) ) return false;
	fseek( csvfile, 0, SEEK_END );
	fprintf( csvfile, "\n" );
	for ( i = 1; i <= (int) ordmax; i++ ) {
	fprintf( csvfile, "all files ( *.* )%s", csvsep );
	fprintf( csvfile, "*%s", csvsep );
	fprintf( csvfile, "%i%s", i - 1, csvsep );
	fprintf( csvfile, "%u%s", tbsum, csvsep );
	fprintf( csvfile, "%u%s", tvsum, csvsep );
	fprintf( csvfile, "%u%s", kmax, csvsep );
	fprintf( csvfile, "%i%s", bpv, csvsep );
	fprintf( csvfile, "%u%s", gumin, csvsep );
	fprintf( csvfile, "%u%s", gumax, csvsep );
	fprintf( csvfile, "%i%s", gsmin, csvsep );
	fprintf( csvfile, "%i%s", gsmax, csvsep );
	fprintf( csvfile, "%u%s", urngmax, csvsep );
	fprintf( csvfile, "%u%s", srngmax, csvsep );
	fprintf( csvfile, "%u%s", gumdiff, csvsep );
	fprintf( csvfile, "%u%s", gsmdiff, csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );	
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%s", csvsep );
	fprintf( csvfile, "%f%s", ( ( double ) pbsum[i] ) / ( ( double ) tbsum ) * 100.0, csvsep );
	fprintf( csvfile, "%u%s", pbsum[i], csvsep );
	fprintf( csvfile, "\n" );
	}
	return true;
}


/* -----------------------------------------------
	checks if file exists
	----------------------------------------------- */
	
bool checkfile( char* filename )
{	
	if ( access( filename, F_OK ) == -1 ) {
		sprintf( errormessage, "file \"%s\" doesn't exist", filename );
		return false;
	}
	
	if ( access( filename, R_OK ) == -1 ) {
		sprintf( errormessage, "file \"%s\" is not readable", filename );
		return false;
	}
	
	totalbytes = getfilesize( filename );
	tbsum += totalbytes;
	
	if ( totalbytes < bpv ) {
		sprintf( errormessage, "file \"%s\" contains no values", filename );
		return false;
	}
	
	return true;
}


/* -----------------------------------------------
	returns size of "filename"
	----------------------------------------------- */
	
int getfilesize( char* filename )
{
	struct stat sta;
	stat( filename, &sta );
	
	return sta.st_size;
}


/* -----------------------------------------------
	Shows help in case of wrong input
	----------------------------------------------- */
	
void showhelp( void )
{	
	printf( "\n" );
	printf( "Usage: hcheck [-b?] [-csv] [filename]\n" );
	printf( "\n" );
	printf( "    [-b?] -> bytes per value: must be > 0 and <= 4 ( optional )\n" );
	printf( "    [-o?] -> max order of entropy: must be > 0 and <= 4 ( optional )\n" );
	printf( "   [-csv] -> write csv-data to file \"%s\" ( optional )\n", csvfilename );
	printf( " [-quiet] -> reduce output to neccesary ( optional )\n" );
	printf( "\n" );
	printf( "Examples: \"hcheck -b2 test.jpg\"\n" );
	printf( "          \"hcheck -csv *.*\"\n" );	
}


/* -----------------------------------------------
	Shows results for one file
	----------------------------------------------- */
	
void showresults( void )
{
	printf( "\n" );
	printf( "--- file \"%s\" - %i. order model ---\n", filename, order - 1 );
	printf( "\n" );
	printf( "total bytes              :  %u%s byte\n", totalbytes, FSIZTEST );
	printf( "bytes per value          :  %u%s byte\n", bpv, FSIZTEST );
	if ( order == 1 ) {		
		printf( "total read values        :  %u values\n", totalvalues );
		printf( "no of different values   :  %u values\n", k );
		printf( "min value found          :  %u/%i ( uval/sval )\n", umin, convert_uval( smin ) );
		printf( "max value found          :  %u/%i ( uval/sval )\n", umax, convert_uval( smax ) );
		printf( "range of values          :  %u/%u ( uval/sval )\n", urng, srng );
		printf( "max difference found     :  %u/%i ( uval/sval )\n", umdiff, smdiff );
		printf( "max frequency found      :  %u/%i/%u ( uval/sval/freq )\n", pmax, convert_uval( pmax ), freqtbl[pmax] );
		printf( "average value ( unsigned ) :  %.2f\n", mu );
		printf( "average value ( signed )   :  %.2f\n", ms );
		printf( "decision content         :  %.3f bit\n", h[0] );
		printf( "entropy                  :  %.3f bit\n", h[order] );
		printf( "redundancy               :  %.3f bit\n", R( order ) );
	}
	else{
		printf( "total read combinations  :  %u combinations\n", totalvalues );
		printf( "no of diff. combinations :  %u combinations\n", k );	
		printf( "max frequency ( unsigned ) :  " ); printcmb( stdout, pmax, false ); printf( "\n" );
		printf( "max frequency ( signed )   :  " ); printcmb( stdout, pmax, true ); printf( "\n" );
		printf( "max frequency ( freq )     :  %u\n", freqtbl[pmax] );
		printf( "decision content         :  %.3f bit\n", h[0] );
		printf( "conditional entropy      :  %.3f bit\n", h[order] );
		printf( "combination entropy      :  %.3f bit\n", hc[order] );
		printf( "redundancy               :  %.3f bit\n", R( order ) );
	}
	printf( "possible comp. ratio     :  %.2f%s percent\n", pcomp, DSYMTEST );
	printf( "possible comp. size      :  %u%s byte\n", pcomb, DSYMTEST );
	printf( "\n" );
}


/* -----------------------------------------------
	Shows final results after everything is done
	----------------------------------------------- */
	
void showfinalresults( void )
{
	int i;
	
	printf( "\n\n" );
	printf( "--> final results ( %i file( s ) ) ( %u. order model ) <--\n", filecount, ordmax - 1 );
	printf( "\n" );
	printf( "total bytes read         :  %u byte\n", tbsum );
	for ( i = 1; i <= (int) ordmax; i++ ) {
	printf( "possible comp. size      :  %u byte ( %i. order )\n", pbsum[i], i - 1 );
	}
	printf( "\n" );
	printf( "min decision content     :  %.3f bit ( %s )\n", h0min, h0minfn );
	printf( "max decision content     :  %.3f bit ( %s )\n", h0max, h0maxfn );
	printf( "minimum cond. entropy    :  %.3f bit ( %s )\n", hmin, hminfn );
	printf( "maximum cond. entropy    :  %.3f bit ( %s )\n", hmax, hmaxfn );
	printf( "min possible comp. size  :  %.2f percent ( %s )\n", pcmin, pcminfn );
	printf( "max possible comp. size  :  %.2f percent ( %s )\n", pcmax, pcmaxfn );
	printf( "average pos. comp. size  :  %.2f percent\n\n", ( pcsum / filecount ) );
	if ( wcsv ) {
	printf( "results written to file  :  %s\n\n", csvfilename );
	}
	printf( "\n" );
}


/* -----------------------------------------------
	Shows message if something goes wrong
	----------------------------------------------- */

void showerrormessage( void )
{
	printf( "\n" );
	printf( "--- error while treating file \"%s\" ---\n", filename );
	printf( "\n" );
	printf( "%s", errormessage );
	printf( "\n" );
}


/* -----------------------------------------------
	Prints combination of values
	----------------------------------------------- */

void printcmb( FILE* outf, unsigned int uval, bool sgn )
{
	int i;
	
	if ( sgn )	fprintf( outf, "%u", get_bytes( uval, sgn, 0, bpv ) );
	else fprintf( outf, "%i", get_bytes( uval, sgn, 0, bpv ) );
	for ( i = 1; i < (int) order; i++ ) {
		if ( sgn )	fprintf( outf, "/%u", get_bytes( uval, sgn, i, bpv ) );
		else fprintf( outf, "/%i", get_bytes( uval, sgn, i, bpv ) );
	}
}
