/*
This file contains special classes for bitwise
reading and writing of arrays
*/

#include <stdlib.h>
#include <string.h>
#include "bitops.h"


/* -----------------------------------------------
	constructor for abitreader class
	----------------------------------------------- */	

abitreader::abitreader( unsigned char* array, int size )
{
	cbyte = 0;	
	cbit = 8;
	eof = false;
	
	data = array;
	lbyte = size;	
}

/* -----------------------------------------------
	destructor for abitreader class
	----------------------------------------------- */	

abitreader::~abitreader( void )
{
}

/* -----------------------------------------------
	reads n bits from abitreader
	----------------------------------------------- */	

unsigned int abitreader::read( int nbits )
{
	unsigned int retval = 0;
	
	// safety check for eof
	if (eof) return 0;
	
	while ( nbits >= cbit ) {
		nbits -= cbit;
		retval |= ( RBITS( data[cbyte], cbit ) << nbits );		
		cbit = 8;
		if ( ++cbyte >= lbyte ) {
			eof = true;
			return retval;
		}
	}
	
	if ( nbits > 0 ) {		
		retval |= ( MBITS( data[cbyte], cbit, (cbit-nbits) ) );
		cbit -= nbits;		
	}
	
	return retval;
}

/* -----------------------------------------------
	to skip padding from current byte
	----------------------------------------------- */

unsigned char abitreader::unpad( unsigned char fillbit )
{
	if ( ( cbit == 8 ) || eof ) return fillbit;
	else {
		fillbit = read( 1 );
		while ( cbit != 8 ) read( 1 );
	}
	
	return fillbit;
}

/* -----------------------------------------------
	get current position in array
	----------------------------------------------- */	

int abitreader::getpos( void )
{
	return cbyte;
}


/* -----------------------------------------------
	constructor for abitwriter class
	----------------------------------------------- */	

abitwriter::abitwriter( int size )
{
	fillbit = 1;
	adds    = 65536;
	cbyte   = 0;
	cbit    = 8;
	
	error = false;
	fmem  = true;
	
	dsize = ( size > 0 ) ? size : adds;
	data = ( unsigned char* ) malloc ( dsize );
	if ( data == NULL ) {
		error = true;
		return;
	}
	
	// fill buffer with zeroes
	memset( data, 0, dsize * sizeof( char ) );
	// for ( int i = 0; i < dsize; i++ ) data[i] = 0;
}

/* -----------------------------------------------
	destructor for abitwriter class
	----------------------------------------------- */	

abitwriter::~abitwriter( void )
{
	// free memory if pointer was not given out
	if ( fmem )	free( data );
}

/* -----------------------------------------------
	writes n bits to abitwriter
	----------------------------------------------- */	

void abitwriter::write( unsigned int val, int nbits )
{
	// safety check for error
	if ( error ) return;
	
	// test if pointer beyond flush treshold
	if ( cbyte > ( dsize - 5 ) ) {
		dsize += adds;
		data = (unsigned char*) realloc( data, dsize );
		if ( data == NULL ) {
			error = true;
			return;
		}
		memset( ( data + cbyte + 1 ), 0, ( dsize - ( cbyte + 1 ) ) * sizeof( char ) );
		// for ( int i = cbyte + 1; i < dsize; i++ ) data[i] = 0;
	}
	
	// write data
	while ( nbits >= cbit ) {
		data[cbyte] |= ( MBITS32(val, nbits, (nbits-cbit)) );		
		nbits -= cbit;		
		cbyte++;
		cbit = 8;
	}
	
	if ( nbits > 0 ) {		
		data[cbyte] |= ( (RBITS32(val, nbits)) << (cbit - nbits) );
		cbit -= nbits;		
	}	
}

/* -----------------------------------------------
	pads data using fillbit
	----------------------------------------------- */
	
void abitwriter::pad( unsigned char fillbit )
{
	while ( cbit < 8 )
		write( fillbit, 1 );
}

/* -----------------------------------------------
	gets data array from abitwriter
	----------------------------------------------- */	

unsigned char* abitwriter::getptr( void )
{
	// data is padded here
	pad( fillbit );
	// forbid freeing memory
	fmem = false;
	// realloc data
	data = (unsigned char*) realloc( data, cbyte );
	
	return data;
}

/* -----------------------------------------------
	gets size of data array from abitwriter
	----------------------------------------------- */	

int abitwriter::getpos( void )
{
	return cbyte;
}


/* -----------------------------------------------
	constructor for abytewriter class
	----------------------------------------------- */	

abytewriter::abytewriter( int size )
{
	adds  = 65536;
	cbyte = 0;
	
	error = false;
	fmem  = true;
	
	dsize = ( size > 0 ) ? size : adds;
	data = (unsigned char*) malloc( dsize );
	if ( data == NULL ) {
		error = true;
		return;
	}
}

/* -----------------------------------------------
	destructor for abytewriter class
	----------------------------------------------- */	

abytewriter::~abytewriter( void )
{
	// free data if pointer is not read
	if ( fmem )	free( data );
}

/* -----------------------------------------------
	writes 1 byte to abytewriter
	----------------------------------------------- */	

void abytewriter::write( unsigned char byte )
{
	// safety check for error
	if ( error ) return;
	
	// test if pointer beyond flush threshold
	if ( cbyte >= ( dsize - 2 ) ) {
		dsize += adds;
		data = (unsigned char*) realloc( data, dsize );
		if ( data == NULL ) {
			error = true;
			return;
		}
	}
	
	// write data
	data[ cbyte++ ] = byte;
}

/* -----------------------------------------------
	writes n byte to abytewriter
	----------------------------------------------- */
	
void abytewriter::write_n( unsigned char* byte, int n )
{
	// safety check for error
	if ( error ) return;
	
	// make sure that pointer doesn't get beyond flush threshold
	while ( ( cbyte + n ) >= ( dsize - 2 ) ) {
		dsize += adds;
		data = (unsigned char*) realloc( data, dsize );
		if ( data == NULL ) {
			error = true;
			return;
		}
	}
	
	// copy data from array
	while ( n-- > 0 )
		data[ cbyte++ ] = *(byte++);
}

/* -----------------------------------------------
	gets data array from abytewriter
	----------------------------------------------- */

unsigned char* abytewriter::getptr( void )
{
	// forbid freeing memory
	fmem = false;
	// realloc data
	data = (unsigned char*) realloc( data, cbyte );
	
	return data;
}

/* -----------------------------------------------
	peeks into data array from abytewriter
	----------------------------------------------- */
	
unsigned char* abytewriter::peekptr( void )
{
	return data;
}

/* -----------------------------------------------
	gets size of data array from abytewriter
	----------------------------------------------- */	

int abytewriter::getpos( void )
{
	return cbyte;
}

/* -----------------------------------------------
	reset without realloc
	----------------------------------------------- */	
	
void abytewriter::reset( void )
{
	// set position of current byte
	cbyte = 0;
}

