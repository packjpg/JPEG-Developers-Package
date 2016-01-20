#define RBITS( c, n )		( c & ( 0xFF >> (8 - n) ) )
#define LBITS( c, n )		( c >> (8 - n) )
#define MBITS( c, l, r )	( RBITS( c,l ) >> r )
#define RBITS16( c, n )		( c & ( 0xFFFFFFFF >> (16 - n) ) )
#define LBITS16( c, n )		( c >> (16 - n) )
#define MBITS16( c, l, r )	( RBITS16( c,l ) >> r )
#define RBITS32( c, n )		( c & ( 0xFFFFFFFF >> (32 - n) ) )
#define LBITS32( c, n )		( c >> (32 - n) )
#define MBITS32( c, l, r )	( RBITS32( c,l ) >> r )
#define BITN( c, n )		( (c >> n) & 0x1 )
#define BITLEN( l, v )		for ( l = 0; ( v >> l ) > 0; l++ )
#define FDIV2( v, p )		( ( v < 0 ) ? -( (-v) >> p ) : ( v >> p ) )


/* -----------------------------------------------
	class to read arrays bitwise
	----------------------------------------------- */

class abitreader
{
public:
	abitreader( unsigned char* array, int size );
	~abitreader( void );	
	unsigned int read( int nbits );
	unsigned char unpad( unsigned char fillbit );
	int getpos( void );	
	bool eof;
	
private:
	unsigned char* data;
	int lbyte;
	int cbyte;
	int cbit;
};


/* -----------------------------------------------
	class to write arrays bitwise
	----------------------------------------------- */

class abitwriter
{
public:
	abitwriter( int size );
	~abitwriter( void );	
	void write( unsigned int val, int nbits );
	void pad ( unsigned char fillbit );
	unsigned char* getptr( void );
	int getpos( void );
	bool error;	
	unsigned char fillbit;
	
private:
	unsigned char* data;
	int dsize;
	int adds;
	int lbyte;
	int cbyte;
	int cbit;
	bool fmem;
};


/* -----------------------------------------------
	class to write arrays bytewise
	----------------------------------------------- */

class abytewriter
{
public:
	abytewriter( int size );
	~abytewriter( void );	
	void write( unsigned char byte );
	void write_n( unsigned char* byte, int n );
	unsigned char* getptr( void );
	unsigned char* peekptr( void );
	int getpos( void );
	void reset( void );
	bool error;	
	
private:
	unsigned char* data;
	int dsize;
	int adds;
	int lbyte;
	int cbyte;
	bool fmem;
};
