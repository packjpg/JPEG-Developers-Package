/* lpaq1.cpp file compressor, July 24, 2007.
(C) 2007, Matt Mahoney, matmahoney@yahoo.com

    LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    Visit <http://www.gnu.org/copyleft/gpl.html>.

To compress:   N input output  (N=0..9, uses 3 + 3*2^N MB memory)
To decompress: d input output  (requires same memory)

For example:

  lpaq1 9 foo foo.lpq

compresses foo to foo.lpq using 1.5 GB memory.

  lpaq1 d foo.lpq foo

decompresses to foo, also using 1.5 GB memory.  Option 0 uses
6 MB memory.  The maximum file size is 2 GB.

lpaq1 is a "lite" version of PAQ, about 35 times faster than paq8l
at the cost of some compression (but similar to high end BWT and
PPM compressors).  It is a context-mixing compressor combining 7 
contexts: orders 1, 2, 3, 4, 6, a lowercase unigram word context 
(for ASCII text), and a "match" order, which predicts the next
bit in the last matching context.  The independent bit predictions of
the 7 models are combined using one of 80 neural networks (selected by
a small context), then adjusted using 2 SSE stages (order 0 and 1)
and arithmetic coded.

Prediction is bitwise.  This means that an order-n context consists
of the last n whole bytes plus any of the 0 to 7 previously coded
bits of the current byte starting with the most significant bit.
The unigram word context consists of a hash of the last (at most) 11
consecutive letters (A-Z, a-z) folded to lower case.  The context
does not include any nonalphabetic characters nor any characters
preceding the last nonalphabetic character.

The first 6 contexts (orders 1..4, 6, word) are used to index a
hash table to look up a bit-history represented by an 8-bit state.
The representable states are the same as in paq8l.  A state can
either represent all histories up to 4 bits long, or a pair of
0,1 counts plus a flag to indicate the most recent bit.  The counts
are bounded by (41,0), (40,1), (12,2), (5,3), (4,4) and likewise
for 1,0.  When a count is exceeded, the opposite count is reduced to 
approximately preserve the count ratio.  The last bit flag is present
only for states whose total count is less than 16.  There are 253
possible states.

A bit history is mapped to a probability using an adaptive table
(StateMap).  This differs from paq8l in that each table entry includes
a count so that adaptation is rapid at first.  Each table entry
has a 22-bit probability (initially p = 0.5) and 10-bit count (initially 
n = 0) packed into 32 bits.  After bit y is predicted, n is incremented
up to the limit (1023) and the probability is adjusted by 
p := p + (y - p)/(n + 0.5).  This model is stationary: 
p = (n1 + 0.5)/(n + 1), where n1 is the number of times y = 1 out of n.

The "match" model (MatchModel) looks up the current context in a
hash table, first using a longer context, then a shorter one.  If
a match is found, then the following bits are predicted until there is
a misprediction.  The prediction is computed by mapping the predicted
bit, the length of the match (1..15 or quantized by 4 in 16..62, max 62),
and the last whole byte as context into a StateMap.  If no match is found,
then the order 0 context (last 0..7 bits of the current byte) is used
as context to the StateMap.

The 7 predictions are combined using a neural network (Mixer) as in
paq8l, except it is a single level network without MMX code.  The
inputs p_i, i=0..6 are first stretched: t_i = log(p_i/(1 - p_i)), 
then the output is computed: p = squash(SUM_i t_i * w_i), where
squash(x) = 1/(1 + exp(-x)) is the inverse of stretch().  The weights
are adjusted to reduce the error: w_i := w_i + L * t_i * (y - p) where
(y - p) is the prediction error and L ~ 0.002 is the learning rate.
This is a standard single layer backpropagation network modified to
minimize coding cost rather than RMS prediction error (thus dropping
the factors p * (1 - p) from learning).

One of 80 neural networks are selected by a context that depends on
the 3 high order bits of the last whole byte plus the context order
(quantized to 0, 1, 2, 3, 4, 6, 8, 12, 16, 32).  The order is
determined by the number of nonzero bit histories and the length of
the match from MatchModel.

The Mixer output is adjusted by 2 SSE stages (called APM for adaptive
probability map).  An APM is a StateMap that accepts both a discrte
context and an input probability, pr.  pr is stetched and quantized
to 24 levels.  The output is interpolated between the 2 nearest
table entries, and then only the nearest entry is updated.  The entries
are initialized to p = pr and n = 6 (to slow down initial adaptation)
with a limit n <= 255.  The APM differs from paq8l in that it uses
the new StateMap rapid initial adaptation, does not update both
adjacent table entries, and uses 24 levels instead of 33.  The two
stages use a discrete order 0 context (last 0..7 bits) and a hashed order-1
context (14 bits).  Each output is averaged with its input weighted
by 1/4.

The output is arithmetic coded.  The code for a string s with probability
p(s) is a number between Q and Q+p(x) where Q is the total probability
of all strings lexicographically preceding s.  The number is coded as
a big-endian base-256 fraction.  A header is prepended as follows:

- "pQ" 2 bytes must be present or decompression gives an error.
- 1 (0x01) version number (other values give an error).
- memory option N as one byte '0'..'9' (0x30..0x39).
- file size as a 4 byte big-endian number.
- arithmetic coded data.

Two thirds of the memory (2 * 2^N MB) is used for a hash table mapping
the 6 regular contexts (orders 1-4, 6, word) to bit histories.  A lookup
occurs every 4 bits.  The input is a byte-oriented context plus possibly
the first nibble of the next byte.  The output is an array of 15 bit
histories (1 byte each) for all possible contexts formed by appending
0..3 more bits.  The table entries have the format:

 {checksum, "", 0, 1, 00, 10, 01, 11, 000, 100, 010, 110, 001, 101, 011, 111}

The second byte is the bit history for the context ending on a nibble
boundary.  It also serves as a priority for replacement.  The states
are ordered by increasing total count, where state 0 represents the
initial state (no history).  When a context is looked up, the 8 bit
checksum (part of the hash) is compared with 3 adjacent entries, and
if there is no match, the entry with the lowest priority is cleared
and the new checksum is stored.

The hash table is aligned on 64 byte cache lines.  A table lookup never
crosses a 64 byte address boundary.  Given a 32-bit hash h of the context,
8 bits are used for the checksum and 17 + N bits are used for the
index i.  Then the entries i, i XOR 1, and i XOR 2 are tried.  The hash h
is actually a collision-free permuation, consisting of multiplying the
input by a large odd number mod 2^32, a 16-bit rotate, and another multiply.

The order-1 context is mapped to a bit history using a 64K direct
lookup table, not a hash table.

One third of memory is used by MatchModel, divided equally between 
a rotating input buffer of 2^(N+19) bytes and an index (hash table)
with 2^(N+17) entries.  Two context hashes are maintained, a long one,
h1, of length ceil((N+17)/3) bytes and a shorter one, h2, of length 
ceil((N+17)/5) bytes, where ceil() is the ceiling function.  The index
does not use collision detection.  At each byte boundary, if there is 
not currently a match, then the bytes before the current byte are
compared with the location indexed by h1.  If less than 2 bytes match,
then h2 is tried.  If a match of length 1 or more is found, the
match is maintained until the next bit mismatches the predicted bit.
The table is updated at h1 and h2 after every byte.

To compile (g++ 3.4.5, upx 3.00w):
  g++ -Wall lpaq1.cpp -O2 -Os -march=pentiumpro -fomit-frame-pointer 
      -s -o lpaq1.exe
  upx -qqq lpaq1.exe

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#define NDEBUG  // remove for debugging
#include <assert.h>
#include "lpaqlib.h"

// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

// Error handler: print message if any, and exit
void quit(const char* message=0) {
  if (message) printf("%s\n", message);
  exit(1);
}

// Create an array p of n elements of type T
template <class T> void alloc(T*&p, int n) {
  p=(T*)calloc(n, sizeof(T));
  if (!p) quit("out of memory");
}

///////////////////////////// Squash //////////////////////////////

// return p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
int squash(int d) {
  static const int t[33]={
    1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
    1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
    4050,4068,4079,4085,4089,4092,4093,4094};
  if (d>2047) return 4095;
  if (d<-2047) return 0;
  int w=d&127;
  d=(d>>7)+16;
  return (t[d]*(128-w)+t[(d+1)]*w+64) >> 7;
}

//////////////////////////// Stretch ///////////////////////////////

// Inverse of squash. stretch(d) returns ln(p/(1-p)), d scaled by 8 bits,
// p by 12 bits.  d has range -2047 to 2047 representing -8 to 8.  
// p has range 0 to 4095 representing 0 to 1.

  static int stretch_t[4096];
  static int dt[1024];  // i -> 16K/(i+3)

///////////////////////// state table ////////////////////////

// State table:
//   nex(state, 0) = next state if bit y is 0, 0 <= state < 256
//   nex(state, 1) = next state if bit y is 1
//
// States represent a bit history within some context.
// State 0 is the starting state (no bits seen).
// States 1-30 represent all possible sequences of 1-4 bits.
// States 31-252 represent a pair of counts, (n0,n1), the number
//   of 0 and 1 bits respectively.  If n0+n1 < 16 then there are
//   two states for each pair, depending on if a 0 or 1 was the last
//   bit seen.
// If n0 and n1 are too large, then there is no state to represent this
// pair, so another state with about the same ratio of n0/n1 is substituted.
// Also, when a bit is observed and the count of the opposite bit is large,
// then part of this count is discarded to favor newer data over old.

static const U8 State_table[256][2]={
{  1,  2},{  3,  5},{  4,  6},{  7, 10},{  8, 12},{  9, 13},{ 11, 14}, // 0
{ 15, 19},{ 16, 23},{ 17, 24},{ 18, 25},{ 20, 27},{ 21, 28},{ 22, 29}, // 7
{ 26, 30},{ 31, 33},{ 32, 35},{ 32, 35},{ 32, 35},{ 32, 35},{ 34, 37}, // 14
{ 34, 37},{ 34, 37},{ 34, 37},{ 34, 37},{ 34, 37},{ 36, 39},{ 36, 39}, // 21
{ 36, 39},{ 36, 39},{ 38, 40},{ 41, 43},{ 42, 45},{ 42, 45},{ 44, 47}, // 28
{ 44, 47},{ 46, 49},{ 46, 49},{ 48, 51},{ 48, 51},{ 50, 52},{ 53, 43}, // 35
{ 54, 57},{ 54, 57},{ 56, 59},{ 56, 59},{ 58, 61},{ 58, 61},{ 60, 63}, // 42
{ 60, 63},{ 62, 65},{ 62, 65},{ 50, 66},{ 67, 55},{ 68, 57},{ 68, 57}, // 49
{ 70, 73},{ 70, 73},{ 72, 75},{ 72, 75},{ 74, 77},{ 74, 77},{ 76, 79}, // 56
{ 76, 79},{ 62, 81},{ 62, 81},{ 64, 82},{ 83, 69},{ 84, 71},{ 84, 71}, // 63
{ 86, 73},{ 86, 73},{ 44, 59},{ 44, 59},{ 58, 61},{ 58, 61},{ 60, 49}, // 70
{ 60, 49},{ 76, 89},{ 76, 89},{ 78, 91},{ 78, 91},{ 80, 92},{ 93, 69}, // 77
{ 94, 87},{ 94, 87},{ 96, 45},{ 96, 45},{ 48, 99},{ 48, 99},{ 88,101}, // 84
{ 88,101},{ 80,102},{103, 69},{104, 87},{104, 87},{106, 57},{106, 57}, // 91
{ 62,109},{ 62,109},{ 88,111},{ 88,111},{ 80,112},{113, 85},{114, 87}, // 98
{114, 87},{116, 57},{116, 57},{ 62,119},{ 62,119},{ 88,121},{ 88,121}, // 105
{ 90,122},{123, 85},{124, 97},{124, 97},{126, 57},{126, 57},{ 62,129}, // 112
{ 62,129},{ 98,131},{ 98,131},{ 90,132},{133, 85},{134, 97},{134, 97}, // 119
{136, 57},{136, 57},{ 62,139},{ 62,139},{ 98,141},{ 98,141},{ 90,142}, // 126
{143, 95},{144, 97},{144, 97},{ 68, 57},{ 68, 57},{ 62, 81},{ 62, 81}, // 133
{ 98,147},{ 98,147},{100,148},{149, 95},{150,107},{150,107},{108,151}, // 140
{108,151},{100,152},{153, 95},{154,107},{108,155},{100,156},{157, 95}, // 147
{158,107},{108,159},{100,160},{161,105},{162,107},{108,163},{110,164}, // 154
{165,105},{166,117},{118,167},{110,168},{169,105},{170,117},{118,171}, // 161
{110,172},{173,105},{174,117},{118,175},{110,176},{177,105},{178,117}, // 168
{118,179},{110,180},{181,115},{182,117},{118,183},{120,184},{185,115}, // 175
{186,127},{128,187},{120,188},{189,115},{190,127},{128,191},{120,192}, // 182
{193,115},{194,127},{128,195},{120,196},{197,115},{198,127},{128,199}, // 189
{120,200},{201,115},{202,127},{128,203},{120,204},{205,115},{206,127}, // 196
{128,207},{120,208},{209,125},{210,127},{128,211},{130,212},{213,125}, // 203
{214,137},{138,215},{130,216},{217,125},{218,137},{138,219},{130,220}, // 210
{221,125},{222,137},{138,223},{130,224},{225,125},{226,137},{138,227}, // 217
{130,228},{229,125},{230,137},{138,231},{130,232},{233,125},{234,137}, // 224
{138,235},{130,236},{237,125},{238,137},{138,239},{130,240},{241,125}, // 231
{242,137},{138,243},{130,244},{245,135},{246,137},{138,247},{140,248}, // 238
{249,135},{250, 69},{ 80,251},{140,252},{249,135},{250, 69},{ 80,251}, // 245
{140,252},{  0,  0},{  0,  0},{  0,  0}};  // 252
#define nex(state,sel) State_table[state][sel]

//////////////////////////// StateMap, APM //////////////////////////

// A StateMap maps a context to a probability.  Methods:
//
// Statemap sm(n) creates a StateMap with n contexts using 4*n bytes memory.
// sm.p(y, cx, limit) converts state cx (0..n-1) to a probability (0..4095).
//     that the next y=1, updating the previous prediction with y (0..1).
//     limit (1..1023, default 1023) is the maximum count for computing a
//     prediction.  Larger values are better for stationary sources.

class StateMap {
protected:
  const int N;  // Number of contexts
  int cxt;      // Context of last prediction
  U32 *t;       // cxt -> prediction in high 22 bits, count in low 10 bits
  inline void update(int y, int limit) {
    assert(cxt>=0 && cxt<N);
    U32 *p=&t[cxt], p0=p[0];
    int n=p0&1023, pr=p0>>10;  // count, prediction
    if (n<limit) ++p0;
    //else p0=pr*1024+limit;
    p0+=(y-pr>>3)*dt[n]&0xfffffc00;
    p[0]=p0;
  }
public:
  StateMap(int n=256);

  // update bit y (0..1), predict next bit in context cx
  int p(int y, int cx, int limit=1023) {
    assert(y>>1==0);
    assert(cx>=0 && cx<N);
    assert(limit>0 && limit<1024);
    update(y, limit);
    return t[cxt=cx]>>20;
  }
};

StateMap::StateMap(int n): N(n), cxt(0) {
  alloc(t, N);
  for (int i=0; i<N; ++i)
    t[i]=1<<31;
}

// An APM maps a probability and a context to a new probability.  Methods:
//
// APM a(n) creates with n contexts using 96*n bytes memory.
// a.pp(y, pr, cx, limit) updates and returns a new probability (0..4095)
//     like with StateMap.  pr (0..4095) is considered part of the context.
//     The output is computed by interpolating pr into 24 ranges nonlinearly
//     with smaller ranges near the ends.  The initial output is pr.
//     y=(0..1) is the last bit.  cx=(0..n-1) is the other context.
//     limit=(0..1023) defaults to 255.

class APM: public StateMap {
public:
  APM(int n);
  int pp(int y, int pr, int cx, int limit=255) {
    assert(y>>1==0);
    assert(pr>=0 && pr<4096);
    assert(cx>=0 && cx<N/24);
    assert(limit>0 && limit<1024);
    update(y, limit);
    pr=(stretch_t[pr]+2048)*23;
    int wt=pr&0xfff;  // interpolation weight of next element
    cx=cx*24+(pr>>12);
    assert(cx>=0 && cx<N-1);
    cxt=cx+(wt>>11);
    pr=(t[cx]>>13)*(0x1000-wt)+(t[cx+1]>>13)*wt>>19;
    return pr;
  }
};

APM::APM(int n): StateMap(n*24) {
  for (int i=0; i<N; ++i) {
    int p=((i%24*2+1)*4096)/48-2048;
    t[i]=(U32(squash(p))<<20)+6;
  }
}

//////////////////////////// Mixer /////////////////////////////

// Mixer m(N, M) combines models using M neural networks with
//     N inputs each using 4*M*N bytes of memory.  It is used as follows:
// m.update(y) trains the network where the expected output is the
//     last bit, y.
// m.add(stretch(p)) inputs prediction from one of N models.  The
//     prediction should be positive to predict a 1 bit, negative for 0,
//     nominally -2K to 2K.
// m.set(cxt) selects cxt (0..M-1) as one of M neural networks to use.
// m.p() returns the output prediction that the next bit is 1 as a
//     12 bit number (0 to 4095).  The normal sequence per prediction is:
//
// - m.add(x) called N times with input x=(-2047..2047)
// - m.set(cxt) called once with cxt=(0..M-1)
// - m.p() called once to predict the next bit, returns 0..4095
// - m.update(y) called once for actual bit y=(0..1).

  int* mxr_tx;         // N inputs

#if 0			// ATTENTION !  CHANGE this to 1 if you start to use
                        //              <mixer max inputs>!=7  in your version.

inline void train(int *w, int err, int n) {
  int *t = mxr_tx;
  for (int i=0; i<n; ++i)
    w[i]+=t[i]*err+0x8000>>16;
}

inline int dot_product(int *w, int n) {
  int *t = mxr_tx;
  int sum=0;
  for (int i=0; i<n; ++i)
    sum+=t[i]*w[i];
  return sum;
}

#else
inline void train(int *w, int err) {
    int *t = mxr_tx;
    w[0]+=t[0]*err+0x8000>>16;
    w[1]+=t[1]*err+0x8000>>16;
    w[2]+=t[2]*err+0x8000>>16;
    w[3]+=t[3]*err+0x8000>>16;
    w[4]+=t[4]*err+0x8000>>16;
    w[5]+=t[5]*err+0x8000>>16;
    w[6]+=t[6]*err+0x8000>>16;
}
inline int dot_product(int *w) {
    int *t = mxr_tx;
    int sum =t[0]*w[0];
	sum+=t[1]*w[1];
	sum+=t[2]*w[2];
	sum+=t[3]*w[3];
	sum+=t[4]*w[4];
	sum+=t[5]*w[5];
	sum+=t[6]*w[6];
  return sum;
}
#endif

class Mixer {
//const int MI;	// max inputs
#define MI 7
  const int MC;	// max contexts
  int* wx;         // MI*M weights
  int cxt;         // context
//int nx;          // Number of inputs in tx, 0 to MI
  int pr;          // last result (scaled 12 bits)
public:
  Mixer(int m);

  // Adjust weights to minimize coding cost of last prediction
  void update(int y) {
    int err=((y<<12)-pr)*7;
    assert(err>=-32768 && err<32768);
    train(&wx[cxt*MI], err);
  }

  // Input x (call up to MI times)
#if 0
  void add(int nx, int x) {
    assert(nx<MI);
    mxr_tx[nx]=stretch_t[x);
  }
#define m_add(a,b) m.add(a,b)
#else
#define m_add(a,b) mxr_tx[a]=stretch_t[b];
#endif

  // Set a context
  void set(int cx) {
    assert(cx>=0 && cx<MC);
    cxt=cx;
  }

  // predict next bit
  int p() {
    return pr=squash(dot_product(&wx[cxt*MI])>>16);
  }
};

Mixer::Mixer(int m):	MC(m), cxt(0), pr(2048) {
//MI=n;
  assert(MC>0);
  alloc(mxr_tx, MI);
  alloc(wx, MI*MC);
}

//////////////////////////// HashTable /////////////////////////

// A HashTable maps a 32-bit index to an array of B bytes.
// The first byte is a checksum using the upper 8 bits of the
// index.  The second byte is a priority (0 = empty) for hash
// replacement.  The index need not be a hash.

// HashTable<B> h(n) - create using n bytes  n and B must be 
//     powers of 2 with n >= B*4, and B >= 2.
// h[i] returns array [1..B-1] of bytes indexed by i, creating and
//     replacing another element if needed.  Element 0 is the
//     checksum and should not be modified.

template <int B>
class HashTable {
  U8* t;  // table: 1 element = B bytes: checksum priority data data
  const int N;  // size in bytes
public:
  HashTable(int n);
  U8* operator[](U32 i);
};

template <int B>
HashTable<B>::HashTable(int n): t(0), N(n) {
  assert(B>=2 && (B&B-1)==0);
  assert(N>=B*4 && (N&N-1)==0);
  alloc(t, N+B*4+64);
  t+=64-int(((long)t)&63);  // align on cache line boundary
}

template <int B>
inline U8* HashTable<B>::operator[](U32 i) {
  i*=123456791;
  i=i<<16|i>>16;
  i*=234567891;
  int chk=i>>24;
  i=i*B&N-B;
  if (t[i]==chk) return t+i;
  if (t[i^B]==chk) return t+(i^B);
  if (t[i^B*2]==chk) return t+(i^B*2);
  if (t[i+1]>t[i+1^B] || t[i+1]>t[i+1^B*2]) i^=B;
  if (t[i+1]>t[i+1^B^B*2]) i^=B^B*2;

#if 0			// ATTENTION !  CHANGE this to 1 if you start to use
                        //              HashTable with B!=16 in your version.
  memset(t+i, 0, B);
  t[i]=chk;
  return t+i;
#else
	{
		U8* p=t+i;
		*(U32*)p = chk;
		*(U32*)(p+ 4)=0;
		*(U32*)(p+ 8)=0;
		*(U32*)(p+12)=0;
		return p;
	}
#endif
}

//////////////////////////// MatchModel ////////////////////////

// MatchModel(n) predicts next bit using most recent context match.
//     using n bytes of memory.  n must be a power of 2 at least 8.
// MatchModel::p(y, m) updates the model with bit y (0..1) and writes
//     a prediction of the next bit to Mixer m.  It returns the length of
//     context matched (0..62).

  static int c0=1;  // last 0-7 bits with leading 1
  U8* buf;    // input buffer
  U32 h1, h2; // context hashes
  int pos;    // number of bytes in buf
  int N, HN; // last hash table index, n/8-1

class MatchModel {
  enum {MAXLEN=62};   // maximum match length, at most 62
  int* ht;    // context hash -> next byte in buf
  int match;  // pointer to current byte in matched context in buf
  int len;    // length of match
  int bcount; // number of bits in c0 (0..7)
  StateMap sm;  // len, bit, last byte -> prediction
public:
  MatchModel(int n);  // n must be a power of 2 at least 8.
  int p(int y, Mixer& m);  // update bit y (0..1), predict next bit to m
};

MatchModel::MatchModel(int n): ht(0), match(0), len(0), bcount(0), sm(56<<8) {
  N=n/2-1;
  HN=n/8-1;
  pos=h1=h2=0;
  assert(n>=8 && (n&n-1)==0);
  alloc(buf, N+1);
  alloc(ht, HN+1);
}

int MatchModel::p(int y, Mixer& m) {

  // update context
  ++bcount;
  if (bcount==8) {
    bcount=0;

    // find or extend match
    if (len>0) {
      ++match;
      match&=N;
      if (len<MAXLEN) ++len;
    }
    else {
      match=ht[h1];
      if (match!=pos) {
        int i;
        while (len<MAXLEN && (i=match-len-1&N)!=pos
               && buf[i]==buf[pos-len-1&N])
          ++len;
      }
    }
    if (len<2) {
      len=0;
      match=ht[h2];
      if (match!=pos) {
        int i;
        while (len<MAXLEN && (i=match-len-1&N)!=pos
               && buf[i]==buf[pos-len-1&N])
          ++len;
      }
    }
    // update index
    ht[h1]=pos;
    ht[h2]=pos;
  }

  // predict
  int cxt=c0;
  if (len>0 && (buf[match]+256>>8-bcount)==cxt) {
    int b=buf[match]>>7-bcount&1;  // next bit
    if (len<16) cxt=len*2+b;
    else cxt=(len>>2)*2+b+24;
    cxt=cxt*256+buf[pos-1&N];
  }
  else
    len=0;

  m_add(0, sm.p(y, cxt));
  return len;
}

//////////////////////////// Predictor /////////////////////////

// A Predictor estimates the probability that the next bit of
// uncompressed data is 1.  Methods:
// Predictor(n) creates with 3*n bytes of memory.
// p() returns P(1) as a 12 bit number (0-4095).
// update(y) trains the predictor with the actual bit (0 or 1).

int MEM=0;  // Global memory usage = 3*MEM bytes (1<<20 .. 1<<29)

class Predictor {
  int pr;  // next prediction
public:
  Predictor();
  int p() const {assert(pr>=0 && pr<4096); return pr;}
  void update(int y);
};

Predictor::Predictor(): pr(2048) {}

void Predictor::update(int y) {
  static U8 t0[0x10000];  // order 1 cxt -> state
  static HashTable<16> t(MEM*2);  // cxt -> state
  static int c4=0;  // last 4 bytes
  static U8 *cp[6]={t0, t0, t0, t0, t0, t0};  // pointer to bit history
  static int bcount=0;  // bit count
  static StateMap sm[6];
  static APM a1(0x100), a2(0x4000);
  static U32 h[6];
  static Mixer m(80);
  static MatchModel mm(MEM);  // predicts next bit by matching context
  assert(MEM>0);
  // update model
  assert(y==0 || y==1);
  *cp[0]=nex(*cp[0], y);
  *cp[1]=nex(*cp[1], y);
  *cp[2]=nex(*cp[2], y);
  *cp[3]=nex(*cp[3], y);
  *cp[4]=nex(*cp[4], y);
  *cp[5]=nex(*cp[5], y);
  m.update(y);

  // update context
  c0+=c0+y;
  if (c0>=256) {
    h1=h1*(3<<3)+c0&HN;
    h2=h2*(5<<5)+c0&HN;
    buf[pos++]=c0;
    pos&=N;

    c0-=256;
    c4=c4<<8|c0;
    h[0]=c0<<8;  // order 1
    h[1]=(c4&0xffff)<<5|0x57000000;  // order 2
    h[2]=(c4<<8)*3;  // order 3
    h[3]=c4*5;  // order 4
    h[4]=h[4]*(11<<5)+c0*13&0x3fffffff;  // order 6
    if (c0>=65 && c0<=90) c0+=32;  // lowercase unigram word order
    if (c0>=97 && c0<=122) h[5]=(h[5]+c0)*(7<<3);
    else h[5]=0;
    cp[1]=t[h[1]]+1;
    cp[2]=t[h[2]]+1;
    cp[3]=t[h[3]]+1;
    cp[4]=t[h[4]]+1;
    cp[5]=t[h[5]]+1;
    c0=1;
    bcount=0;
  }
  else {
        ++bcount;
        if (bcount==4) {
          cp[1]=t[h[1]+c0]+1;
          cp[2]=t[h[2]+c0]+1;
          cp[3]=t[h[3]+c0]+1;
          cp[4]=t[h[4]+c0]+1;
          cp[5]=t[h[5]+c0]+1;
        }
        else {
          int j=y+1<<(bcount&3)-1;
          cp[1]+=j;
          cp[2]+=j;
          cp[3]+=j;
          cp[4]+=j;
          cp[5]+=j;
        }
  }
  cp[0]=t0+h[0]+c0;

  // predict
  int y22=y<<22, len=mm.p(y22, m), order;
  if (len==0)
	order=(*cp[1]!=0)+(*cp[2]!=0)+(*cp[3]!=0)+(*cp[4]!=0);
   else order=5+(len>=8)+(len>=12)+(len>=16)+(len>=32);
  m_add(1, sm[0].p(y22, *cp[0]));
  m_add(2, sm[1].p(y22, *cp[1]));
  m_add(3, sm[2].p(y22, *cp[2]));
  m_add(4, sm[3].p(y22, *cp[3]));
  m_add(5, sm[4].p(y22, *cp[4]));
  m_add(6, sm[5].p(y22, *cp[5]));
  m.set(order+10*(h[0]>>13));
  pr=m.p();
  pr=pr+3*a1.pp(y22, pr, c0)>>2;
  pr=pr+3*a2.pp(y22, pr, c0^h[0]>>2)>>2;
}

//////////////////////////// Encoder ////////////////////////////

// An Encoder does arithmetic encoding.  Methods:
// Encoder(COMPRESS, f) creates encoder for compression to archive f, which
//     must be open past any header for writing in binary mode.
// Encoder(DECOMPRESS, f) creates encoder for decompression from archive f,
//     which must be open past any header for reading in binary mode.
// code(i) in COMPRESS mode compresses bit i (0 or 1) to file f.
// code() in DECOMPRESS mode returns the next decompressed bit from file f.
// compress(c) in COMPRESS mode compresses one byte.
// decompress() in DECOMPRESS mode decompresses and returns one byte.
// flush() should be called exactly once after compression is done and
//     before closing f.  It does nothing in DECOMPRESS mode.

typedef enum {COMPRESS, DECOMPRESS} Mode;
class Encoder {
private:
  Predictor predictor;
  const Mode mode;       // Compress or decompress?
  iostream* str_arch;    // Compressed data file
  U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
  U32 x;                 // Decompress mode: last 4 input bytes of archive

  // Compress bit y or return decompressed bit
  int code(int y=0) {
    int p=predictor.p();
	int c=0;
    assert(p>=0 && p<4096);
    p+=p<2048;
    U32 xmid=x1 + (x2-x1>>12)*p + ((x2-x1&0xfff)*p>>12);
    assert(xmid>=x1 && xmid<x2);
    if (mode==DECOMPRESS) y=x<=xmid;
    y ? (x2=xmid) : (x1=xmid+1);
    predictor.update(y);
    while (((x1^x2)&0xff000000)==0) {  // pass equal leading bytes of range
      if (mode==COMPRESS) {
	    c=x2>>24;
		str_arch->write(&c,1,1);
	  }
      x1<<=8;
      x2=(x2<<8)+255;
	  if (mode==DECOMPRESS) {
	    str_arch->read(&c,1,1);
	    if (str_arch->chkeof()) c=EOF;
	    x=(x<<8)+(c&255);		
	  } // EOF is OK
    }
    return y;
  }

public:
  Encoder(Mode m, iostream* f);
  void flush();  // call this when compression is finished

  // Compress one byte
  void compress(int c) {
    assert(mode==COMPRESS);
    for (int i=7; i>=0; --i)
      code((c>>i)&1);
  }

  // Decompress and return one byte
  int decompress() {
    int c=0;
    for (int i=0; i<8; ++i)
      c+=c+code();
    return c;
  }
};

Encoder::Encoder(Mode m, iostream* f):
    mode(m), str_arch(f), x1(0), x2(0xffffffff), x(0) {
  if (mode==DECOMPRESS) {  // x = first 4 bytes of archive
    int c=0;
    for (int i=0; i<4; ++i) {
	  str_arch->read(&c,1,1);
      x=(x<<8)+(c&255);
	}
  }
  int pi=0;
  for (int x=-2047; x<=2047; ++x) {  // invert squash()
    int i=squash(x);
    for (int j=pi; j<=i; ++j)
      stretch_t[j]=x;
    pi=i+1;
  }
  stretch_t[4095]=2047;

  for (int i=0; i<1024; ++i)
    dt[i]=16384/(i+i+3);
}

void Encoder::flush() {
  if (mode==COMPRESS) {
    int c=(x1>>24)&255;  // Flush first unequal byte of range
	str_arch->write(&c,1,1);
  }
}

void lpaq1(iostream* str_in, iostream* str_out, Mode mode) {
  MEM=1<<(5+20); // set memory usage (same as level 5 in original LPAQ1)
  // c0=1;  // last 0-7 bits with leading 1
  buf=NULL;    // input buffer
  // h1=h2=pos=0; // context hashes / number of bytes in buf
  // N=HN=0;
  if (mode==COMPRESS) {
    long size = str_in->getsize();
	if (size<0 || size>=0x7FFFFFFF) quit("input file too big");
	str_out->write(&size,4,1);
	
	// Compress	
    Encoder e(COMPRESS, str_out);
    int c=0;
	while (str_in->read(&c,1,1) > 0)
      e.compress(c);
    e.flush();
  }
  else {
    long size = 0;
	str_in->read(&size,4,1);
	if (size<0) quit("Bad file size");
	
	// Decompress
	Encoder e(DECOMPRESS, str_in);
	int c=0;
    while (size-->0) {
	  c=e.decompress();
	  str_out->write(&c,1,1);
	}
  }
}

//////////////////////////// User Interface ////////////////////////////
/*
// main test function
int main(int argc, char** argv) {
	iostream* str_in;
	iostream* str_out;
	
	// check number of arguments
	if ( argc != 4 ) { printf( "usage: %s [mode c/d] [input file] [output file]\n", argv[ 0 ] ); return 0; }
	
	// init & check streams
	str_in = new iostream( (void*) argv[2], 0, 0, 0 );
	str_out = new iostream( (void*) argv[3], 0, 0, 1 );
	if ( str_in->chkerr() ) { printf( "error opening input \"%s\"\n", argv[2] ); return 0; }
	if ( str_out->chkerr() ) { printf( "error writing output \"%s\"\n", argv[3] ); return 0; }
	
	// (de)compress
	if ( argv[ 1 ][ 0 ] == 'd' ) {
		printf( "decompressing \"%s\" to \"%s\"", argv[2], argv[3] );
		lpaq1( str_in, str_out, DECOMPRESS );
	}
	else {
		printf( "compressing \"%s\" to \"%s\"", argv[2], argv[3] );
		lpaq1( str_in, str_out, COMPRESS );
	}
	
	printf( "\ndone!\nexiting...\n" );
	
	delete( str_in );
	delete( str_out );
	
	return 0;
}
*/

// compress function (wrapper around lpaq1())
void lpaq1_compress( iostream* str_in, iostream* str_out ) {
	lpaq1( str_in, str_out, COMPRESS );
}

// decompress function (wrapper around lpaq1())
void lpaq1_decompress( iostream* str_in, iostream* str_out ) {
	lpaq1( str_in, str_out, DECOMPRESS );
}


