#define MAXCHAR 1024

/* min & max macro */
#define MIN(A,B) (A<B?A:B)
#define MAX(A,B) (A>B?A:B)

/* atom type macro */
#define COMP(ityp,jtyp) ((ityp\
<jtyp)?(jtyp+1)*jtyp/2+ityp:(i\
typ+1)*ityp/2+jtyp)

/* memory macros */
#define GROW(A,oldsize,newsize\
) memory->grow(A,oldsize,newsi\
ze,#A,__LINE__,__FILE__,__FUNC\
TION__)

#define CREATE1D(A,d0) memory-\
>create(A,d0,#A,__LINE__,__FIL\
E__,__FUNCTION__)

#define CREATE2D(A,d0,d1) memo\
ry->create(A,d0,d1,#A,__LINE__\
,__FILE__,__FUNCTION__)

/* 3 dimensional shorthand macros */
#define V3ZERO(V) (V[0]=0,V[1]\
=0,V[2]=0)

#define M3ZERO(A) (A[0][0]=0,\
A[0][1]=0,A[0][2]=0,A[1][0]=0\
,A[1][1]=0,A[1][2]=0,A[2][0]=\
0,A[2][1]=0,A[2][2]=0)

#define M3DET(A) A[0][0]*A[1]\
[1]*A[2][2]+A[1][0]*A[2][1]*A\
[0][2]+A[2][0]*A[0][1]*A[1][2\
]-A[0][2]*A[1][1]*A[2][0]-A[1\
][2]*A[2][1]*A[0][0]-A[2][2]*A\
[0][1]*A[1][0]

#define M3EQV(A,B) (B[0][0]=\
A[0][0],B[0][1]=A[0][1],B[0]\
[2]=A[0][2],B[1][0]=A[1][0],\
B[1][1]=A[1][1],B[1][2]=A[1]\
[2],B[2][0]=A[2][0],B[2][1]=\
A[2][1],B[2][2]=A[2][2])

#define M3INV_TRI_LOWER(A,B)\
(B[0][0]=1.0/A[0][0],B[0][1]\
=0,B[0][2]=0,B[1][0]=-A[1][0\
]/(A[0][0]*A[1][1]),B[1][1]=\
1.0/A[1][1],B[1][2]=0,B[2][0\
]=(A[1][0]==0 || A[2][1]==0)\
? (A[2][0]==0 ? 0 : (-A[1][1\
]*A[2][0])/(A[0][0]*A[1][1]*\
A[2][2])) : (A[2][0]==0 ? (A\
[1][0]*A[2][1])/(A[0][0]*A[1\
][1]*A[2][2]) :(A[1][0]*A[2]\
[1]-A[1][1]*A[2][0])/(A[0][0\
]*A[1][1]*A[2][2])),B[2][1]=\
-A[2][1]/(A[1][1]*A[2][2]),B\
[2][2]=1.0/A[2][2])

#define M3MUL_TRI_LOWER(A,B,\
C) (C[0][0]=A[0][0]*B[0][0],\
C[0][1]=0,C[0][2]=0,C[1][0]=\
A[1][0]*B[0][0]+A[1][1]*B[1]\
[0],C[1][1]=A[1][1]*B[1][1],\
C[1][2]=0,C[2][0]=A[2][0]*B[\
0][0]+A[2][1]*B[1][0]+A[2][2\
]*B[2][0],C[2][1]=A[2][1]*B[\
1][1]+A[2][2]*B[2][1],C[2][2\
]=A[2][2]*B[2][2])

